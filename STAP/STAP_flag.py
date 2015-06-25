#! /usr/bin/env python

import os
import re
import ephem
import pyfits
import numpy as np
import shlex
import subprocess
import argparse
import cPickle as pickle
from milk.supervised.randomforest import rf_model
from Utils.Constants import FilenamesSub, sextractor_output
from Utils.RealBogus import RealBogus
from Utils.Catalog import SExtractorDetection as SExDet, match_coords
from Utils.DetectionMerger import DetectionMerger
from Utils.TrackableException import TrackableException as STAP_Error
from LocalSettings.Skymapper.Skymapper import load_rbf_v084 as load_classifier
from STAP_comm import print_cmd_line
from STAP.STAP_tools.headerprop import headerprop


def classify(newname, refname, subname,
             interact=False, noapply=False, timeout=None):
    """Real/Bogus classification of transient detections in a subtracted image
    
    Examine the list of detected sources and apply a machine classifier
    (could be anything -- neural net, boosted decision tree, etc.) to select
    sources more likely to be real astrophysical transients than artifacts.
    """
    print_cmd_line("STAP_flag.py", newname, refname, subname,
                   interact=False, noapply=noapply, timeout=timeout)
    return classify_randomforest(newname, refname, subname, interact=interact,
                                 noapply=noapply, timeout=timeout)

def classify_randomforest(newname, refname, subname,
                          interact=False, noapply=False, timeout=None):
    """Implements a random forest classifier for classify()"""

    # ------------------------------------------------------------------------
    # Global image quantities (get these from user, headers, or calculate)
    # ------------------------------------------------------------------------
    
    # Pickled machine learning classifier; currently we're using milk
    # RS 2015/03/14:  Update -- use floating point random forest classifier
    rbclass = load_classifier()

    # RS 2014/02/19:  Before ANYTHING else is done, propagate the runtag and
    # field ID into the subtraction's SExtractor output FITS header.  We'll
    # need it if we're going to use SExtractor output to make light curves
    # using the DetectionMerger apparatus, and if this stage crashes then
    # the files will be screwed up for later.  We may still want to know what
    # was found on a subtraction even if it crashed.
    sub = FilenamesSub(new=subname)
    substarsname = sextractor_output(subname)
    with pyfits.open(sextractor_output(subname), mode='update') as substarptr:
        substarhdr = substarptr[0].header
        substarhdr.update('RUNTAG', sub.runtag, 'Pipeline run tag')

    # RS 2014/02/19:  We also need to propagate all these keywords which
    # would normally be found in a .cands file by the end of this routine.
    headerprop(newname, substarsname,
               ['ZPCALSRC', 'ZPCALTYP', 'ZPMAG',    'ZPMAGERR',
                'FAINTLIM', 'FLXLIM50', 'FLXLIM95', 'MAGLIM50', 'MAGLIM95'])

    # Ok.  Things we actually *use* from the above:
    # -- FITS pixel data for SUB only
    # -- FITS headers for NEW, REF, SUB
    # -- SExtractor detections from (NEW, REF, SUB)
    newstars, newstarhdr = SExDet.read_fits_file(sextractor_output(newname))
    refstars, refstarhdr = SExDet.read_fits_file(sextractor_output(refname))
    substars, substarhdr = SExDet.read_fits_file(sextractor_output(subname))
    subfitsptr = pyfits.open(sub.sub_fits_name, mode='update')
    subfitshdr, subfitspix = subfitsptr[0].header, subfitsptr[0].data
    reffitsptr = pyfits.open(refname, mode='update')
    reffitshdr, reffitspix = reffitsptr[0].header, reffitsptr[0].data

    # RS 2012/06/25:  Even in cases where the subtraction died, we'll want
    # to make a DS9 region file so it's easy for the user to see what was
    # detected on the SUB.  If we pass the "flag" stage, we'll overwrite
    # this later with new regions color-coded Real or Bogus.
    subds9  = sub.sub_reg_name
    with open(subds9,"w") as rf:
        rf.write("global color=yellow width=2\n")
        rf.write("image\n")
        for s in substars:
            rf.write("box({0:7.2f},{1:7.2f},30,30,0)\n"
                     .format(float(s.x),float(s.y)))

    # Check for obvious bad subtractions
    if len(substars) == 0:
        raise STAP_Error("SExtractor detected no objects on subtraction")
    elif len(substars) > 500:
        raise STAP_Error("SExtractor detected too many objects on subtraction")

    # Compute a few quantities we'll need later
    newfwhm, newbkgsig = newstarhdr["SEEING"], newstarhdr["BKGSIG"]
    reffwhm, refbkgsig = refstarhdr["SEEING"], refstarhdr["BKGSIG"]
    subfwhm, subbkgsig = newstarhdr["SEEING"], substarhdr["BKGSIG"]
    # RS 2012/06/29:  If this is meant to be the limiting flux, it's wrong.
    # Unfortunately we've trained the classifier now this way so it'll have
    # to wait until the next iteration to be fixed.
    newflim = np.pi*(newfwhm**2/4)*subbkgsig
    refflim = np.pi*(reffwhm**2/4)*subbkgsig
    subflim = np.pi*(subfwhm**2/4)*subbkgsig

    # ------------------------------------------------------------------------
    #                           The main event
    # ------------------------------------------------------------------------

    # Initialize a list of candidates.
    candlist = [ ]
    
    # Take a few median quantities for sources found in the REF and NEW images
    medaref  = np.median([r.a     for r in refstars])
    mederef  = np.median([r.e     for r in refstars])
    medthref = np.median([r.theta for r in refstars])
    medanew  = np.median([r.a     for r in newstars])
    medenew  = np.median([r.e     for r in newstars])
    medthnew = np.median([r.theta for r in newstars])
    
    # Find the boundaries of the active region on the SUB.  Regions outside
    # of the overlap of NEW and REF will have pixel values ~1e-30.  They will
    # be rectangular, so to save time just examine slices through the middle.
    # RS 2013/10/01:  NB that these regions may not be zero on the SUB!
    # Look at the REF, since swarp will set inactive pixels to zero.
    epsover = max(1e-10,10*subfitshdr["MASKVAL"])
    xslice = reffitspix[subfitshdr["NAXIS2"]/2,:]
    yslice = reffitspix[:,subfitshdr["NAXIS1"]/2]
    minx = min([i for i in range(len(xslice)) if abs(xslice[i]) > epsover])
    maxx = max([i for i in range(len(xslice)) if abs(xslice[i]) > epsover])
    miny = min([j for j in range(len(yslice)) if abs(yslice[j]) > epsover])
    maxy = max([j for j in range(len(yslice)) if abs(yslice[j]) > epsover])
    # Update the SUB header with a few keywords.
    subfitshdr.update('RUNTAG', sub.runtag, 'Pipeline run tag')
    subfitshdr.update('SUBOVER', '[{0:d}:{1:d},{2:d}:{3:d}]'
                      .format(minx,maxx,miny,maxy),
                      'Region of overlap between NEW and REF')
    subfitshdr.update('FIELD_ID', sub.field, 'SkyMapper field ID')
    # Add MINRA, MAXRA, MINDEC, MAXDEC keywords to header.  We'll use xy2sky
    # as an external, since I don't know of any Python WCS libraries that can
    # cope with TNX headers (other than pyraf...  ugh).
    minx_str, maxx_str = str(minx), str(maxx)
    miny_str, maxy_str = str(miny), str(maxy)
    p = subprocess.Popen(['xy2sky', sub.sub_fits_name,
                          minx_str, miny_str, maxx_str, maxy_str],
                          stdout=subprocess.PIPE)
    stdoutstr, stderrstr = p.communicate()
    line1, line2 = stdoutstr.split('\n')[:2]
    coo1 = ephem.Equatorial(*(line1.split()[:2]))
    coo2 = ephem.Equatorial(*(line2.split()[:2]))
    minra, mindec = min(coo1.ra, coo2.ra), min(coo1.dec, coo2.dec)
    maxra, maxdec = max(coo1.ra, coo2.ra), max(coo1.dec, coo2.dec)
    coomin = ephem.Equatorial(minra, mindec)
    coomax = ephem.Equatorial(maxra, maxdec)
    subfitshdr.update('MINRA', str(coomin.ra))
    subfitshdr.update('MINDEC', str(coomin.dec))
    subfitshdr.update('MAXRA', str(coomax.ra))
    subfitshdr.update('MAXDEC', str(coomax.dec))
    # Update the FITS header on disk.
    subfitsptr.flush()
    
    # Match coordinates of REF and NEW detections to SUB detections.
    ref_neighbors, dref = match_coords(substars,refstars,30.0)
    new_neighbors, dnew = match_coords(substars,newstars,30.0)
    
    # Loop over detections in the SUB...
    for s in substars:

        # RS 2013/10/01:  If this detection isn't even in the SUB active area,
        # just ignore it.
        if not ((minx <= s.x < maxx) or (miny <= s.y < maxy)):
            continue

        # Look for a nearest neighbor in the REF.
        # RS 2011/10/26:  We expect these to be associated with host galaxies,
        # or possibly with bright stars which might be poorly subtracted.
        idx = ref_neighbors[s.id-1]
        if idx != None:  rref = refstars[idx]
        else:            rref = None
    
        # Look for a nearest neighbor in the NEW.
        # RS 2011/11/24:  I added these because of cases where hotpants
        # convolves the NEW to match the REF; cosmic rays will then be actual
        # size in the NEW, but may look like a PSF in the SUB.
        idx = new_neighbors[s.id-1]
        if idx != None:  rnew = newstars[idx]
        else:            rnew = None
    
        # Initialize a Candidate object with the information we have so far
        cand = RealBogus(subsex=s, refsex=rref, newsex=rnew,
                         subimg=subfitspix, bkgsig=subbkgsig)
    
        # Now fill a few other quantities depending on global properties of
        # the image data and/or sets of sources found in the NEW, REF and SUB:
    
        # Ratio of SUB seeing to REF seeing
        cand.Rfwhm = subfwhm/reffwhm
        # Surface density of detections on the subtraction, detections/pix^2
        cand.goodcn = len(substars)/((maxx-minx)*(maxy-miny))
        # Did we convolve the NEW to the REF, 0/1?
        cand.subconv = (subfitshdr["CONVOL00"] == "IMAGE")*1
    
        # Fill some other quantities that depend on quantities filled above
        cand.apsig4   = cand.f4sub/cand.df4sub
        cand.apsig8   = cand.f8sub/cand.df8sub
        cand.normrms  = np.sqrt(cand.asub*cand.bsub)/subfwhm
        cand.normfwhm = cand.fwhmsub/subfwhm
        cand.Ranew    = cand.asub/medanew
        cand.Renew    = cand.esub/medenew
        cand.Dthnew   = cand.thsub-medthnew
        cand.Raref    = cand.asub/medaref
        cand.Reref    = cand.esub/mederef
        cand.Dthref   = cand.thsub-medthref
        if cand.refsrc: cand.Rfref = cand.f4sub/cand.f4ref
        else:           cand.Rfref = cand.f4sub/refflim
        if cand.newsrc: cand.Rfnew = cand.f4sub/cand.f4new
        else:           cand.Rfnew = cand.f4sub/newflim
    
        # Now run the classifier, yay!
        if not noapply:
            cand.rbscore = rbclass.score(cand)
    
        # Add it to the end!
        candlist.append(cand)

    # Write a FITS table with the information on all the candidates.
    RealBogus.write_fits_file(sub.sub_cands_name,candlist,header=subfitshdr)

    # Too many objects on this subtraction?  It's probably bad.
    reals = [rb for rb in candlist if rb.rbscore > RealBogus.rbscore_thresh]
    if len(reals) > 50:
        raise STAP_Error("Implausibly many Real objects on subtraction")
    
    # Convert the table data to a ds9 regions file which we can overlay on the
    # subtraction just for giggles.
    with open(subds9,"w") as rf:
        rf.write("global color=red width=2\n")
        rf.write("image\n")
        for c in candlist:
            if c.rbscore > RealBogus.rbscore_thresh:
                rf.write("box({0:7.2f},{1:7.2f},30,30,45) # color=green\n"
                         .format(float(c.xsub),float(c.ysub)))
            else:
                rf.write("circle({0:7.2f},{1:7.2f},10)\n"
                         .format(float(c.xsub),float(c.ysub)))
    
    # If interactive mode is selected, run a ds9 window!
    if interact:
        subkern = sub.hotpants_kernel_reg
        cmd_ext = "ds9 %s %s\n" \
               % ("-zscale -rotate 90 -zoom to 0.25 -geometry 900x600",
                 ("-fits %s -regions %s -regions %s" % (subfits,subds9,subkern)))
        print cmd_ext
        subprocess.call(shlex.split(cmd_ext))

    # Return some hopefully useful metadata
    rad2deg = 180/ephem.pi
    return { 'minra': minra*rad2deg, 'mindec': mindec*rad2deg,
             'maxra': maxra*rad2deg, 'maxdec': maxdec*rad2deg }


def main():
    """Wrapper allowing STAP.classify() to be run from the command line"""

    parser = argparse.ArgumentParser(description='Classify detected transient sources in a subtraction')
    parser.add_argument('newname', help='filename of NEW image')
    parser.add_argument('refname', help='filename of REF image ')
    parser.add_argument('subname', help='filename of SUB image')
    parser.add_argument('--timeout',type=float, default=30,
                        help='maximum running time allowed for external call (default: %(default)ss)')
    parser.add_argument('--interact',action='store_true',
                        help='run interactive scanning session')
    parser.add_argument('--noapply', action="store_true", default=False,
                        help='don\'t apply milk (for training purposes)')
    args = parser.parse_args()
    classify(args.newname, args.refname, args.subname,
             timeout=args.timeout, interact=args.interact, noapply=args.noapply)

if __name__ == "__main__":
    main()
