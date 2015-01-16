#!/usr/bin/env python

# ============================================================================
# RS 2011/04/29:  First-pass source classifier for possibly interesting
# sources found in subtractions.
# ============================================================================

import os
import re
import pyfits
import numpy
import argparse
import cPickle as pickle
from milk.supervised.randomforest import rf_model
from Utils.RealBogus import RealBogus
from Utils.Catalog import SExtractorDetection, match_coords
from Utils.DetectionMerger import DetectionMerger

# Parse the arguments
parser = argparse.ArgumentParser \
    (description='Extract X_IMAGE, Y_IMAGE from sextractor binary FITS output.')
parser.add_argument('newfname', help='NEW filename')
parser.add_argument('reffname', help='REF filename')
parser.add_argument('subfname', help='SUB filename')
parser.add_argument('--interactive', action='store_true',
                    help='interactive scanning session, True/False?')
parser.add_argument('--noapply', action="store_true", default=False,
                    help='don\'t apply milk (for training purposes)')
args = parser.parse_args()


# ----------------------------------------------------------------------------
#                         Minor class redefinition
# ----------------------------------------------------------------------------
    
    
class rf_fuzzy_model (rf_model):
    # ------------------------------------------------------------------------
    # This class just extends the milk random forest model in such a way as
    # to allow us to access the number of votes in favor of Real or Bogus.
    # Taking the continuum into account in the training process itself would
    # probably take more radical redesign than I'm willing to do right now,
    # but we can at least look at the internal level of agreement.
    # ------------------------------------------------------------------------

    def load_rf (self, binary_forest):
        self.forest = binary_forest.forest
        self.names = binary_forest.names

    def apply_fuzzy (self, features):
        rf = len(self.forest)
        votes = 1.0*sum(t.apply(features) for t in self.forest)
        return votes/rf


# ----------------------------------------------------------------------------
# Global image quantities (need to get these from user, headers, or calculate)
# ----------------------------------------------------------------------------


# Pickled machine learning classifier; currently we're using milk
pklfname = "{0}/{1}".format(os.environ["SUBETCPATH"],"randomforest.pkl")
pklfile = open(pklfname,"rb")
rbmodel = pickle.load(pklfile)
pklfile.close()

# Information about the input images
newinfo = { "fitsfn" : args.newfname }
refinfo = { "fitsfn" : args.reffname }
subinfo = { "fitsfn" : args.subfname }
for fd in (newinfo,refinfo,subinfo):
    # NEW, REF, SUB FITS image header + data
    fitsname = fd["fitsfn"]
    fitsptr = pyfits.open(fitsname,mode="update")
    fd["fitsptr"] = fitsptr
    fd["fitshdr"] = fitsptr[0].header
    fd["fitspix"] = fitsptr[0].data
    # NEW, REF, SUB SExtractor header + binary table
    starname = fitsname.replace(".fits",".fits.stars")
    starptr = pyfits.open(starname)
    fd["starfn"] = starname
    fd["starhdr"] = starptr[0].header
    fd["starcols"] = cnlist = starptr[1].columns.names
    # Some sanity checks and failsafes
    if starptr[1].data == None:
        print "FATAL:  No stars found in {0}!".format(fitsname)
        exit(1)
    elif len(starptr[1].data) > 500 and fitsname == args.subfname:
        nstars = len(starptr[1].data)
        print "WARNING:  Too many ({0}) stars found in SUB".format(nstars)
        print "We're probably going to time out here.  I'll give it a go,"
        print "but you should consider raising your SExtractor threshold,"
        print "looking for structured background in this image, etc."
    fd["starrows"] = [SExtractorDetection(r,cnlist) for r in starptr[1].data]
    # Convert the arrays to dictionaries for easier access
    # Standard deviation of background noise in ADU
    fd["bkgsig"] = fd["starhdr"]["BKGSIG"]
    # Seeing in pixels; use NEW seeing for SUB in case SUB's bad
    fd["fwhm"] = fd["starhdr"]["SEEING"]
    # Limiting flux in aperture of 1 FWHM
    fd["flim"] = numpy.pi*(fd["fwhm"]**2/4)*fd["bkgsig"]


# ----------------------------------------------------------------------------
#                         Some helpful subroutines
# ----------------------------------------------------------------------------


def median(array):
    # Finds the median of an array of floats.
    N = len(array)
    arrtmp = array
    arrtmp.sort()
    if N % 2 == 1: return arrtmp[(N-1)/2]
    else:          return (arrtmp[(N-1)/2]+arrtmp[(N+1)/2])/2


# ----------------------------------------------------------------------------
#                             The main event
# ----------------------------------------------------------------------------


# Initialize a list of candidates.
candlist = [ ]

# Make sure we have detections in the SUB!
if subinfo["starrows"] is None:
    print "No candidate found in file %s"%(substar)
    exit()

# Take a few median quantities for sources found in the REF and NEW images
medaref  = median([r.a     for r in refinfo["starrows"]])
mederef  = median([r.e     for r in refinfo["starrows"]])
medthref = median([r.theta for r in refinfo["starrows"]])
medanew  = median([r.a     for r in newinfo["starrows"]])
medenew  = median([r.e     for r in newinfo["starrows"]])
medthnew = median([r.theta for r in newinfo["starrows"]])

# Calculate photometric corrections for and between the two images
newapcor, refapcor, zeropoint = DetectionMerger.phot_calibrate \
    (new_sex_objects=newinfo["starrows"], ref_sex_objects=refinfo["starrows"])
newinfo["fitshdr"].update("APERCOR4",newapcor,
                          "Aperture correction ~ MAG_APER[6]-MAG_APER[2]")
newinfo["fitshdr"].update("ZP2REF",zeropoint,
                          "Zeropoint for 4-pixel aperture relative to REF")
refinfo["fitshdr"].update("APERCOR4",refapcor,
                          "Aperture correction ~ MAG_APER[6]-MAG_APER[2]")
refinfo["fitshdr"].update("ZP2REF",0.0,
                          "Zeropoint for 4-pixel aperture relative to REF")
subinfo["fitshdr"].update("APERCOR4",newapcor,
                          "Aperture correction ~ MAG_APER[6]-MAG_APER[2]")
subinfo["fitshdr"].update("ZP2REF",zeropoint,
                          "Zeropoint for 4-pixel aperture relative to REF")

# Find the boundaries of the active region on the SUB.  The regions outside
# of the overlap of NEW and REF will have pixel values ~1e-30.  These regions
# will be rectangular, so to save time just examine slices through the middle.
epsover = max(1e-10,10*subinfo["fitshdr"]["MASKVAL"])
xslice = subinfo["fitspix"][subinfo["fitshdr"]["NAXIS2"]/2,:]
yslice = subinfo["fitspix"][:,subinfo["fitshdr"]["NAXIS1"]/2]
minx = min([i for i in range(len(xslice)) if abs(xslice[i]) > epsover])
maxx = max([i for i in range(len(xslice)) if abs(xslice[i]) > epsover])
miny = min([j for j in range(len(yslice)) if abs(yslice[j]) > epsover])
maxy = max([j for j in range(len(yslice)) if abs(yslice[j]) > epsover])
# Also find the field ID.  The temporary fix for this is to grep on the
# SUB filename; we'll do this better soonish.
field_id = re.search('_([\d\-\+]+)_[uvgriz]_',args.subfname).group(1)
# Then update the SUB header with these keywords.
subinfo["fitshdr"].update('SUBOVER','[{0:d}:{1:d},{2:d}:{3:d}]'.format
                          (minx,maxx,miny,maxy),
                          'Region of overlap between NEW and REF')
subinfo["fitshdr"].update('FIELD_ID',field_id,'SkyMapper field ID')
subinfo["fitsptr"].flush()

# Match coordinates of REF and NEW detections to SUB detections.
ref_neighbors, dref = match_coords(subinfo["starrows"],refinfo["starrows"],5.0)
new_neighbors, dnew = match_coords(subinfo["starrows"],newinfo["starrows"],5.0)

# Loop over detections in the SUB...
for s in subinfo["starrows"]:
    # Look for a nearest neighbor in the REF.
    # RS 2011/10/26:  We expect these to be associated with host galaxies,
    # or possibly with bright stars which might be poorly subtracted.
    idx = ref_neighbors[s.id-1]
    if idx != None:  rref = refinfo["starrows"][idx]
    else:            rref = None

    # Look for a nearest neighbor in the NEW.
    # RS 2011/11/24:  I added these because of cases where hotpants convolves
    # the NEW to match the REF.  In these cases, things like cosmic rays will
    # be actual size in the NEW, but may look like a PSF in the SUB.
    idx = new_neighbors[s.id-1]
    if idx != None:  rnew = newinfo["starrows"][idx]
    else:            rnew = None

    # Initialize a Candidate object with the information we have so far
    cand = RealBogus(subsex=s, refsex=rref, newsex=rnew,
                     subimg=subinfo["fitspix"], bkgsig=subinfo["bkgsig"])

    # Now fill a few other quantities depending on global properties of the
    # image data and/or sets of sources found in the NEW, REF and SUB:

    # Ratio of SUB seeing to REF seeing
    cand.Rfwhm = subinfo["fwhm"]/refinfo["fwhm"]
    # Surface density of detections on the subtraction, detections/pix^2
    cand.goodcn = len(subinfo["starrows"])/((maxx-minx)*(maxy-miny))
    # Did we convolve the NEW to the REF, 0/1?
    cand.subconv = (subinfo["fitshdr"]["CONVOL00"] == "IMAGE")*1

    # Fill some other quantities that depend on quantities filled above
    cand.apsig4   = cand.f4sub/cand.df4sub
    cand.apsig8   = cand.f8sub/cand.df8sub
    cand.normrms  = numpy.sqrt(cand.asub*cand.bsub)/subinfo["fwhm"]
    cand.normfwhm = cand.fwhmsub/subinfo["fwhm"]
    cand.Ranew    = cand.asub/medanew
    cand.Renew    = cand.esub/medenew
    cand.Dthnew   = cand.thsub-medthnew
    cand.Raref    = cand.asub/medaref
    cand.Reref    = cand.esub/mederef
    cand.Dthref   = cand.thsub-medthref
    if cand.refsrc: cand.Rfref = cand.f4sub/cand.f4ref
    else:           cand.Rfref = cand.f4sub/refinfo["flim"]
    if cand.newsrc: cand.Rfnew = cand.f4sub/cand.f4new
    else:           cand.Rfnew = cand.f4sub/newinfo["flim"]

    # Now run the classifier, yay! -- the following is milk syntax
    # wait... don't do this if we're generating a new set to scan
    if not args.noapply:
        cand.rbscore = 100.0*rbmodel.apply(cand.milk_features())

    # Add it to the end!
    candlist.append(cand)

# Write a FITS table with the information on all the candidates.
subcand = subinfo["fitsfn"].replace(".fits",".fits.cands")
RealBogus.write_fits_file(subcand,candlist,header=subinfo["fitshdr"])

# Convert the table data to a ds9 regions file which we can overlay on the
# subtraction just for giggles.
subds9  = subinfo["fitsfn"].replace(".fits",".fits.reg")
subkern = subinfo["fitsfn"].replace(".fits",".kern.reg")
rf = open(subds9,"w")
rf.write("global color=red width=2\n")
rf.write("image\n")
for c in candlist:
    if c.rbscore > 40:
        rf.write("box(%7.2f,%7.2f,30,30,45) # color=green\n" % (c.xsub,c.ysub))
    else:
        rf.write("circle(%7.2f,%7.2f,10)\n" % (c.xsub,c.ysub))
rf.close()

# If interactive mode is selected, run a ds9 window!
if args.interactive:
    import shlex
    import subprocess
    cmd_ext = "ds9 %s %s\n" \
           % ("-zscale -rotate 90 -zoom to 0.25 -geometry 900x600",
             ("-fits %s -regions %s -regions %s" % (subfits,subds9,subkern)))
    print cmd_ext
    cmd = shlex.split(cmd_ext)
    subprocess.call(cmd)
