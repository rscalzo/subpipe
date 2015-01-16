#!/usr/bin/env python

import argparse
import os
import re
import ephem
import pyfits
import numpy as np
from Utils import Constants
from Utils.Catalog import SExtractorDetection
from Utils.TrackableException import TrackableException as STAP_Error
from Utils.TrackableException import ExternalFailure
from Utils.Photometry import calc_seeing, calc_apcor, calc_elong
from STAP_comm import STAP_callexternal, print_cmd_line

def SEx(imname, seeing=5.0, do_apcor=False, do_seeing=False, rerun=False,
        config="STAP_SEx.sex", SEx_switches=None, timeout=None):
    """
    Runs SExtractor on an image.  For input named imname.fits, stores output
    in FITS binary table format in imname.fits.stars.

    imname:     name of SExtractor input FITS image
    seeing:     default expected FWHM of objects and of detection filter
                (will attempt to pull from SEEING keyword in FITS header)
    rerun:      option to run SExtractor twice, once to measure the seeing
                and then once again to re-detect objects and do photometry
    do_seeing:  option to calculate the seeing from detections
    do_apcor:   option to calculate aperture corrections from detections
    timeout:    optional timeout, in seconds (useful in pipeline operations).
    """

    # TODO:  Need to fix display of SEx_switches below.
    print_cmd_line("STAP_SEx.py", imname, seeing=seeing, rerun=rerun,
                   config=config, do_apcor=do_apcor, do_seeing=do_seeing,
                   SEx_switches=SEx_switches, timeout=timeout)
    if rerun:  do_seeing = True

    # RS 2014/07/02:  Verify which pyfits we're using...
    print "Using pyfits", pyfits.__version__
    
    # RS 2011/08/18:  The convolution kernel, star/galaxy neural net weights,
    # etc. are copied synchronously by the main script into the working path
    # (/ramdisk for pipeline mode running).
    
    # RS 2012/02/17:  Parse the CCD keyword from the filename for later.
    # Right now this isn't being filled elsewhere (but someday it will).
    iminfo = Constants.Filenames(imname)
    ccd = iminfo.ccd

    # RS 2011/10/20:  If the SEEING keyword exists in the image header,
    # it probably came from a previous run of SExtractor.  Use it as a more
    # reliable estimate of the appropriate detection filter FWHM.
    with pyfits.open(imname, mode='update') as hdulist:
        imghdr = hdulist[0].header
        if "SEEING" in imghdr:  seeing = float(imghdr["SEEING"])
        # RS 2014/02/05:  Fix the EQUINOX keyword, which needs to be a float
        # and not a string (new TAROS update is filling EQUINOX='2000.0').
        # RS 2014/11/07:  Weirdly, the EQUINOX keyword doesn't seem to exist
        # for old data and pyfits isn't letting us create it via assignment,
        # so this behavior of dictionaries doesn't map exactly onto headers.
        # Make a call to the "update" method instead.
        imghdr.update("EQUINOX", 2000.0, "Equinox of coordinates")
        # RS 2012/06/19:  Extract a few other keywords too...
        exptime = float(imghdr["EXPTIME"])
        airmass = float(imghdr["AIRMASS"])
        filtname = imghdr["FILTNAME"]
        minx, maxx = 0, float(imghdr['NAXIS1'])
        miny, maxy = 0, float(imghdr['NAXIS2'])
    
    # RS 2011/04/21:  Choose a good convolution kernel.
    if seeing > 0.0:  convfile = "gauss_1.5_3x3.conv"
    if seeing > 2.0:  convfile = "gauss_2.0_3x3.conv"
    if seeing > 2.5:  convfile = "gauss_2.5_5x5.conv"
    if seeing > 3.0:  convfile = "gauss_3.0_5x5.conv"
    if seeing > 3.5:  convfile = "gauss_3.0_7x7.conv"
    if seeing > 4.0:  convfile = "gauss_4.0_7x7.conv"
    if seeing > 4.5:  convfile = "gauss_5.0_9x9.conv"
    
    # RS 2012/02/27:  Generate the right filename for the star catalog.
    # RS 2012/06/29:  If some other file with this name exists, blow it away.
    starsfile = Constants.sextractor_output(iminfo.basefname)
    
    # RS 2011/04/21:  List of switches to override sextractor defaults.
    # This maps from the SkyMapper parameters, with whatever names we
    # give them and/or extract from the FITS headers, to the parameter
    # names sextractor uses.
    switches = \
        [
        # Pixel scale (arcsec/pix)
        "-PIXEL_SCALE {0:.2f}".format(Constants.Imager.pixel_scale),
        # FWHM of PSF (seeing, arcsec)
        "-SEEING_FWHM {0:.2f}".format(Constants.Imager.pixel_scale*seeing),
        # Upper saturation limits for pixels
        "-SATUR_LEVEL {0:.1f}".format(Constants.Imager.saturate_adu),
        # Detector gain
        "-GAIN {0:.3f}".format(Constants.Imager.gain),
        # Convolution kernel for extraction
        "-FILTER_NAME " + convfile,
        # Output catalog name
        "-CATALOG_NAME " + starsfile,
        ]
    if SEx_switches: switches = switches + SEx_switches.split(" ")

    # RS 2011/04/19:  Now assumes this script is in the $PATH.
    # RS 2011/08/18:  Use new switches.
    print "SExtractor:  Using convolution kernel", convfile
    cmd = "sex {0} -c {1} {2}".format(imname, config, " ".join(switches))

    # RS 2012/02/27:  We'll want to save the output here...
    status, stdoutstr = STAP_callexternal(cmd, combinestderr=True,
                                          getstdout=True, timeout=timeout)
    if status != 0:  raise ExternalFailure(cmd=cmd, exit_code=status)
    
    # RS 2011/09/14:  Similar to Brian's runsexnewfits.pl, we extract the
    # background mean and sigma from the SExtractor output.
    med, sig, thr = None, None, None
    for line in stdoutstr.split('\n'):
        # Use regular expressions to match the log output.
        mm = re.match(r'.*Background:\s+(\S+).*RMS:\s+(\S+).*Threshold:\s+(\S+)',line)
        if mm:
            med, sig, thr = mm.groups()
            med, sig, thr = float(med), float(sig), float(thr)
            break
    if med is None:
        print ""
        raise STAP_Error(msg="Couldn't find median in SExtractor output")

    # If sextractor didn't find anything, now would be a good time to die.
    stars, starhdr = SExtractorDetection.read_fits_file(starsfile)
    if len(stars) == 0:
        raise STAP_Error(msg="SExtractor detected no objects")
    else:
        print "SExtractor:  Successfully sextracted", len(stars), "objects."
        print "Background:  ", med, "  RMS:  ", sig, "  Thresh:  ", thr


    # RS 2013/10/05:  Use an external call to xy2sky to get the borders of
    # the image in RA and DEC.
    cmd = "xy2sky {0} {1} {2} {3} {4}".format(imname, minx, miny, maxx, maxy)
    status, stdoutstr = STAP_callexternal(cmd, combinestderr=True,
                                          getstdout=True, timeout=timeout)
    if status == 0:
        line1, line2 = stdoutstr.split('\n')[:2]
        coo1 = ephem.Equatorial(*(line1.split()[:2]))
        coo2 = ephem.Equatorial(*(line2.split()[:2]))
        minra, mindec = min(coo1.ra, coo2.ra), min(coo1.dec, coo2.dec)
        maxra, maxdec = max(coo1.ra, coo2.ra), max(coo1.dec, coo2.dec)
        coomin = ephem.Equatorial(minra, mindec)
        coomax = ephem.Equatorial(maxra, maxdec)
        minra_str, mindec_str = str(coomin.ra), str(coomin.dec)
        maxra_str, maxdec_str = str(coomax.ra), str(coomax.dec)
    else:
        minra_str = maxra_str = "00:00:00.00"
        mindec_str = maxdec_str = "+00:00:00.0"

    # If we don't already have a seeing estimate, we should calculate one.
    # RS 2011/10/26:  I'd argue that we use the SExtractor star/galaxy
    # neural net here to select starlike objects, but that gives weird
    # results for normal images.  The RMS (from A_IMAGE*B_IMAGE) seems to
    # fluctuate a lot.  FWHM_IMAGE, on the other hand, is very stable and
    # has a clear mode in the distribution near the (apparent) seeing.
    # So use that for all high-S/N (apsig > 20) detections on the image.
    # Make a cut on elongation (e < 1.2) to weed out really weird sources.
    # Also, don't use any sources flagged as bad (FLAGS = 0).
    tmp_seeing = None
    if do_seeing:  tmp_seeing = calc_seeing(stars, verbose=True)
    if tmp_seeing != None:  seeing = tmp_seeing
    print "Using seeing = {0:.2f} pix ({1:.2f}\") for downstream.".format(
            seeing, seeing*Constants.Imager.pixel_scale)

    # RS 2012/11/27:  Added elongation here.
    elong = calc_elong(stars, verbose=True)

    # RS 2012/06/20:  Include aperture correction here as well.  Right now
    # it's just a crappy weighted sum, but in the future we should use GPs.
    if do_apcor:
        apcor, apcor_err, apcor_chi2 = calc_apcor(stars, verbose=True)

    # Write out the new keyword values to both the SExtractor file header and
    # the main FITS image header (for easy reference).  Some of these keywords
    # will hopefully be filled by the main survey pipeline someday.
    keywords_to_update = \
    [
        ("IMAGEID", ccd, "CCD number in the mosaic"),
        ("EXPTIME", exptime, "Exposure time (sec)"),
        ("FILTNAME", filtname, "Centred Filter"),
        ("BKGMED", med, "SExtractor-calculated background level in ADU"),
        ("BKGSIG", sig, "SExtractor-calculated background sigma in ADU"),
        ("BKGTHR", thr, "SExtractor-calculated pixel threshold in ADU"),
        ("SEEING", seeing, "SExtractor-calculated median seeing in pixels"),
        ("ELONG", elong, "SExtractor-calculated median elongation (a/b)"),
    ]

    # RS 2014/02/19:  Add field boundaries in RA and DEC.
    if status==0:
        keywords_to_update += \
        [
            ('MINRA', minra_str, "Lower boundary of image in RA"),
            ('MINDEC', mindec_str, "Lower boundary of image in RA"),
            ('MAXRA', maxra_str, "Upper boundary of image in RA"),
            ('MAXDEC', maxdec_str, "Upper boundary of image in DEC"),
        ]

    # RS 2014/02/19:  Add aperture correction info, if these were performed.
    if do_apcor and apcor != None:
        for i in range(len(apcor)):
            if apcor[i] in [None, np.nan, np.inf]:  continue
            num = "{0:02d}".format(i+1)
            keywords_to_update += \
            [
                ("APCORR"+num, apcor[i],
                 "Aperture correction from ap #{0} (mag)".format(i+1)),
                ("APCERR"+num, apcor_err[i],
                 "Error on aperture correction #{0}".format(i+1)),
            ]

    # RS 2014/02/19:  Update the image header with this info.
    with pyfits.open(imname, mode="update") as hdulist:
        for kw in keywords_to_update:  hdulist[0].header.update(*kw)
        hdulist.flush()

    # RS 2014/02/19:  Propagate these keywords into the header as well.
    imgcards = imghdr.ascardlist()
    for kw in ['FIELD_ID', 'SUBFIELD', 'DATE-OBS']:
        if kw in imghdr:
            keywords_to_update.append(
                (imgcards[kw].key, imgcards[kw].value, imgcards[kw].comment))

    # RS 2014/02/19:  Update the SExtractor output FITS header with this info.
    with pyfits.open(starsfile, mode="update") as hdulist:
        for kw in keywords_to_update:  hdulist[0].header.update(*kw)
        # RS 2013/10/05:  Need to stick this in there because the .stars file
        # is generated by SExtractor, not by FITSRecord.write_fits_file().
        hdulist[1].header.update("CLASSNAM", "SExtractorDetection",
                                 "Class name of table row objects")
        hdulist.flush()

    # RS 2012/02/27:  If we're going to re-run with the new seeing estimate,
    # go ahead and do that now.  (Assuming it's worthwhile.)
    if rerun and tmp_seeing != None:
        SEx(imname, seeing=seeing, rerun=False, config=config,
            SEx_switches=SEx_switches, timeout=timeout)
    

def main():
    """Wrapper allowing STAP.SEx() to be run from the command line"""

    parser = argparse.ArgumentParser(description='Run SExtractor on an image')
    parser.add_argument('imname', 
                        help='filename of input image')
    parser.add_argument('--timeout', type=float, default=30,
                        help='maximum running time allowed for external call (default: %(default)s s)')
    parser.add_argument('--config', default="STAP_SEx.sex",
                        help='SExtractor configuration file (default: %(default)s)')
    parser.add_argument('--seeing', type=float, default=5.0,
                        help='FWHM in pix of detection filter (default: %(default)s)')
    parser.add_argument('--do_seeing', action='store_true', default=False,
                        help='calculate seeing from detections? (default: %(default)s)')
    parser.add_argument('--do_apcor', action='store_true', default=False,
                        help='calculate aperture corrections? (default: %(default)s)')
    parser.add_argument('--rerun', action='store_true', default=False,
                        help='run SExtractor twice to measure seeing? (default: no)')
    parser.add_argument('--SEx_switches',default=None,
                        help='additional SExtractor swithces')
    args = parser.parse_args()
    SEx(args.imname, seeing=args.seeing, rerun=args.rerun, config=args.config,
        do_seeing=args.do_seeing, do_apcor=args.do_apcor,
        timeout=args.timeout, SEx_switches=args.SEx_switches)
    
    
if __name__ == "__main__":
    main()
