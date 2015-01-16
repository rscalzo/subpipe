#! /usr/bin/env python

import os
import sys
import pyfits
import argparse
import numpy as np
from STAP_comm import print_cmd_line
from Utils.TrackableException import TrackableException as STAP_Error
from Utils.Catalog import APASSObject, SExtractorDetection, SkymapperCalibStar
from Utils.Catalog import match_coords, CatalogEntry
from Utils.Photometry import calc_zp, make_SM_cal, merge_cal
from Utils.Photometry import calc_fluxlim, calc_fluxlim_img


def zeropoint(starsfname, calstarsfname="calstars.fits", imgfname=None):
    """Applies a calibration to zeropoint an array of SExtractorDetections"""

    print_cmd_line("STAP_zeropoint.py", starsfname,
                   calstarsfname=calstarsfname, imgfname=imgfname)

    # Read the objects in; use the file extension to guess what they are.
    # (Find a better way of doing this in the future.)
    fnext = os.path.splitext(starsfname)
    stars, hdr = SExtractorDetection.read_fits_file(starsfname)
    exptime, seeing, bkgsig = hdr['EXPTIME'], hdr['SEEING'], hdr['BKGSIG']
    filter = hdr['FILTNAME']
    filterr = filter + "_err"
    calsrc_kw, calsrc = "CALSRC_" + filter.upper(), None
    caltyp_kw, caltyp = "CALTYP_" + filter.upper(), None
    calstars, xcalstars, hdr, zp, zp_err = [ ], [ ], None, 24.0, 5.0
    print "Image seeing = {0:.2f} pix, background = {1:.1f} ADU".format(
        seeing, bkgsig)
    print "Fraction of PSF flux in 4-pixel aperture = {:.3f}".format(
        1.0 - np.exp(-0.5*(4.0/seeing)**2))

    # Get ready to match stars to the catalog
    ralist, declist = [s.ra for s in stars], [s.dec for s in stars]
    ra_min, dec_min = np.min(ralist), np.min(declist)
    ra_max, dec_max = np.max(ralist), np.max(declist)
    ra_mid, dec_mid = (ra_max + ra_min)/2, (dec_max + dec_min)/2
    cD = np.cos(np.pi*0.5*(dec_max - dec_min)/180.0)
    R = np.max([cD*(ra_max - ra_min)/2, (dec_max - dec_min)/2])

    # option 1:  Pre-existing SkyMapper frame zeropointed to APASS DR6
    if os.path.exists(calstarsfname):
        print "Reading FITS file", calstarsfname
        calstars, hdr = SkymapperCalibStar.read_fits_file(calstarsfname)
        print "Read {0} stars from {1}".format(len(calstars), calstarsfname)
        if calsrc_kw in hdr and caltyp_kw in hdr:
            calsrc, caltyp = hdr[calsrc_kw], hdr[caltyp_kw]
            zp, zp_err = calc_zp(stars, calstars, filter, verbose=True)
    # option 2:  APASS DR6 (best external catalog we have right now)
    if calsrc is None or zp_err > 1.0:
        print "Querying APASS DR6 for calibration stars..."
        calsrc, caltyp = "APASS DR6", "EXTERNAL"
        catstars = APASSObject.pgquery(ra_mid, dec_mid, R, verbose=True)
        print "{0} stars retrieved from APASS".format(len(catstars))
        zp, zp_err = calc_zp(stars, catstars, filter, verbose=True)
        xcalstars = make_SM_cal(stars, filter, zp=zp)
    # option 3:  "auto-calibrate" on SkyMapper data (differential LCs only)
    if zp_err > 1.0:
        print "No good calibration available for this filter; auto-filling"
        calsrc, caltyp = starsfname, "INTERNAL"
        xcalstars = make_SM_cal(stars, filter, exptime=exptime)
        print len(xcalstars), "auto-calibrated sources filled in"
        zp, zp_err = calc_zp(stars, xcalstars, filter, verbose=True)

    # If we've got something good enough, write it to disk.
    if zp_err < 1.0:
        calstars = merge_cal(xcalstars, calstars)
        caltyp_comment = "Type of calibration for filter " + filter
        SkymapperCalibStar.write_fits_file(calstarsfname, calstars,
            header=hdr, keywords=[(calsrc_kw, calsrc, None),
                                  (caltyp_kw, caltyp, caltyp_comment)])
    # If not, well...
    else:
        print "Total calibration fail; writing FAILED to image header"
        calsrc, caltyp = "None", "FAILED"
        zp, zp_err = 24.0, 5.0

    # Calculate a limiting flux and magnitude as well.
    faintlim = calc_fluxlim(stars)
    fluxlim50 = calc_fluxlim_img(bkgsig, CL=0.50)
    fluxlim95 = calc_fluxlim_img(bkgsig, CL=0.95)
    maglim50 = -2.5*np.log10(fluxlim50) + zp
    maglim95 = -2.5*np.log10(fluxlim95) + zp
    print "Limiting flux from stars:", faintlim
    print "Limiting flux (50%) from image background:", fluxlim50
    print "Limiting flux (95%) from image background:", fluxlim95
    print "50% CL limiting magnitude:", maglim50
    print "95% CL limiting magnitude:", maglim95

    # Finally, put this back into the FITS header of the .stars file.
    # If imgfname has been specified, i.e., the actual FITS image file,
    # add the keywords to its primary header as well.
    keywords_to_update = \
    [
        ('ZPMAG', zp, 'Zeropoint to standard aperture (mag)'),
        ('ZPMAGERR', zp_err, 'Error on zeropoint (mag)'),
        ('ZPCALSRC', calsrc, None),
        ('ZPCALTYP', caltyp, 'Type of zeropoint calibration'),
        ('FAINTLIM', faintlim, 'Flux of faintest star in image'),
        ('FLXLIM50', fluxlim50, '50% completeness flux from BKGSIG'),
        ('FLXLIM95', fluxlim95, '95% upper limit flux from BKGSIG'),
        ('MAGLIM50', maglim50, '50% completeness magnitude'),
        ('MAGLIM95', maglim95, '95% upper limit magnitude'),
    ]
    fnames_to_update = [starsfname]
    if imgfname != None:  fnames_to_update.append(imgfname)
    for fname in fnames_to_update:
        with pyfits.open(fname, mode='update') as hdulist:
            for kw in keywords_to_update:  hdulist[0].header.update(*kw)

    print "Done!"
    sys.stdout.flush()


def main():
    """Wrapper allowing STAP.getmag() to be run from the command line"""

    parser = argparse.ArgumentParser(
             description='Calculates a magnitude for an object in a SUB')
    parser.add_argument('starsfname',
                        help='SExtractor detections on image')
    parser.add_argument('--imgfname', default=None,
                        help='Original image file from which detections came')
    parser.add_argument('--calstarsfname', default="calstars.fits",
                        help='Skymapper calibration star FITS table filename (default: %(default)s)')
    args = parser.parse_args()
    zeropoint(args.starsfname, calstarsfname=args.calstarsfname, imgfname=args.imgfname)

if __name__ == "__main__":
    main()
