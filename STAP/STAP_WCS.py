#! /usr/bin/env python

import os
import re
import math
import pyfits
import argparse
from Utils.Constants import Imager
from Utils.TrackableException import TrackableException, ExternalFailure
from STAP_comm import STAP_callexternal, print_cmd_line

def WCS(imname, outname, astronet=False, timeout=None):
    """
    Solves for the World Coordinate System (WCS) of a single-CCD image
    extracted from a SkyMapper mosaic exposure.  Uses an external system call
    to Brian's WCS code, driven by $SUBPIPEHOME/bin/SM-WCS_perchip.py.

    imname:    name of input image file in FITS format
    outname:   name of new file containing input image with WCS in header
    astronet:  use astrometry.net instead of Brian's WCS code;
               assumes astrometry.net's solve_field is in the user's $PATH.
    timeout:   optional timeout, in seconds (useful in pipeline operations).
    """

    print_cmd_line("STAP_WCS.py", imname, outname,
                   astronet=astronet, timeout=timeout)

    if astronet:
        cmd = "solve-mosaic_single.py {0} {1}".format(imname,outname)
    else:
        cmd = "{0}/bin/SM-WCS-perchip.py {1} --outname {2}".format(
               os.environ['BRIANWCS'],imname,outname)

    status, stdoutstr = STAP_callexternal(cmd, combinestderr=True,
                                          getstdout=True, timeout=timeout)
    print stdoutstr
    if status != 0:
        os.system('rm -fr {0}'.format(outname))
        raise ExternalFailure(cmd=cmd, exit_code=status)

    # Check to ensure the solution makes sense
    print "Initiating sanity checks..."
    with pyfits.open(outname, mode='update') as hdulist:
        hdr = hdulist[0].header
        cd11, cd12 = hdr['CD1_1'], hdr['CD1_2']
        cd21, cd22 = hdr['CD2_1'], hdr['CD2_2']
        # Is the plate scale right?
        psx = 3600*math.sqrt(cd11**2 + cd12**2)
        psy = 3600*math.sqrt(cd21**2 + cd22**2)
        print "   psx, psy =", psx, psy, "arcsec/pix"
        if (abs(1.0-psx/Imager.pixel_scale) > 0.05 or
            abs(1.0-psx/Imager.pixel_scale) > 0.05):
            os.system('rm -fr {0}'.format(outname))
            raise TrackableException("WCS solution doesn't make sense")
        # Are the axes orthogonal?
        ctheta = (cd11*cd21 + cd12*cd22)/(psx*psy)
        theta = math.acos(ctheta)*180/math.pi
        print "   ctheta =", ctheta, "theta =", theta, "deg"
        if abs(ctheta) > 0.01:
            os.system('rm -fr {0}'.format(outname))
            raise TrackableException("WCS solution doesn't make sense")
        # What's the position angle?
        if not astronet:
            pa = math.atan2(cd12, cd11)
            print "   pa =", pa*180/math.pi, "deg"
            if abs(math.sin(pa)) > 0.02:
                os.system('rm -fr {0}'.format(outname))
                raise TrackableException("WCS solution doesn't make sense")
    print "All checks done, WCS makes sense."

def main():
    """Wrapper allowing STAP.WCS() to be run from the command line"""

    parser = argparse.ArgumentParser(
             description='Solve WCS of input image and update the header')
    parser.add_argument('imname', 
                        help='filename of input image')
    parser.add_argument('outname', 
                        help='filename of output image')
    parser.add_argument('--timeout',type=float, default=300,
                        help='maximum running time allowed for external call '
                            +'(default: %(default)ss)')
    parser.add_argument('--astronet',action='store_true', default=False,
                        help='use astrometry.net (default is Brian\'s WCS)')
    args = parser.parse_args()
    WCS(args.imname, args.outname,
        astronet=args.astronet, timeout=args.timeout)

if __name__ == "__main__":
    main()
