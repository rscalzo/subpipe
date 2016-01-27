#! /usr/bin/env python

import re
import sys
import pyfits
import argparse
from Utils import Constants
from Utils.Catalog import SExtractorDetection
from Utils.TrackableException import TrackableException as STAP_Error
from Utils.TrackableException import ExternalFailure
from STAP_comm import STAP_callexternal, print_cmd_line

def hotpants(newname, refname, subname, sexstamps=False, timeout=None):
    """Runs hotpants to produce a subtraction

    Andrew Becker's "hotpants" code detects sources in two WCS-matched images,
    convolves the two images to the same PSF, and subtracts them to reveal
    transient and variable sources.

    newname:    name of FITS image from which to subtract (discovery img NEW)
    refname:    name of FITS image to subtract (host galaxy template REF)
    subname:    name of output FITS image (NEW - REF)
    sexstamps:  optional flag to use the SExtractor output from the NEW to
                extract "stamps" from which to compute the convolution kernel.
    timeout:  optional timeout, in seconds (useful in pipeline operations).
    """

    print_cmd_line("STAP_hotpants.py", newname, refname, subname,
                   sexstamps=sexstamps, timeout=timeout)

    # RS 2012/02/28:  Set up the collection of filenames for this subtraction.
    kernel_reg_fname = subname.replace(".fits",".kern.reg")

    # RS 2011/04/21:  List of switches to hotpants.  They are legion.
    # In general, switches "-t*" refer to the REF and "-i*" to the NEW.
    
    #set some default values
    newsat=refsat=Constants.Imager.saturate_adu
    newgain=refgain=Constants.Imager.gain
    newread=refread=Constants.Imager.read_noise
    #update with values in header
    newhdr=pyfits.getheader(newname)
    refhdr=pyfits.getheader(refname)
    if 'SATURATE' in newhdr:
        if newhdr['SATURATE']>0: newsat=newhdr['SATURATE']
    if 'SATURATE' in refhdr:
        if refhdr['SATURATE']>0: refsat=refhdr['SATURATE']
    if 'GAIN' in newhdr:
        if newhdr['GAIN']>0: newgain=newhdr['GAIN']
    if 'GAIN' in refhdr:
        if refhdr['GAIN']>0: refgain=refhdr['GAIN']
    if 'RDNOISE' in newhdr:
        if newhdr['RDNOISE']>0: newread=newhdr['RDNOISE']
    if 'RDNOISE' in refhdr:
        if refhdr['RDNOISE']>0: refread=refhdr['RDNOISE']
    #lower valid count
    #default
    tl=-100
    il=-100
    #estimate from background-10*sig
    if 'BKGMED' in refhdr and 'BKGSIG' in refhdr:
        tl=refhdr['BKGMED']-10*refhdr['BKGSIG']
    if 'BKGMED' in newhdr and 'BKGSIG' in newhdr:
        il=newhdr['BKGMED']-10*newhdr['BKGSIG']
        
    switches = \
    [
        # Upper saturation limits for pixels to subtract
        "-tu {0:.1f}".format(refsat),
        "-iu {0:.1f}".format(newsat),
        # Upper saturation limits for pixels to us in the PSF matching kernel
        "-tuk {0:.1f}".format(refsat*0.9),
        "-iuk {0:.1f}".format(newsat*0.9),
        # Lower valid data count
        "-tl {0:.1f}".format(tl),
        "-il {0:.1f}".format(il),
        # Detector gain
        "-tg {0:.1f}".format(refgain),
        "-ig {0:.1f}".format(newgain),
        # Detector read noise
        "-tr {0:.1f}".format(refread),
        "-ir {0:.1f}".format(newread),
        # Outer limit half-width of kernel in pixels
        "-r {0:d}".format(Constants.Imager.rpsf),
        # Output kernel file (we're saving these)
        "-savexy " + kernel_reg_fname,
        # Normalize to NEW image
        "-n i",
        # Use histogram convolution merit, not sigma or variance
        #"-fom h",
        # Polynomial order of spatial variation of kernel
        "-ko 2",
        # Polynomial order of spatial variation of sky background
        "-bgo 2",
        # Force convolution on template
        #"-c t",
        # Use as many stamps to minimize gap
        "-nsx 40",
        "-nsy 80",
        # Lower threshold for substamps
        "-ft 10",
    ]

    # Extract seeing from NEW image, and use it to scale the kernel
    kernfpars = [(6, 0.5), (4, 1.0), (2, 2.0)]
    with pyfits.open(newname) as hdulist:
        seeing = hdulist[0].header['SEEING']
    """
    with pyfits.open(refname) as hdulist:
        refseeing = hdulist[0].header['SEEING']
    kernswitch = "-ng {0}".format(len(kernfpars))
    if refseeing > seeing:
        raise STAP_Error("REF seeing worse than NEW seeing")
    kwidth = 0.5*(seeing**2 - refseeing**2)**0.5
    for kfp in kernfpars:
        kernswitch += " {0} {1:.3f}".format(kfp[0], kwidth*kfp[1])
    switches.append(kernswitch)
    """
    
    if sexstamps:
        # RS:  If -sexstamps flag is set, read a list of sources sextracted
        # from the NEW image and use these to build the kernel.  Otherwise,
        # just let hotpants divide the image up into substamps itself.
        newstarsfits = Constants.sextractor_output(newname)
        newstarstxt = newstarsfits + ".txt"
        objects = SExtractorDetection.read_fits_file(newstarsfits)[0]
        if len(objects) == 0:
            raise STAP_Error("-sexstamps was passed to hotpants,"
                             " but no SExtractor detections")
        with open(newstarstxt,"w") as stamp_coords:
            for r in objects:
                #if r.flag ==0 and fwhm > 1 and r.fwhm > seeing*0.5:    
                stamp_coords.write("{0:8.3f} {1:8.3f}\n".format
                                   (float(r.x),float(r.y)))
        switches.append("-ssf " + newstarstxt)
    else:
        # Number of distinct regions in which to build independent kernels
        # These flags will solve for a different kernel in each of eight square
        # regions 512x512 pix in the original NEW.
        # switches = switches + ["-nrx 2", "-nry 4"]
        # These flags will solve for a different kernel in each amp of the NEW.
        switches = switches + ["-nrx 1", "-nry 1"]
    
    # RS 2011/04/18:  Now assumes this is in the $PATH.
    # RS 2011/04/21:  Now uses switches that I gave to Luke a few months ago.
    cmd = "hotpants -inim {0} -tmplim {1} -outim {2} {3}".format(
           newname, refname, subname, " ".join(switches))

    # FY: problem of missing logs!!!!
    # RS 2011/04/28:  Solved.  Both wcsremap and hotpants write most of their
    # log output to stderr.  We could hack this but I don't think it's worth it.
    status, stdoutstr = STAP_callexternal(cmd, combinestderr=True,
                                          getstdout=True, timeout=timeout)
    if status != 0:  raise ExternalFailure(cmd=cmd, exit_code=status)
    # RS 2012/04/03:  Sadly, hotpants is too dumb to know when it's failed for
    # some failure conditions.  Parse the output to make sure it didn't die.
    if not re.search("Convolving...", stdoutstr):
        raise STAP_Error(msg="hotpants reported SUCCESS, but didn't convolve")

def main():
    """Wrapper allowing STAP.hotpants() to be run from the command line"""

    parser = argparse.ArgumentParser(description='Solve WCS of input image and update the header')
    parser.add_argument('newname', help='filename of new image')
    parser.add_argument('refname', help='filename of reference')
    parser.add_argument('subname', help='filename of difference image')
    parser.add_argument('--timeout',type=float, default=300,
                        help='maximum running time allowed for external call (default: %(default)ss)')
    # RS 2011/04/28:  added option to use sextractor to define substamps
    parser.add_argument('--sexstamps',action='store_true',default=False,
                        help='sextract substamps, True/False?')
    args = parser.parse_args()
    hotpants(args.newname, args.refname, args.subname,
             sexstamps=args.sexstamps, timeout=args.timeout)

if __name__ == "__main__":
    main()
