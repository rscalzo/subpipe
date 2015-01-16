#! /usr/bin/env python2.6

import sys
import numpy as np
import pyfits
import argparse
from Utils.TrackableException import ExternalFailure
from STAP_comm import STAP_callexternal, print_cmd_line

def quantiles(dist, qq):
    """Finds the quantiles qq of a distribution dist"""
    if len(dist) < 1:  return 0.0*np.array(qq)
    sortdist = list(dist)
    sortdist.sort()
    sortdist = np.array(sortdist)
    # interpolate across fractional bins
    fi = np.array(qq)*(len(sortdist)-1)
    i = np.array([np.int(f) for f in fi])
    return (fi-i)*sortdist[i+1] + (1-(fi-i))*sortdist[i]

def flatcombine(imglist, outname, swarp_config="STAP_flats.swarp",
                sex_config="STAP_SEx.sex", weightmap=None, timeout=None,mask=None):
    """Median-combines a set of images to form a superflat

    imglist:  list of names of FITS images to be combined into the flat
    outname:  name of output FITS image (combined "superflat")
    timeout:  optional timeout, in seconds (useful in pipeline operations)
    sex_config:   optional name of SExtractor configuration file
    swarp_config:   optional name of SWarp configuration file
    """

    print_cmd_line("STAP_flatcombine.py", *imglist, outname=outname,
                   sex_config=sex_config, swarp_config=swarp_config,
                   weightmap=weightmap, timeout=timeout)

    # Prep, just for now:  Add FLATNORM = 1/MEDPIX to all the headers.
    print "Adding FLATNORM to header, in case it's not there"
    for fname in imglist:
        with pyfits.open(fname, mode="update") as hdulist:
            header = hdulist[0].header
            if "MEDPIX" in header:
                header.update("FLATNORM", 1.0/header["MEDPIX"],
                              "SWarp flux scaling for flat field = 1/FLATNORM")

    # First use swarp to combine the images into a median stack.
    imgnames = " ".join(imglist)
    cmd = "swarp {0} -c {1} -IMAGEOUT_NAME {2}".format(
           imgnames, swarp_config, outname)
    status, stdoutstr = STAP_callexternal(cmd, combinestderr=True,
                                          getstdout=True, timeout=timeout)
    if status != 0:  raise ExternalFailure(cmd=cmd, exit_code=status)


    if mask:
        maskdata=pyfits.getdata(mask)
        bad=np.where(maskdata==0)
        if len(bad[0])>0:
            with pyfits.open(outname,mode="update") as hdulist:
                data=hdulist[0].data
                data[bad]=1
                hdulist[0].header.update('MASK',mask.split('/')[-1],
                                         'Bad pixel mask applied')
                hdulist.flush()

    # If we're just making flat fields, we don't need to do anything else.
    # The following only happens if we're making bad pixel masks.
    if not weightmap:  return

    # Run sextractor on each image to generate a background-subtracted image
    # with any sources (if they exist) zeroed out.
    objfree_imglist = [ ]
    for fname in imglist:
        objfree_fname = fname.replace(".fits","_objfree.fits")
        cmd = "sex {0} -c {1} -CHECKIMAGE_NAME {2} -CHECKIMAGE_TYPE -OBJECTS"\
              " -MASK_TYPE CORRECT -FILTER_NAME gauss_5.0_9x9.conv".format(
              fname, sex_config, objfree_fname)
        status, stdoutstr = STAP_callexternal(cmd, combinestderr=True,
                                              getstdout=True, timeout=timeout)
        if status != 0:  raise ExternalFailure(cmd=cmd, exit_code=status)
        objfree_imglist.append(objfree_fname)

    # Now combine these images with SWarp using a chi-square combine.
    objfree_imgnames = " ".join(objfree_imglist)
    cmd = "swarp {0} -c {1} -IMAGEOUT_NAME {2} -COMBINE_TYPE CHI2".format(
           objfree_imgnames, swarp_config, "chisq.fits")
    status, stdoutstr = STAP_callexternal(cmd, combinestderr=True,
                                          getstdout=True, timeout=timeout)
    if status != 0:  raise ExternalFailure(cmd=cmd, exit_code=status)

    # To create the mask, run sextractor to detect all single pixels which
    # deviate from the average of their neighbors by more than 5 sigma.
    # That is, this finds pixels which have wider spread (larger chi-square)
    # than a typical pixel in that part of the image.
    cmd = "sex {0} -c {1} -CHECKIMAGE_NAME {2} -CHECKIMAGE_TYPE OBJECTS "\
          " -FILTER N -DETECT_THRESH 5 -DETECT_MINAREA 1".format(
          "chisq.fits", sex_config, "maskobjs.fits")
    status, stdoutstr = STAP_callexternal(cmd, combinestderr=True,
                                          getstdout=True, timeout=timeout)
    if status != 0:  raise ExternalFailure(cmd=cmd, exit_code=status)

    # Finally, convert the detection image to a weight map (type = uint8)
    # so that 1 = good pixel, 0 = should be masked/interpolated/ignored.
    print "Integerizing and switching polarity of mask...  ",
    sys.stdout.flush()
    with pyfits.open("maskobjs.fits") as hdulist:
        hdulist[0].data = np.array(1 - 1*(hdulist[0].data > 1e-30),
                                   dtype=np.uint8)
        hdulist.writeto(weightmap, clobber=True)
    print "done."

def main():
    """Wrapper allowing STAP.flatcombine() to be run from the command line"""

    parser = argparse.ArgumentParser(description='Solve WCS of input image and update the header')
    parser.add_argument('imgnames', nargs="+",
                        help='filenames of images to combine')
    parser.add_argument('--outname', default="combined.fits",
                        help='filename of median combined flat')
    parser.add_argument('--mask',
                        help='filename of a badpix map, All bad pixel will be set to 1.')
    parser.add_argument('--weightmap',
                        help='filename of weight map (not-bad pixels)')
    parser.add_argument('--swarp_config',default="STAP_flats.swarp",
                        help='SWarp configuration file (default: %(default)s)')
    parser.add_argument('--sex_config',default="STAP_SEx.sex",
                        help='SExtractor configuration file (default: %(default)s)')
    parser.add_argument('--timeout',type=float, default=300,
                        help='maximum running time allowed for external call (default: %(default)ss)')
    
    args = parser.parse_args()
    flatcombine(args.imgnames, args.outname,
                sex_config=args.sex_config, swarp_config=args.swarp_config,
                weightmap = args.weightmap, timeout=args.timeout,mask=args.mask)

if __name__ == "__main__":
    main()
