#! /usr/bin/env python

import argparse
import pyfits
import numpy as np
import subprocess
import shlex
import os

def roundint(imgname, outname=None, gzip=False, mask=None, maskflag=0, timeout=None):
    """Round off pixel values in floating-point images to 16-bit integers

    Compute the BSCALE and BZERO keywords in the FITS header according
    to the noise properties.  Will also gzip if asked.
    """

    mode = "update"
    if outname != None:  mode = "readonly"
    with pyfits.open(imgname, mode=mode) as hdulist:
        header, pix = hdulist[0].header, hdulist[0].data
        # pyfits has a feature to calculate BSCALE and BZERO automatically,
        # but it's broken in our version.  So, we'll do it ourselves.
        if 'SOFTNAME' in header and header['SOFTNAME'] == 'HOTPanTS':
            #negative value can be as large as positive.
            bscale, bzero = 1.0, 0.0
        else:
            maxpix, minpix, medpix = np.max(pix), np.min(pix), np.median(pix)
            drange=2**17
            if (maxpix-minpix+1)>drange:
                print "The dynamic range is too large, trim to sensible values"
                minallow=-100.
                if minpix<minallow:minpix=minallow
                maxallow=minpix+drange-1.
                if maxpix>maxallow:maxpix=maxallow
                #bad1=pix<minallow
                #bad2=pix>maxallow
                #print "under/overflow pixels",bad1.sum(),bad2.sum()
            bscale, bzero = 1.0*(maxpix-minpix+1)/(2**16), (maxpix+minpix)/2
        print "calculated:  bscale, bzero =", bscale, bzero
        # Calculate the image pixel values before and after, so we can do
        # a sanity check on how the image changed by converting to int.
        before = np.array(pix)
        hdulist[0].scale('int16', 'old', bscale=bscale, bzero=bzero)
        after = np.array(hdulist[0].data)*bscale + bzero
        diff = (before - after) #.ravel()
        # Disregard obviously saturated or over-/underflowed pixels
        bad = (np.abs(diff) > 1.01)
        diff = diff[np.invert(bad)]
        diffmin, diffmax, diffrms = diff.min(), diff.max(), diff.std()
        diffmean, diffmed = np.mean(diff), np.median(diff)
        print "diff:  min = {0:.6f}, max = {1:.6f}, " \
            "mean = {2:.6f}, med = {3:.6f}, rms = {4:.6f}".format(
            diffmin, diffmax, diffmean, diffmed, diffrms)
        print "# of bad translations:  {0}/{1}".format(bad.sum(), bad.shape[0]*bad.shape[1])
        if not mask is None and bad.sum()>0:
            #update the mask for bad transformations
            with pyfits.open(mask,mode='update') as maskhdu:
                bad0=(maskhdu[0].data==maskflag)
                maskhdu[0].data[bad]=maskflag
                bad1=(maskhdu[0].data==maskflag)
                print "pixel mask updated, # of bad pixels changed from %d to %d"%(bad0.sum(),bad1.sum())
        # Write the result out to disk.  If the user gave a separate name for
        # the output file, write it there.  Else just clobber the original.
        if outname:
            if gzip:  outname = outname.replace('.gz','')
            hdulist.writeto(outname,clobber=True)
        else:
            hdulist.flush()
            outname = imgname

    # ...and compress it, if the user so desires.
    if gzip:
        size_before = os.path.getsize(outname)*1.0/(2**20)
        print "BEFORE size of file:  {0:.1f} MB".format(size_before)
        print "gzip returned:  ", subprocess.call(['gzip', outname])
        outname += ".gz"
        size_after = os.path.getsize(outname)*1.0/(2**20)
        print "AFTER size of file:   {0:.1f} MB".format(size_after)


def main():
    """Wrapper allowing STAP.hotpants() to be run from the command line"""
    parser = argparse.ArgumentParser(description=
            'Round off pixel values in floating-point images to 16-bit integers.')
    parser.add_argument('imgname', help='filename of input image')
    parser.add_argument('--outname', default=None,
                        help='optional filename of output image')
    parser.add_argument('--gzip', action='store_true', default=False,
                        help='gzip the output image? (default=no)')
    parser.add_argument('--timeout',type=int, default=15,
                        help='maximum running time allowed (default: %(default)s)')
    parser.add_argument('--maskname', default=None,
                        help='optional bad pixel mask for possible update')
    parser.add_argument('--maskflag', type=int, default=0,
                        help='flag for pixels with bad transformation (default: %(default))')
    args = parser.parse_args()
    roundint(args.imgname,
             outname=args.outname, gzip=args.gzip, mask=args.maskname,
             maskflag=args.maskflag, timeout=args.timeout)

if __name__ == "__main__":
    main()
