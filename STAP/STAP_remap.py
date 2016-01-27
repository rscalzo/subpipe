#! /usr/bin/env python

import os
import pyfits
import argparse
from Utils.TrackableException import TrackableException as STAP_Error
from Utils.TrackableException import ExternalFailure
from STAP_comm import STAP_callexternal, print_cmd_line
import numpy as np

def remap(newname, refname, outname,
          module="swarp", config="STAP_SWarp.swarp", timeout=None,
          weightmap=None, weightmap_out=None):
    """
    Resamples one image to the coordinate system of another using SWarp
    (default) or wcsremap.

    newname:  name of FITS image containing the WCS onto which to resample;
              in subtraction pipeline use, the NEW (discovery image)
    refname:  name of FITS image to be to be resampled;
              in subtraction pipeline use, the REF (host galaxy template)
    outname:  name of output FITS image (REF resampled to NEW)
    module:   optional name of module to use, either "swarp" or "wcsremap"
    config:   optional name of SWarp configuration file, if using swarp
    timeout:  optional timeout, in seconds (useful in pipeline operations).
    """

    # First check to see whether we have the necessary WCS information in
    # the header to be able to remap the image.  Otherwise we're toast.
    try:
        with pyfits.open(refname) as hdulist:
            refhdr = hdulist[0].header
            if "WAT1_001" not in refhdr and "A_ORDER" not in refhdr:
                raise STAP_Error(msg="No WCS information in header")
    except:
        if 'coadd' in refname or 'combined' in refname:
            #coadded image using swarp would have different WCS
            pass
        else:
            raise STAP_Error(msg="Error reading file %s"%refname)
    # If the raw REF is gzip-compressed, ungzip it
    refbase, refext = os.path.splitext(refname)
    if refext == ".gz":
        cmd = "gunzip {0}".format(refname)
        status, stdoutstr = STAP_callexternal(cmd, combinestderr=True,
                                              getstdout=True, timeout=timeout)
        if status != 0:
            raise ExternalFailure(cmd=cmd, exit_code=status)
        refname = refbase

    print_cmd_line("STAP_remap.py", newname, refname, outname,
                   module=module, config=config, timeout=timeout,
                   weightmap=weightmap, weightmap_out=weightmap_out)

    if module == "wcsremap":
        cmd = "wcsremap -template {0} -source {1} -outIm {2}".format(
               newname, refname, outname)
    elif module == "swarp":
        # RS 2011/10/28:  This part is a bit fiddly -- we have to copy the NEW's
        # FITS header to a separate file before we can remap.  Use pyfits.
        # FY - really only want the WCS solution of new
        # important to preserve the SATURATE, GAIN, EXPTIME, RDNOISE, NCOADD
        refhdr=pyfits.getheader(refname)
        with pyfits.open(newname) as newptr:
            for key in ['EXPTIME','SATURATE','GAIN','RDNOISE','NCOADD','SEEING']:
                newptr[0].header.update(key,refhdr[key])
            with open(outname.replace(".fits",".head"),"w") as newrefhdr:
                newrefhdr.write(str(newptr[0].header.ascardlist()))
        # Here's the actual swarp command we're going to use.
        cmd = "swarp {0} -c {1} -IMAGEOUT_NAME {2}".format(
               refname, config, outname)
        if weightmap:
            weightmap_temp='weight.fits'
            cmd = "{0} -WEIGHT_IMAGE {1} -WEIGHT_TYPE MAP_WEIGHT -WEIGHTOUT_NAME {2}".format(cmd,weightmap,weightmap_temp)
            
    else:
        raise STAP_Error("Unknown -module switch to STAP.remap")

    status, stdoutstr = STAP_callexternal(cmd, combinestderr=True,
                                          getstdout=True, timeout=timeout)
    
    # print stdoutstr
    if status != 0:
        raise ExternalFailure(cmd=cmd, exit_code=status)
    
    if module == "swarp" and weightmap:
        if not weightmap_out: weightmap_out=outname.replace('.fits','_mask.fits')
        with pyfits.open(weightmap_temp) as hdulist:
            hdulist[0].data = np.array(1*(hdulist[0].data > 0),
                                       dtype=np.uint8)
            hdulist.writeto(weightmap_out, clobber=True)

def main():
    """Wrapper allowing STAP.remap() to be run from the command line"""

    parser = argparse.ArgumentParser(description='Solve WCS of input image and update the header')
    parser.add_argument('newname', 
                        help='filename of new image')
    parser.add_argument('refname', 
                        help='filename of reference')
    parser.add_argument('outname', 
                        help='filename of remapped reference')
    parser.add_argument('--timeout',type=float, default=300,
                        help='maximum running time allowed for external call (default: %(default)ss)')
    parser.add_argument('--module',default="swarp",
                        help='which software to use: swarp (default) or wcsremap?')
    parser.add_argument('--config',default="STAP_SWarp.swarp",
                        help='SWarp configuration file (default: %(default)s)')
    parser.add_argument('--weightmap', default=None,
                        help='optinal input of weight map for reference')
    parser.add_argument('--weightmap_out', default=None,
                        help='name of remapped weight map. Default is outname.mask')
    
    args = parser.parse_args()
    remap(args.newname, args.refname, args.outname,
          module=args.module, config=args.config, timeout=args.timeout,
          weightmap=args.weightmap, weightmap_out=args.weightmap_out)

if __name__ == "__main__":
    main()
