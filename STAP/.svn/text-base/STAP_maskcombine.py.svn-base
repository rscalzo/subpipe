#! /usr/bin/env python

import pyfits
import os
import numpy as np
import argparse
from STAP_comm import print_cmd_line

def maskcombine(masknames,outname,combinetype="AND",verbose=False):
    """Combine input masks using a specified method (default AND)
    
    masknames: list of names of FITS masks to be combined
    outname:  name of output FITS mask
    combinetype: combine the masks using OR or AND
    """

    print_cmd_line("STAP_maskcombine.py", *masknames,
                   outname=outname,combinetype=combinetype, verbose=verbose)
    
    combinetype=combinetype.upper()
    if ("OR" not in combinetype) and ("AND" not in combinetype):
        print "Input combinetype {0} is not recognized. Must be OR or AND".format(combinetype)
        return(1)
    
    try:
        image0=pyfits.getdata(masknames[0])
        head0=pyfits.getheader(masknames[0])
        key='MASK00'
        head0.update(key,os.path.basename(masknames[0]))
        if verbose:
            (imx,imy)=image0.shape
            nbad=len(np.where(image0 ==0)[0])/float(imx)/float(imy)*100.
            print "Input {0} {1:.2f} percent bad".format(key,nbad)
        if len(masknames)<2:
            print "Less than two input masks, no combine"
        else:
            for ind,file in enumerate(masknames[1:]):
                try:
                    image1=pyfits.getdata(file)
                    if "OR" in combinetype:
                        image0=image0 | image1
                    else:
                        image0=image0 & image1
                    key='MASK{0:02d}'.format(ind+1)
                    head0.update(key,os.path.basename(file))
                    if verbose:
                        nbad=len(np.where(image1 ==0)[0])/float(imx)/float(imy)*100.
                        print "Input {0} {1:.2f} percent bad".format(key,nbad)
                except:
                    print "Problem adding mask",file

        if verbose:
            nbad=len(np.where(image0 ==0)[0])/float(imx)/float(imy)*100.
            print "Output mask {0:.2f} percent bad".format(nbad)

        pyfits.writeto(outname,data=image0,header=head0,clobber=True)
        
    except:
        #make a dummy 2048 x 4096 mask
        image0=np.ones((4096,2048),dtype=np.uint8)
        print "Problem combining mask, make a fake 2048 x 4096 mask with 1 everywhere"
        pyfits.writeto(outname,data=image0,clobber=True)


def main():
    """Wrapper allowing STAP.maskcombine() to be run from the command line"""

    parser = argparse.ArgumentParser(description='Combine masks')
    parser.add_argument('masknames', nargs="+",
                        help='filenames of masks to combine')
    parser.add_argument('--outname', default="maskcombined.fits",
                        help='filename of combined mask')
    parser.add_argument('--combinetype', default="AND",
                        help='combine masks using OR or AND')
    parser.add_argument('--verbose', action='store_true', default=False,
                        help='print percentage masked area before and after combine')
    
    args = parser.parse_args()
    maskcombine(args.masknames, args.outname,combinetype=args.combinetype, verbose=args.verbose)

if __name__ == "__main__":
    main()
