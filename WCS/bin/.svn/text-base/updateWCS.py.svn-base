#!/usr/bin/env python

# ============================================================================
# RS 2011/11/02:  updateWCS.py -- convert Brian's WCS image to IRAF TNX
# ----------------------------------------------------------------------------
# updateWCS.py - put a IRAF TNX WCS onto an image
# Brian Schmidt, Australian National University April 2006 in PERL     
# moved to Python Aug 2011
# Taken over by RS 2011/11/02, cleaned up, debugged to match behavior of
# updateWCS.pl, but using pyRAF to make it pipelineable.
# ============================================================================

from pyraf import iraf
import glob
import re
import numpy as np
import os
import argparse
import pyfits

# Environment setup.  May need to refactor for pipeline use.
startingdir=os.getcwd()

# Defaults for some parameters...
RAref='INDEF'
DECref='INDEF'
lngrefunits='INDEF'
latrefunits='INDEF'
refpoint='coords'
proj='tan'
crossterms='half'

# Parse command-line arguments.
parser = argparse.ArgumentParser \
        (description='ADD IRAF TNX WCS to FITS File given series of inputs')
parser.add_argument('filename',
        help='input filename (single extension image.fits)')
parser.add_argument('-o', '--order', dest="iraforder", type=int, default=2,
        help='order of IRAF fit (2 is default and is a plane tangent)')
parser.add_argument('-s', '--sig', dest="sig", type=float, default=2.5,
        help='sigma reject (2.5)')
parser.add_argument('-i', '--interactive', action="store_true", default=False,
        help='make fit interactive')
parser.add_argument('-c', '--clean', type=float, default=None,
        help='clean to this level in arcsecs')
parser.add_argument('-r', '--refpt',
        help='RA DEC of CRVALs')
parser.add_argument('-m', '--matchlist', default='default',
        help='RA,DEC match to x,y list')

args=parser.parse_args()
# filename=os.path.splitext(args.imname)[0]
if args.refpt:
    refptarr = args.refpt.split(' ')
    RAref=refptarr[0]
    DECref=refptarr[1]
if RAref!='INDEF':
    lngrefunits='hours'
    latrefunits='degrees'
    refpoint='user'
if args.interactive:  interactkw = "yes"
else:                 interactkw = "no"

resultfilename = args.filename.replace(".fits",".fits.wcs")
matchlist=args.matchlist
if matchlist == 'default':
    matchlist = args.filename.replace(".fits",".fits.wcsmatch")

if args.iraforder > 2:
    proj='tnx'

# delete existing WCS 
# fix potential problem with deleting keywords in place
os.system("delwcs -n %s" % args.filename)
newfile=args.filename.replace(".fits","e.fits")
os.system("mv %s %s"%(newfile,args.filename))

# Find the dimensions of the image.
fitsptr = pyfits.open(args.filename)
xaxis = fitsptr[0].header["NAXIS1"]
yaxis = fitsptr[0].header["NAXIS2"]
fitsptr.close()

# Clean up any previous IRAF output.  IRAF will crash if it's left lying
# around, since imcoords.ccmap doesn't have a "clobber" option.
if os.path.exists(resultfilename):  os.remove(resultfilename)

# Create a pyRAF image.
iraf.image(_doprint=0)
iraf.image.imcoords(_doprint=0)

# Definition of a function we're about to use below.  It turns out that though
# Brian's original updateWCS.pl uses two calls to ccmap, they have the same
# parameters, and there are a lot of them so I'm abbreviating the call.
def iraf_ccmap_call ():
    iraf.image.imcoords.ccmap \
    (
        matchlist,                      # List of matched stars w/residuals
        'wcs',
        solution='wcsout',
        images=args.filename,           # Image to run
        results=resultfilename,
        xcolumn=1,   ycolumn=2,         # Which columns in matchlist are x,y?
        lngcolumn=3, latcolumn=4,       # and which are RA, DEC?
        insystem='j2000',               # Epoch of input RA, DEC
        refpoint=refpoint,              # User-defined reference RA, DEC?
        lngrefunits=lngrefunits,        # If so:  refpoint RA units
        latrefunits=latrefunits,        #         refpoint DEC units 
        lngref=RAref,                   #         refpoint RA
        latref=DECref,                  #         refpoint DEC
        xmin=1, xmax=xaxis,             # Pixel coords of image area where
        ymin=1, ymax=yaxis,             #    WCS is valid
        lngunit='degrees',              # RA units in matchlist
        latunit='degrees',              # DEC units in matchlist
        project=proj,                   # Type of projection (TAN or TNX)
        fitgeom='general',              # Plate solution geometry
        function='polynomial',          # Function to fit residuals
        xxorder=args.iraforder,         # Polynomial coefficient settings...
        xyorder=args.iraforder,
        xxterms=crossterms,
        yxorder=args.iraforder,
        yyorder=args.iraforder,
        yxterms=crossterms,
        reject=2.5,                     # Sigma-clip residuals at this sigma
        update='yes',                   # Update WCS in input FITS header
        pixsyst='logical',              # Work in single-image pixel coords
        interactive=interactkw          # Do this interactively? (def=False)
    )


if args.clean:
    # "Clean to this type of residual."
    # RS:  I've punted on this because I can't figure out which residuals
    # file Brian was using to do this kind of rejection.  My guess is that
    # it's done automatically in IRAF and the --clean option is legacy.
    print "Sigma-clipping not yet supported in this version, sorry."
    iraf_ccmap_call()
else:
    # "Take what we get from IRAF."
    iraf_ccmap_call()
