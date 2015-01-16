#!/usr/bin/env python

# ============================================================================
# RS 2011/11/07:  Brian's WCS Solution, Per Chip
# ----------------------------------------------------------------------------
# This script is meant to parallelize the execution of WCS.pl and updateWCS.py
# in Brian's original SM-WCS.py.  These steps take a long time when done in
# series, but really involve doing the same thing to each chip in the mosaic.
# Brian's original comments:
#    matches stars on each chip against ucac3
#    fits a ZPN + distortion matrix to a grid
#    then puts a TNX-5order onto each image
#    this provides the best astrometry possible within
#    the confines of the WCS system that can be
#    easily interpreted by Sextractor, SWARP, DS9 etc.
#    Based on CASU Code written by Mike Irwin and Jim Lewis        
#    Brian Schmidt, Australian National University April 2006 in PERL     
#    moved to Python May 2011
#    21Jul11 - added TNX to be calculated from same CRVAL at center of array
#    21Jul11 - made fit go from 2mass rather than ucac2 - slow, but better
# ============================================================================

import argparse
import pyfits
import os
import sys
import subprocess
import ephem
import shlex
import math

# RS 2011/11/02:  We already had to create an environment variable $BRIANWCS
# in which to store global configuration files (in $BRIANWCS/etc/).  We might
# as well not clutter up $PATH with the executables in $BRIANWCS/bin/.

brianwcspath = os.environ["BRIANWCS"]
if not os.path.exists(brianwcspath):
    print "SM-WCS.py:  Whoops!  Please setenv BRIANWCS and try again."
    exit(-1)


############################################################
# MAIN PROGRAM
#

startingdir=os.getcwd()

parser = argparse.ArgumentParser(description='Take SkyMapper FITS files and put rough WCS on it')
parser.add_argument('imname',
                    help='input filename (assumed to be a single chip, not a mosaic')
parser.add_argument('--outname',
                    help='optional output filename, image with WCS added')
parser.add_argument('--ucac2',action='store_true',
                    help='use ucac2 rather than 2mass')
parser.add_argument('--noresid',action='store_true',
                    help='straight zpn rather than zpn and distortions')
parser.add_argument('--crval', type=float, nargs=2, default = [-1,-1],
                    help='directly provide crval1 and crval2 in degrees')

args=parser.parse_args()

# Helpful subroutines to run things

class DoageErr:
    def __init__(self,code,msg=""):
        self.code = code
        self.msg = msg

def doage(cmd):
    retcode = subprocess.call(shlex.split(cmd))
    if retcode != 0:  raise DoageErr(retcode)


# Pre-digest some of these flags into text keyword arguments for scripts.
if args.noresid:  fitresid = "";
else:             fitresid = "-fitresid "
if args.ucac2:    astrocat = "-ucac2 "
else:             astrocat = "-2mass "
if args.outname:
    doage ("cp %s %s" % (args.imname,args.outname))
    file = args.outname
else:
    file = args.imname

# Originally enclosed in the (1..32) chip loop:

# "First do a WCS match, but do not apply."  If the match file doesn't exist,
# we know this stage must have failed, so die horribly.
maskname=args.imname.replace('.fits','_mask.fits')
if os.path.exists(maskname):
    cmd = "%s/bin/WCS.pl %s -usewcs %s -donotapply -force -mask %s" %(brianwcspath,file,astrocat, maskname)
else:
    cmd = "%s/bin/WCS.pl %s -usewcs %s -donotapply -force" %(brianwcspath,file,astrocat)
    
head=pyfits.getheader(file)
if head['FILTNAME'] in 'v':
    cmd = cmd + " -thresh 2"

doage(cmd)

matchfile = file.replace(".fits",".fits.wcsmatch")
if not os.path.exists(matchfile):
    #try it again
    cmd = cmd + " -tol 0.015 -roterr 5"
    doage(cmd)
    matchfile = file.replace(".fits",".fits.wcsmatch")
    if not os.path.exists(matchfile):
        sys.stderr.write('SM-WCS-perchip.py:  WCS.pl failed on image '+ file +'\n')
        exit(-1)

# "This fits a grid."  Again, die if the output doesn't exist.
doage ("%s/bin/brian_fitwcs %s %s %s -outputgrid 2048 4096 50 50 -outfile" %
       (brianwcspath,file,matchfile,fitresid))

gridfile = file.replace(".fits",".fits.wcsmatch.grid")
if not os.path.exists(gridfile):
    sys.stderr.write('SM-WCS-perchip.py:  brian_fitwcs failed on image '+ file +'\n')
    exit(-1)

# "Finally use the grid produced by above to put a 5 TNX into the header
#  via IRAF.  Yuk!"  NB:  The below assumes that RACENT and DECCENT used in
# the original "roughwcs" call have been stored in the FITS header as
# CRVAL1 and CRVAL2, which is true in the current version of the code.
# (It makes sense since the center of the focal plane is the tangent point.)
# They'll be stored in decimal degrees; convert to sexigesimal.

if args.crval[0] > 0:
    RACENT = args.crval[0]
    DECCENT = args.crval[1]
    print "CRVAL from command line:  %f %f" % (RACENT,DECCENT),
else:
    fitsptr = pyfits.open(file)
    RACENT = fitsptr[0].header["CRVAL1"]
    DECCENT = fitsptr[0].header["CRVAL2"]
    fitsptr.close()
    print "CRVAL from FITS header:  %f %f" % (RACENT,DECCENT),
deg2rad = math.pi/180.0
cent = ephem.Equatorial(deg2rad*RACENT, deg2rad*DECCENT, epoch=2000)
print "= %s %s" % (cent.ra,cent.dec)
doage ("%s/bin/updateWCS.py %s --matchlist %s --order 5 --refpt '%s %s'" %
       (brianwcspath,file,gridfile,cent.ra,cent.dec))
