#!/usr/bin/env python

# ============================================================================
# solve-mosaic by Luke Shingles (lukes@mso.anu.edu.au)
# RS 2011/04/19:  Cleaned up and added to local git for subtraction pipeline.
# ----------------------------------------------------------------------------
# solve-mosaic takes a list of FITS files as arguments. For each FITS file,
# a new directory will be created containing every extension as an individual
# FITS file, for example:
#    python solve-mosaic.py myfile.fits
# will create a new directory called myfile containing
#    myfile_01.fits, myfile_02.fits, etc.
# After each extraction, solve-mosaic will run astrometry.net's solve-field
# (path may need to be changed) with RA and DEC coordinates from the header
# of myfile.fits and a default search radius of 100 arcminutes
# (set with -radius option).
# The running time of solve-field is logged to myfile.log.
# ============================================================================

import argparse
import numpy
import pyfits
import os, sys
import subprocess

parser = argparse.ArgumentParser(description='Split FITS files into extensions and run solve-field on each.')
parser.add_argument('imname', help='input filename')
parser.add_argument('outname', help='output filename')
parser.add_argument('-radius', default=100, action='store',
		    type=int, help='search radius in arcminutes (default: 100)')

args = parser.parse_args()

filename=args.imname
newname=args.outname
[basename, extension] = filename.split(".")

print "Running solve-mosaic_single.py"

hdulist = pyfits.open(filename)
try:
	ra = hdulist[0].header['CRVAL1']
	dec = hdulist[0].header['CRVAL2']
	print "use crval1/crval2 from header"
except:
	ra = hdulist[0].header['RA']
	dec = hdulist[0].header['DEC']
	print "use RA/DEC from header"
hdulist.close()

# RS 2011/04/19:  Now assumes solve-field is in the user's $PATH.
# RS 2011/04/21:  Now uses --no-fits2fits to avoid filling up /tmp.
cmd=("/usr/local/astrometry-old/bin/solve-field --parity neg --no-plots --no-fits2fits "
     + filename +
     " --scale-low 0.45 --scale-high 0.5 --scale-units app --ra "
     + str(ra) + " --dec " + str(dec) + " --radius "
     + str(args.radius / 60.0) + " --new-fits " + newname + " --overwrite")

print cmd
cmd=cmd.split(" ")
process=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
(stdoutstr,stderrstr) = process.communicate()
status = process.returncode
print stdoutstr
print stderrstr

if os.path.exists("%s.solved"%basename) and status ==0:
	os.remove("%s.axy"%basename)
	os.remove("%s.solved"%basename)
	os.remove("%s.match"%basename)
	os.remove("%s.rdls"%basename)
	os.remove("%s.wcs"%basename)
	os.remove("%s.corr"%basename)
	os.remove("%s-indx.xyls"%basename)
else:
	print "Solving WCS failed"
	sys.exit(1)
	
