#!/usr/bin/env python

# RS:  Reads X_IMAGE and Y_IMAGE from sextractor output binary FITS table
# and returns an ASCII list for input into hotpants.

import pyfits
import argparse
parser = argparse.ArgumentParser \
    (description='Extract X_IMAGE, Y_IMAGE from sextractor binary FITS output.')
parser.add_argument('imname', help='input filename')
args = parser.parse_args()

hdulist = pyfits.open(args.imname)
tbdata = hdulist[1].data
colnames = hdulist[1].columns.names
col_x = colnames.index("X_IMAGE")
col_y = colnames.index("Y_IMAGE")
col_m = colnames.index("MAG_BEST")
col_s = colnames.index("FWHM_IMAGE")
col_d = colnames.index("CLASS_STAR")
print colnames
for r in tbdata:
    print "%8.3f %8.3f %5.2f %6.2f %f" % (r[col_x],r[col_y],r[col_m],r[col_s],r[col_d])
