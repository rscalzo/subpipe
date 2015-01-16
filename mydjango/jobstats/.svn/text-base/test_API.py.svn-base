#!/usr/bin/env python

import sys
import numpy as np
from mydjango.jobstats.models import SkymapperField, SkymapperPointing
from mydjango.jobstats.STAP_API import SMfield_lookup, register_pointing

ra, dec = float(sys.argv[1]), float(sys.argv[2])
print "Using ra, dec = ", ra, dec

myfield = SMfield_lookup(ra, dec)
arcdist = np.sqrt((ra-myfield.ra)**2 + (dec-myfield.dec)**2)
print "Closest field is", myfield,
print "(arc dist = {0:.4f} deg)".format(arcdist)

# fitsfile = sys.argv[3]
# print "Registering pointing ", fitsfile
# register_pointing(fitsfile)
