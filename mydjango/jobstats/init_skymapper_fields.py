#!/usr/bin/env python

# ============================================================================
# RS 2012/05/18:  Initialize SkyMapper fields in a new database
# ----------------------------------------------------------------------------
# This script simply reads a text file containing SkyMapper field centers,
# then inserts them into the SkymapperField table in the Django database.
# It should be run every time we need to reinitialize the site from empty.
# ============================================================================

import ephem
from mydjango.jobstats.models import SkymapperField
from Utils.Constants import PipelinePath as CPP

def init_fields(fname):
    nfields = 0
    with open(fname) as catfile:
        for line in catfile:
            if line[0] == '#':  continue
            cols = line.strip().split()
            id = cols[0]
            ra = ephem.hours(cols[1])
            dec = ephem.degrees(cols[2])
            f = SkymapperField(
                    id=id, ra=ra*180/ephem.pi, dec=dec*180/ephem.pi)
            f.save()
            nfields += 1
    print "read_cat: {0} fields added".format(nfields)

if __name__ == "__main__":
    init_fields(CPP.etc + "/skymapper_field_centers.txt")
