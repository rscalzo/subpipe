#!/usr/bin/env python

"""RS 2014/02/20:  Can't think where else to put this right now..."""

import pyfits

def headerprop(fitsfname1, fitsfname2, kwlist):
    """Propagates keywords from one FITS header to another"""
    kw2upd8 = [ ]
    # Get the keywords from the first header
    hdr1 = pyfits.getheader(fitsfname1)
    cards = hdr1.ascardlist()
    for kw in kwlist:
        if kw in hdr1:
            kw2upd8.append(
                (cards[kw].key, cards[kw].value, cards[kw].comment))
    # Then propagate to the second header
    with pyfits.open(fitsfname2, mode='update') as hdulist:
        hdr2 = hdulist[0].header
        for kw in kw2upd8:  hdr2.update(*kw)
