#!/usr/bin/env python

# ============================================================================
# RS 2011/11/21:  Classifier Support Classes -- RealBogus
# ----------------------------------------------------------------------------
# This package implements a class RealBogus, as well as an I/O API for
# reading/writing FITS tables of RealBogus instances.  I don't want to overdo
# things here, so I'm not going to worry too much about data hiding; I just
# want to make sure certain things need to be changed in at most one place.
# ----------------------------------------------------------------------------
# RS 2012/02/15:  Rewritten to take advantage of new abstract base classes
# for FITS I/O.  Otherwise it's just a huge headache.
# ============================================================================

import re
import pyfits
import numpy as np
from .Record import FITSRecord, KeywordRecord
from .Catalog import CatalogEntry


class RealBogus (FITSRecord,KeywordRecord,CatalogEntry):
    # ========================================================================
    # This class represents a detection on a subtraction together with a set
    # of scores used to drive an automated Real/Bogus classifier.  It uses
    # the same internal attribute strategy as the Record classes, but also
    # inherits from them for its FITS interface.
    # ========================================================================

    version_tag = "v0.8.0a"

    # Candidate attributes.  Most of these have to do with FITS formatting
    # (see below), except for the last column which tells us whether the
    # attribute will be used for scoring in milk.
    _cand_fields = \
    [
        # Name        Type  Units    Format  Default milk?
        [ "id",       "1J", "",      "I3",    0    , False ],
        [ "xsub",     "1E", "pixel", "F7.2",  0.0  , True  ],
        [ "ysub",     "1E", "pixel", "F7.2",  0.0  , True  ],
        [ "asub",     "1E", "pixel", "F5.2",  0.0  , False ],
        [ "bsub",     "1E", "pixel", "F5.2",  0.0  , False ],
        [ "esub",     "1E", "",      "F5.3",  0.0  , True  ],
        [ "thsub",    "1E", "deg",   "F5.1",  0.0  , True  ],
        [ "rasub",    "1D", "deg",   "F10.6", 0.0  , False ],
        [ "decsub",   "1D", "deg",   "F10.6", 0.0  , False ],
        [ "fwhmsub",  "1E", "pixel", "F5.2",  0.0  , True  ],
        [ "f4sub",    "1E", "count", "G10.5", 0.0  , True  ],
        [ "df4sub",   "1E", "count", "G10.5", 0.0  , False ],
        [ "f8sub",    "1E", "count", "G10.5", 0.0  , True  ],
        [ "df8sub",   "1E", "count", "G10.5", 0.0  , False ],
        [ "flagsub",  "1I", "",      "I3",    0    , True  ],
        [ "starsub",  "1E", "",      "F5.3",  0.0  , True  ],
        [ "refsrc",   "1L", "",      "I1",    False, False ],
        [ "xref",     "1E", "pixel", "F7.2",  0.0  , True  ],
        [ "yref",     "1E", "pixel", "F7.2",  0.0  , True  ],
        [ "aref",     "1E", "pixel", "F5.2",  0.0  , False ],
        [ "bref",     "1E", "pixel", "F5.2",  0.0  , False ],
        [ "eref",     "1E", "",      "F5.3",  0.0  , True  ],
        [ "thref",    "1E", "deg",   "F5.1",  0.0  , True  ],
        [ "fwhmref",  "1E", "pixel", "F5.2",  0.0  , True  ],
        [ "f4ref",    "1E", "count", "G10.5", 0.0  , True  ],
        [ "df4ref",   "1E", "count", "G10.5", 0.0  , False ],
        [ "flagref",  "1I", "",      "I3",    0    , True  ],
        [ "starref",  "1E", "",      "F5.3",  0.0  , True  ],
        [ "newsrc",   "1L", "",      "I1",    False, False ],
        [ "xnew",     "1E", "pixel", "F7.2",  0.0  , False ],
        [ "ynew",     "1E", "pixel", "F7.2",  0.0  , False ],
        [ "anew",     "1E", "pixel", "F5.2",  0.0  , False ],
        [ "bnew",     "1E", "pixel", "F5.2",  0.0  , False ],
        [ "enew",     "1E", "",      "F5.3",  0.0  , True  ],
        [ "thnew",    "1E", "deg",   "F5.1",  0.0  , True  ],
        [ "fwhmnew",  "1E", "pixel", "F5.2",  0.0  , True  ],
        [ "f4new",    "1E", "count", "G10.5", 0.0  , True  ],
        [ "df4new",   "1E", "count", "G10.5", 0.0  , False ],
        [ "flagnew",  "1I", "",      "I3",    0    , True  ],
        [ "starnew",  "1E", "",      "F5.3",  0.0  , True  ],
        [ "n2sig3",   "1I", "",      "I3",    0    , True  ],
        [ "n3sig3",   "1I", "",      "I3",    0    , True  ],
        [ "n2sig5",   "1I", "",      "I3",    0    , True  ],
        [ "n3sig5",   "1I", "",      "I3",    0    , True  ],
        [ "nmask",    "1I", "",      "I3",    0    , True  ],
        [ "Rfwhm",    "1E", "",      "F5.3",  0.0  , True  ],
        [ "goodcn",   "1E", "",      "G7.5",  0.0  , True  ],
        [ "subconv",  "1L", "",      "",      0.0  , True  ],
        [ "nndref",   "1E", "pixel", "F7.2",  1e+3 , True  ],
        [ "nndnew",   "1E", "pixel", "F7.2",  1e+3 , True  ],
        [ "apsig4",   "1E", "",      "F8.3",  0.0  , True  ],
        [ "apsig8",   "1E", "",      "F8.3",  0.0  , True  ],
        [ "normrms",  "1E", "",      "F6.3",  0.0  , True  ],
        [ "normfwhm", "1E", "",      "F6.3",  0.0  , True  ],
        [ "Rfref",    "1E", "",      "F6.3",  0.0  , True  ],
        [ "Raref",    "1E", "",      "F6.3",  0.0  , True  ],
        [ "Reref",    "1E", "",      "F6.3",  0.0  , True  ],
        [ "Dthref",   "1E", "deg",   "F6.3",  0.0  , True  ],
        [ "Rfnew",    "1E", "",      "F6.3",  0.0  , True  ],
        [ "Ranew",    "1E", "",      "F6.3",  0.0  , True  ],
        [ "Renew",    "1E", "",      "F6.3",  0.0  , True  ],
        [ "Dthnew",   "1E", "deg",   "F6.3",  0.0  , True  ],
        [ "rbscore",  "1E", "",      "F7.3",  0.0  , False ],
    ]

    # Construct _fits_fields from the above (which has more information).
    # Use the attribute name as the FITS column name.
    _fits_fields = [ (f[0],f[0],f[1],f[2],f[3],f[4]) for f in _cand_fields ]
    _keyword_fields = [ (f[0],f[4]) for f in _cand_fields ]
    _milk_features = [ f[0] for f in _cand_fields ]

    # rbscore threshold over which an object is considered Real
    rbscore_thresh = 40

    def __init__(self, *args, **kwargs):

        # Initialize all attributes with default values.
        if len(args) == 2:  self._init_from_fitsrow(*args)

        # If additional information is provided, this must be because we're
        # initializing a record for the first time.  Basically all the fields
        # in this case are filled through keywords.  We support the following
        # shorthand for filling SExtractor-related keywords and keywords for
        # quantities we can calculate from the SUB pixel data alone.

        if "subsex" in kwargs and kwargs["subsex"] != None:
            # Fill SExtractor quantities from the SUB
            subsex       = kwargs["subsex"]
            self.id      = subsex.id
            self.xsub    = subsex.x
            self.ysub    = subsex.y
            self.asub    = subsex.a
            self.bsub    = subsex.b
            self.esub    = subsex.e
            self.thsub   = subsex.theta
            self.ra      = subsex.ra
            self.rasub   = subsex.ra
            self.dec     = subsex.dec
            self.decsub  = subsex.dec
            self.fwhmsub = subsex.fwhm
            self.f4sub   = subsex.flux[2]
            self.df4sub  = subsex.fluxerr[2]
            self.f8sub   = subsex.flux[4]
            self.df8sub  = subsex.fluxerr[4]
            self.flagsub = subsex.flag
            self.starsub = subsex.star
        if "refsex" in kwargs and kwargs["refsex"] != None:
            # Fill SExtractor quantities from the REF
            refsex       = kwargs["refsex"]
            self.refsrc  = True
            self.xref    = refsex.x
            self.yref    = refsex.y
            self.aref    = refsex.a
            self.bref    = refsex.b
            self.eref    = refsex.e
            self.thref   = refsex.theta
            self.fwhmref = refsex.fwhm
            self.f4ref   = refsex.flux[2]
            self.df4ref  = refsex.fluxerr[2]
            self.flagref = refsex.flag
            self.starref = refsex.star
            self.nndref  = np.sqrt((self.xsub-self.xref)**2 +
                                   (self.ysub-self.yref)**2)
        if "newsex" in kwargs and kwargs["newsex"] != None:
            # Fill SExtractor quantities from the NEW
            newsex       = kwargs["newsex"]
            self.newsrc  = True
            self.xnew    = newsex.x
            self.ynew    = newsex.y
            self.anew    = newsex.a
            self.bnew    = newsex.b
            self.enew    = newsex.e
            self.thnew   = newsex.theta
            self.fwhmnew = newsex.fwhm
            self.f4new   = newsex.flux[2]
            self.df4new  = newsex.fluxerr[2]
            self.flagnew = newsex.flag
            self.starnew = newsex.star
            self.nndnew  = np.sqrt((self.xsub-self.xnew)**2 +
                                   (self.ysub-self.ynew)**2)
        if ("subimg" in kwargs and kwargs["subimg"] != None and
            "bkgsig" in kwargs and kwargs["bkgsig"] != None):
            # Fill quantities derived from SUB pixel data.
            # Only the data in a thumbnail around the candidate is needed.
            pixarr = kwargs["subimg"]
            # Don't forget:  array indices in [0,n-1], pixel coords in [1,n].
            # Also:  NAXIS1 = array index 2, NAXIS2 = array index 1.
            # We can overstep the upper index bound in numpy and it'll crop.
            i, j = self.xsub-1, self.ysub-1
            pixarr3x3 = np.array(pixarr[max(0,j-1):j+2,max(0,i-1):i+2])
            pixarr5x5 = np.array(pixarr[max(0,j-2):j+3,max(0,i-2):i+3])
            # numbers of pixels above threshold
            self.n2sig3 = len(pixarr3x3[pixarr3x3 < -2*kwargs["bkgsig"]])
            self.n3sig3 = len(pixarr3x3[pixarr3x3 < -3*kwargs["bkgsig"]])
            self.n2sig5 = len(pixarr5x5[pixarr5x5 < -2*kwargs["bkgsig"]])
            self.n2sig5 = len(pixarr5x5[pixarr5x5 < -3*kwargs["bkgsig"]])
            # number of masked pixels -- fix when we have bad pixel masks
            self.nmask  = 0

        # The other attributes depend on global properties of the NEW, REF
        # and SUB that are easier to calculate in the routine that's actually
        # filling these tables.  So punt on those, pass them in as keywords.
        self._init_from_keywords(**kwargs)
        self.ra, self.dec = self.rasub, self.decsub


    def milk_features(self):
        # This creates an array of floats for use with milk.  Include every
        # attribute except "rbscore", since that's the answer when training.
        # return np.array([getattr(self,f) for f in self._milk_features])
        return np.array([getattr(self,f[0]) for f in self._cand_fields
                         if f[0] != "rbscore"])

    def milk_label(self):
        return self.rbscore
