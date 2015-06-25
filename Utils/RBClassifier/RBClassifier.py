#!/usr/bin/env python

"""
============================================================================
RS 2015/03/18:  Abstract Real/Bogus Classifier
============================================================================
"""


from ..RealBogus import RealBogus
from ..TrackableException import TrackableException


class RBExcept(TrackableException):
    pass


class RBClassifier(object):
    """Base class for Real/Bogus classifiers"""

    def __init__(self):
        """Constructor"""
        # RS 2015/04/03:  Is RealBogus really a good place to keep the score
        # threshold?  Seems to me like that ought to be a property of each
        # classifier *version*, and should be settable by the user if we want
        # to override and experiment with new thresholds.  Anyway this is in
        # RealBogus.rbscore_thresh for now so let's set it from there.
        self.thresh = RealBogus.rbscore_thresh

    def train(self, candlist):
        """Train a classifier given lists of data
        
        candlist:  list of RealBogus candidate instances with rbscores
            set to provide labels (0 = Bogus, 100 = Real)
        """
        pass

    def nfoldxvalid(self, candlist, n=5):
        """Perform n-fold cross-validation given lists of data
        
        candlist:  list of RealBogus candidate instances with rbscores
            set to provide labels (0 = Bogus, 100 = Real)
        n:  integer number of folds
        """
        pass

    def score(self, cand):
        """Given a RealBogus instance, return a score between 0 and 100"""
        return 0.0

    def predict(self, cand):
        """Given a RealBogus instance, return 0 if Bogus, 100 if Real"""
        return 100.0*(self.score(cand) > self.thresh)

    def load(self, fname):
        """Loads from a pickle or HTF5 file
        
        NB that loading to/from pickle files requires the classifier objects
        to be picklable, which may depend on which library we're using.
        So I'll let the derived classes sort this out.
        """
        pass

    def save(self, fname):
        """Saves to a pickle or HTF5 file

        NB that loading to/from pickle files requires the classifier objects
        to be picklable, which may depend on which library we're using.
        So I'll let the derived classes sort this out.
        """
        pass
