#!/usr/bin/env python

"""
==============================================================================
RS 2015/04/04:  milk Random Forest RBClassifier
==============================================================================
"""


import numpy as np
import cPickle as pickle
from milk import nfoldcrossvalidation
from milk.supervised import randomforest
from milk.supervised.multi import one_against_one
from milk.supervised.classifier import ctransforms
from milk.supervised.normalise import chkfinite, interval_normalise
from milk.supervised.featureselection \
    import sda_filter, featureselector, linear_independent_features
from .RBClassifier import RBClassifier, RBExcept
from ..RealBogus import RealBogus


class RBForestMilk(RBClassifier):
    """The original random forest classifier based on milk"""

    def __init__(self, invert_score=False):
        """Constructor
        
        invert:  bool saying whether a score of 100 is Real or Bogus.
            This setting is needed when using the classifier with older
            pickle files where the model was accidentally trained backwards.
        """
        RBClassifier.__init__(self)
        self._invert_score = invert_score
        self._model = None

    def _setup_learner(self):
        """Initializes a milk "learner" instance"""

        # Use milk's "ctransforms" pipeline functionality to clean up and
        # select features, to make training faster and more reliable
        self._learner = ctransforms(
                chkfinite(),
                featureselector(linear_independent_features),
                sda_filter(),
                randomforest.rf_learner()
            )

    def train(self, candlist):
        """Train a classifier given lists of data
        
        candlist:  list of RealBogus candidate instances with rbscores
            set to provide labels (0 = Bogus, 100 = Real)
        """

        features = np.array([c.milk_features() for c in candlist])
        labels = np.array([1.0*c.milk_label() for c in candlist])
        self._setup_learner()
        self._model = self._learner.train(features, labels, return_label=False)

    def nfoldxvalid(self, candlist, n=5, verbose=False):
        """Perform n-fold cross-validation given lists of data
                        
        candlist:  list of RealBogus candidate instances with rbscores
            set to provide labels (0 = Bogus, 100 = Real)
        n:  integer number of folds
        """

        features = np.array([c.milk_features() for c in candlist])
        labels = np.array([1.0*c.milk_label() for c in candlist])
        self._setup_learner()
        cmat, names, preds = nfoldcrossvalidation(features, labels,
                learner=self._learner, nfolds=n, return_predictions=1)
        if verbose:
            print '    Reals tagged as Reals:  ', cmat[1,1]/float(cmat[1].sum())
            print '    Reals tagged as Bogus:  ', cmat[1,0]/float(cmat[1].sum())
            print '    Bogus tagged as Reals:  ', cmat[0,1]/float(cmat[0].sum())
            print '    Bogus tagged as Bogus:  ', cmat[0,0]/float(cmat[0].sum())
        return cmat

    def score(self, cand):
        """Given a RealBogus instance, return a score between 0 and 100"""
        if self._model is None:
            raise RBClassExcept("model is not trained")
        score = self._model.apply(cand.milk_features())
        if self._invert_score:
            return 100.0*(1.0 - score)
        else:
            return 100.0*score

    def load(self, fname):
        """Loads from a pickle or HTF5 file
        
        NB that loading to/from pickle files requires the classifier objects
        to be picklable, which may depend on which library we're using.
        So I'll let the derived classes sort this out.
        """
        try:
            with open(fname) as pklfile:
                self._model = pickle.load(pklfile)
        except:
            raise RBExcept("couldn't read {0} in {1}.load()"
                           .format(fname, self.__class__.__name__))

    def save(self, fname):
        """Saves to a pickle or HTF5 file

        NB that loading to/from pickle files requires the classifier objects
        to be picklable, which may depend on which library we're using.
        So I'll let the derived classes sort this out.
        """
        try:
            with open(fname, "wb") as pklfile:
                pickle.dump(self._model, pklfile)
        except:
            raise RBExcept("couldn't write {0} in {1}.save()"
                           .format(fname, self.__class__.__name__))
