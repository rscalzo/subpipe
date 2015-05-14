#!/usr/bin/env python

"""
==============================================================================
RS 2015/04/11:  sklearn Random Forest RBClassifier
==============================================================================
"""


import numpy as np
import cPickle as pickle
from sklearn import cross_validation
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import LinearSVC
from sklearn.pipeline import Pipeline
from .RBClassifier import RBClassifier, RBExcept
from ..RealBogus import RealBogus


class RBForestSklearn(RBClassifier):
    """Updated random forest classifier based on sklearn"""

    def __init__(self, invert_score=False):
        """Constructor"""
        RBClassifier.__init__(self)
        self._model = None

    def _setup_model(self):
        """Initializes a milk "learner" instance"""

        # Use sklearn's Pipeline functionality to clean up and select
        # features; here "models" and "learners" are basically the same
        self._model = Pipeline(
            [
                ('feature_selection', LinearSVC(penalty="l1", dual=False)),
                ('classification', RandomForestClassifier(n_estimators=101))
            ])

    def train(self, candlist):
        """Train a classifier given lists of data
        
        candlist:  list of RealBogus candidate instances with rbscores
            set to provide labels (0 = Bogus, 100 = Real)
        """

        features = np.array([c.milk_features() for c in candlist])
        labels = np.array([1.0*c.milk_label() for c in candlist])
        self._setup_model()
        self._model.fit(features, labels)

    def nfoldxvalid(self, candlist, n=5, verbose=False):
        """Perform n-fold cross-validation given lists of data
                        
        candlist:  list of RealBogus candidate instances with rbscores
            set to provide labels (0 = Bogus, 100 = Real)
        n:  integer number of folds
        """

        features = np.array([c.milk_features() for c in candlist])
        labels = np.array([1.0*c.milk_label() for c in candlist])
        self._setup_model()
        scores = cross_validation.cross_val_score(self._model,
                features, labels, cv=n)
        if verbose:
            print '    Accuracy for folds:  ', scores
        return scores

    def score(self, cand):
        """Given a RealBogus instance, return a score between 0 and 100"""
        if self._model is None:
            raise RBClassExcept("model is not trained")
        # self._model.predict predicts either 0 (Bogus) or 100 (Real),
        # based on the class with the highest predicted probability.
        # self._model.predict_proba returns a 2-D array of probabilities
        # which sum to 1 along axis=1.  So we want element [0,1].
        return 100.0*self._model.predict_proba(cand.milk_features())[0,1]

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
