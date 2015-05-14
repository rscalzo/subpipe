#!/usr/bin/env python

"""
RS 2015/05/14:  Local (survey-specific) settings for SkyMapper

This module is meant, eventually, to supersede Constants as a place to put
settings that depend on the specifics of the particular survey.  That will
make it easier to swap them in and out, as originally intended.

As a bonus, we'll be able to import commonly used objects like pickled
Real/Bogus classifiers.  Eventually.
"""

import os
from Utils.RBClassifier.RBForestMilk import RBForestMilk
from Utils.RBClassifier.RBForestSklearn import RBForestSklearn

# All the pickled classifiers live in a data directory underneath this
_datadir = os.path.join(os.path.dirname(__file__), "forestfiles")

# We're going to define loaders for each type in the package, which will
# need to be called -- lazy imports would be a better way, maybe later
def _forest_loader_generator(cls, pklfname, **kwargs):
    def forest_loader():
        pklabsfname = os.path.join(_datadir, pklfname)
        forest = cls(**kwargs)
        forest.load(pklabsfname)
        return forest
    return forest_loader

# v0.8.2 -- milk random forest trained on data from 2011 August 25
#     Original implementation was trained as one-against-one, and has a
#     reversed numerical sense i.e. 100 = Bogus.
load_rbf_v082 = _forest_loader_generator(
        RBForestMilk, "randomforest_v0.8.2.pkl", invert_score=True)

# v0.8.4 alpha -- milk random forest trained on Science Verification data.
load_rbf_v083 = _forest_loader_generator(
        RBForestMilk, "randomforest_v0.8.4a.pkl")

# v0.8.3 -- sklearn random forest trained on Science Verification data.
load_rbf_v084 = _forest_loader_generator(
        RBForestSklearn, "randomforest_v0.8.3b_sklearn.pkl")
