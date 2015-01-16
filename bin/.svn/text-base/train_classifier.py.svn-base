#!/usr/bin/env python

# ============================================================================
# RS 2012/02/01:  Training a Random Forest Real/Bogus Classifier
# ----------------------------------------------------------------------------
# This code trains a random forest classifier using a Python machine learning
# package called milk.
# ============================================================================


import sys
import time
import pyfits
import argparse
import cPickle as pickle
import numpy as np
from Utils.RealBogus import RealBogus
import milk
from milk.supervised import randomforest
from milk.supervised.multi import one_against_one

# Parse the arguments
parser = argparse.ArgumentParser \
    (description='Train a random forest classifier!')
parser.add_argument('candsfiles', nargs="+",
                    help='SUB candidates FITS table filename(s)')
parser.add_argument('--Nfolds', metavar='N', type=int, default=1,
                    help='Run N-fold cross-validation (default = 10)')
args = parser.parse_args()

vtag = RealBogus.version_tag


# ----------------------------------------------------------------------------
#                         Some helpful subroutines
# ----------------------------------------------------------------------------


def pick_rbset (candlist):
    # --------------------------------------------------------------------
    # Picks a useful set of candidates to see what we've got.  We don't
    # have many Real examples, so we need to use all of them, but we
    # should cut down the Bogus stuff in our set to no more than 10x
    # the number of Real things.  Do this randomly (but reproducibly).
    # --------------------------------------------------------------------
    np.random.seed(42)
    sheep = [c for c in candlist if c.milk_label() > RealBogus.rbscore_thresh]
    goats = [c for c in candlist if c.milk_label() < RealBogus.rbscore_thresh]
    ngoat = int((len(candlist)-len(sheep))
                 *min(10.0*len(sheep)/len(candlist),1.0))
    goats = (np.random.permutation(goats))[:ngoat]
    candlist = list(sheep) + list(goats)
    print "pick_rbset:  Picked {0} real and {1} bogus candidates"\
          .format(len(sheep),len(goats))

    # Also write the list of candidates used to disk.
    candfname = "randomforest_sheepgoats_{0}.fits.cands".format(vtag)
    header = pyfits.PrimaryHDU().header
    header.update("COMMENT","File created by train_classifier.py at {0}"
                  .format(time.asctime(time.localtime(time.time()))))
    header.update("VERSION",RealBogus.version_tag,
                  "Version number of RealBogus class record")
    header.update("NUMREAL",  len(sheep), "Number of Real detections")
    header.update("NUMBOGUS", len(goats), "Number of Bogus detections")
    header.update("RBTHRESH", RealBogus.rbscore_thresh,
                  "A detection is Real if rbscore > this")
    RealBogus.write_fits_file(candfname,candlist,header=header)

    # Finally, return the list.
    return candlist


def train_rf (candlist, Nfolds = 1, return_label=False):
    # --------------------------------------------------------------------
    # Trains a random forest classifier based on a list of RealBogus
    # objects.
    # --------------------------------------------------------------------

    print "train_rf:  {0} total objects".format(len(candlist))

    # Translate the candlist into something milk can use internally.
    # Cut on 40% real as a reasonable threshold cf. conversations w/Fang.
    features = np.array([c.milk_features() for c in candlist])
    if return_label:
        labels = np.array([1*(c.milk_label() > RealBogus.rbscore_thresh)
                           for c in candlist])
    else:
        labels = np.array([1.0*c.milk_label() for c in candlist])

    # Before we start we'd like to clean up the features, get rid of linearly
    # independent features, etc., to make training faster and more reliable.
    # (rf is a binary learner, so transform it into a multi-class classifier)
    from milk.supervised.classifier import ctransforms
    from milk.supervised.normalise import chkfinite, interval_normalise
    from milk.supervised.featureselection \
         import sda_filter, featureselector, linear_independent_features
    learner = ctransforms ( chkfinite(),
                            featureselector(linear_independent_features),
                            sda_filter(),
                            one_against_one(randomforest.rf_learner()) )

    if Nfolds > 1:
        # Run N-fold cross-validation with this learner.
        print "Running cross-validation..."
        sys.stdout.flush()
        cmat, names, preds = milk.nfoldcrossvalidation\
            (features, labels, learner=learner,
             nfolds=Nfolds, return_predictions=1)
        print 'cross-validation accuracy:  ', cmat.trace()/float(cmat.sum())
        print 'fraction of Reals kept:  ', cmat[1,1]/float(cmat[1].sum())
        print 'fraction of Bogus rejected:  ', cmat[0,0]/float(cmat[0].sum())
        pklfile = open("randomforest_{0}fxval_{1}.pkl".format(Nfolds,vtag),"wb")
        pickle.dump(cmat,pklfile)
        pickle.dump(names,pklfile)
        pickle.dump(preds,pklfile)
        pklfile.close()
    else:
        # Just run normal training on the whole data set
        print "Training model..."
        sys.stdout.flush()
        model = learner.train(features,labels,return_label=return_label)
        print "Done training; pickling model..."
        sys.stdout.flush()
        pklfile = open("randomforest_model_{0}.pkl".format(vtag),"wb")
        pickle.dump(model,pklfile)
        pklfile.close()


# ----------------------------------------------------------------------------
#                             The main event
# ----------------------------------------------------------------------------


def main():

    # Read in the candidates
    candlist = [ ]
    for candsfname in args.candsfiles:
        print "Reading ", candsfname
        sys.stdout.flush()
        candlist += RealBogus.read_fits_file(candsfname)[0]
    candlist = pick_rbset(candlist)

    # Now train the classifier
    train_rf (candlist, Nfolds=args.Nfolds)


# ----------------------------------------------------------------------------
#                        This is a callable thingy
# ----------------------------------------------------------------------------


if __name__ == "__main__":
    main()
