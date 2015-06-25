#!/usr/bin/env python

"""
============================================================================
RS 2015/05/25:  Views and Plots for RBClassifiers
============================================================================
"""


import numpy as np
import matplotlib.pyplot as pypl


def scorehist(rbclass, sheep, goats, show=True):
    """
    Plot cumulative probability densities of RealBogus scores, to make it
    easy to visualize how placing a cut at a given R/B score helps separate
    Real from Bogus detections.

    Inputs
        rbclass:  RBClassifier instance
        sheep:  list of RealBogus instances of (ground-truth) Real detections
        goats:  list of RealBogus instances of (ground-truth) Bogus detections
    """

    # Plot the histograms first
    sheepscores = np.array([rbclass.score(cand) for cand in sheep])
    goatsscores = np.array([rbclass.score(cand) for cand in goats])
    kw = { 'bins': 50, 'range': (0.0, 100.0), 'alpha': 0.5, 'normed': True }
    sheephist, sheepbins = pypl.hist(sheepscores, color='g', **kw)[:2]
    goatshist, goatsbins = pypl.hist(goatsscores, color='r', **kw)[:2]
    pypl.xlabel("Real/Bogus Score")
    # Then plot cumulative histograms
    x = 0.5*(sheepbins[1:] + sheepbins[:-1])
    sheepvsx = np.cumsum(sheephist)
    goatsvsx = np.cumsum(goatshist)
    efficiency = 1.0*sheepvsx/len(sheepscores)
    rejection = 1.0 - bogusvsx/len(goatsscores)
    pypl.plot(x, efficiency, lw=2, ls='--', c='g',
              label='Efficiency (% of Reals kept)')
    pypl.plot(x, rejection, lw=2, ls='--', c='r',
              label='Rejection (% of Boguses chucked)')
    # Show and return the results!
    if show:
        pypl.legend(loc='best')
        pypl.show()
    return x, efficiency, rejection
