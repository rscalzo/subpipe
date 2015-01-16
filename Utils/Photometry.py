#!/usr/bin/env python

# ============================================================================
# RS 2012/06/20:  Separate library for common photometry algorithms
# ----------------------------------------------------------------------------
# This is where we put things like seeing measurements, aperture corrections,
# and zeropointing algorithms derived from SExtractorDetections.
# ============================================================================

import re
import shlex
import subprocess as sp
import pyfits
import numpy as np
from scipy.stats import scoreatpercentile
from scipy.optimize import curve_fit
from scipy.linalg import lstsq
from scipy.special import erfinv
from Catalog import SExtractorDetection, SkymapperCalibStar, match_coords
from Constants import Imager

# GIGO flags:  divide by zeroes, logs of negative numbers, etc.
dievals = [None, np.inf, -np.inf, np.nan]

# RS 2014/02/19:  Index of SExtractor aperture for zeropointing.
# I've been using the 4-pixel aperture, but the 6-pixel aperture should
# perform better in general given the seing we've been getting.
sexzpidx = 3


def is_good_star(s, snrlim=20.0):
    """Boolean test for whether a single star is useful for calibration
    
    All of the routines below require some idea of what a "good" star is.
    This filter selects stars suitable for calculating aggregate quantities
    pertaining to an exposure.  Our criteria are:
        no SExtractor warning flags,
        aperture fluxes are all well-defined and non-negative,
        signal-to-noise in seeing-matched aperture greater than snrlim.
    Returns a Boolean flag, True if the star passed the cut, False otherwise.
    """
    return (s.flag == 0 and all([f not in dievals for f in s.flux])
                        and all([f > 0 for f in s.flux])
                        and (s.flux[sexzpidx]/s.fluxerr[sexzpidx] > snrlim))

def select_good_stars(stars, snrlim=20.0):
    """Given a list of SExtractorDetections, returns the good ones."""
    return filter(lambda s: is_good_star(s, snrlim), stars)

def where_good_stars(stars, snrlim=20.0):
    """Given a list of SExtractorDetections, returns INDICES of good ones."""
    return np.array([is_good_star(s, snrlim) for s in stars])

def calc_seeing(stars, verbose=False, quartiles=False):
    """Calculate the seeing from a list of SExtractorDetections.

    stars:      list of SExtractorDetections
    verbose:    flag to turn status messages to stdout on/off
    quartiles:  flag to return [25%, 50%, 75%] quantiles, not just median
    """
    # RS 2011/10/26:  I'd argue that we use the SExtractor star/galaxy neural
    # net here to select starlike objects, but that gives weird results for
    # normal images.  The RMS (from A_IMAGE*B_IMAGE) seems to fluctuate a lot.
    # FWHM_IMAGE, on the other hand, is very stable and has a clear mode in
    # its distribution near the (apparent) seeing.  So use that for all
    # high-S/N (apsig > 20), non-flagged (FLAGS == 0) detections.  Might also
    # consider a shape cut at some point if it's possible to differentiate
    # badly-vibrated stars from galaxies.
    good_stars = select_good_stars(stars)
    if len(good_stars) < 5:
        if verbose:
            print "Warning:  not enough good objects to calculate seeing!"
        return None
    # We've got enough good stars, so carry on.
    qq = [float(scoreatpercentile([d.fwhm for d in good_stars], q))
          for q in (25, 50, 75)]
    if verbose:
        print "Seeing (from {0:d} objects):  {1:.2f} pix".format(
              len(good_stars), float(qq[1]))
    if quartiles:  return qq
    return qq[1]

def calc_elong(stars, verbose=False, quartiles=False):
    """Calculate the elongation from a list of SExtractorDetections.

    stars:      list of SExtractorDetections
    verbose:    flag to turn status messages to stdout on/off
    quartiles:  flag to return [25%, 50%, 75%] quantiles, not just median
    """
    # RS:  Just pull the unwindowed a/b of each object and find the median.
    # Windowed parameters would be better than unwindowed, but this is just
    # a first pass and I'm being lazy.  I'll fix later if it's worth it.
    good_stars = select_good_stars(stars)
    if len(good_stars) < 5:
        if verbose:
            print "Warning:  not enough good objects to calculate elongation!"
        #return None
        return -1
    # We've got enough good stars, so carry on.
    qq = [float(scoreatpercentile([d.e for d in good_stars], q))
          for q in (25, 50, 75)]
    if verbose:
        print "Aspect ratio (from {0:d} objects):  {1:.2f}".format(
              len(good_stars), float(qq[1]))
    if quartiles:  return qq
    return qq[1]

def calc_apcor_gauss(r1, r2):
    """Calculates aperture correction from r1 to r2 assuming a Gaussian PSF
            
    r1 and r2 are in units of the image seeing.  The aperture correction is
    given in magnitudes, such that m2 = m1 + apcor.
    """
    return -2.5*np.log10((1 - np.exp(-0.5*r2**2))/(1 - np.exp(-0.5*r1**2)))

def calc_apcor(stars, verbose=False):
    """Calculate the aperture correction(s) from a list of SExtractorDetections.

    stars:    list of SExtractorDetections
    verbose:  flag to turn status messages to stdout on/off
    """
    # RS 2012/06/20:  Include aperture correction here as well.  Right now
    # it's just a crappy weighted sum, but in the future we could use GPs.
    good_stars = np.array(select_good_stars(stars))
    if len(good_stars) < 5:
        if verbose:
            print "Warning:  not enough good objects to correct apertures!"
        return None, None, None
    # We've got enough good stars, so carry on.
    flux = np.array([s.flux for s in good_stars])
    fluxerr = np.array([s.fluxerr for s in good_stars])
    nstars, naps = flux.shape
    # Take flux ratios in different apertures.  They'll be correlated,
    # but let's ignore that just for now.
    Rf = flux[:,:-1]/flux[:,-1:]
    # Rfvar = (fluxerr/flux)**2 + (fluxerr[:,-1:]/flux[:,-1:])**2
    Rfvar = (fluxerr[:,-1:]**2 - fluxerr[:,:-1]**2
            + fluxerr[:,-1:]**2*(1 - flux[:,:-1]/flux[:,-1:]))/flux[:,-1:]**2
    Rferr = np.sqrt(Rfvar)
    # These will be noisy, presumably because not everything is really a
    # point source.  Histograms show a very tight core along with noisy
    # tails, presumably of extended objects on the CCD.  So we can use the
    # median first, and then get rid e.g. of 5-sigma deviants.
    Rfmed = np.median(Rf, axis=0)
    Rfpull = (Rf[:,3] - Rfmed[3])/Rferr[:,3]
    keep = np.arange(nstars)[abs(Rfpull) < 5.0]
    good_stars, Rf, Rferr = good_stars[keep], Rf[keep,:], Rferr[keep,:]
    nstars = len(keep)
    # Now do a weighted mean of the remaining ones.  Assume the fluctuations
    # that are left are statistical, but use the actual dispersion of the
    # flux ratios to take correlations into account.
    wi = 1.0/Rferr**2
    Rfavg = np.average(Rf, weights=wi, axis=0)
    mag_inst, apcor = -2.5*np.log10(Rf), 2.5*np.log10(Rfavg)
    apcor_chi2 = np.average(((Rf-Rfavg)/Rferr)**2, axis=0)
    disp = np.std(mag_inst + apcor, axis=0)
    avgerr = np.mean(1.0857*Rferr, axis=0)
    apcor_err = disp/np.sqrt(nstars - 1)
    if verbose:
        print "Aperture corrections (from {0} objects):".format(nstars)
        for i in np.arange(naps - 1):
            cols = (apcor, apcor_err, apcor_chi2, disp, avgerr)
            print "Ap. corr'n {0} = {1:.3f} +/- {2:.3f}, chisq_nu = {3:6.3f}, " \
                  "dispersion {4:.3f} mag, mean err {5:.3f}".format(
                  i, *[float(f[i]) for f in cols])
    # If any corrections are bad, just set them to stock values
    badidx = np.any([np.isinf(apcor), np.isinf(apcor_err),
                     np.isnan(apcor), np.isnan(apcor_err)], axis=0)
    apcor[badidx], apcor_err[badidx], apcor_chi2[badidx] = 0.0, 99.9, 1.0e+8
    return apcor, apcor_err, apcor_chi2

def calc_fluxlim(stars):
    """Calculates the limiting flux of an array of SExtractorDetections.
    
    This procedure assumes that the surface density of stars on the sky is a
    steeply increasing function of magnitude, so that the limiting magnitude
    is close to the mode of the distribution, which should be close to the
    median because most of these stars will be pretty faint!
    """
    # Find the faintest good object detected, and its signal-to-noise
    good_stars = select_good_stars(stars, snrlim=0.1)
    imin = np.array([s.flux[sexzpidx] for s in good_stars]).argmin()
    smin = good_stars[imin]
    fluxmin, fluxerrmin = smin.flux[sexzpidx], smin.fluxerr[sexzpidx]
    return fluxmin

def calc_fluxlim_img(bkgsig, aperture=6.0, siglim=6.0, CL=0.95):
    """Calculates the limiting flux of an image from header properties.

    This function tries to get accurate point-source limiting flux, in ADU,
    from bulk image properties, like the seeing and background noise.

    seeing:    seeing in pixels
    bkgsig:    standard deviation of the background in ADU
    siglim:    limiting significance in sigma
    aperture:  limiting flux aperture diameter in pixels (assumed circular)
    """
    # Background noise within aperture:  assume 50% completeness at siglim.
    # Calculate how far below threshold this confidence level represents.
    dsig = np.sqrt(2.0) * erfinv(2.0*CL - 1.0)
    return (siglim + dsig) * np.sqrt(np.pi) * (aperture/2.0) * bkgsig

def calc_zp(stars, calstars, filter, verbose=False,
            airmasses=None, apcors=None):
    """Apply a calibration to zeropoint an array of SExtractorDetections.

    stars:      list of SExtractorDetections
    calstars:   list of SkymapperCalibStars
    filter:     one-character name of SkyMapper filter for 'stars' detections
    verbose:    flag to turn status messages to stdout on/off
    airmasses:  list of airmasses for SExtractorDetections
    apcors:     list of aperture corrections for SExtractoDetections
    """
    # If the user supplied a list of airmasses, check it's the same length
    # as the list of stars, otherwise ignore.
    if airmasses is not None and len(airmasses) != len(stars):
        print "WARNING:  len(airmasses) = {0} != {1} = len(stars);".format(
                len(airmasses), len(stars)),
        airmasses = None
    if apcors is not None and len(apcors) != len(stars):
        print "WARNING:  len(apcors) = {0} != {1} = len(stars);".format(
                len(apcors), len(stars)),
        apcors = None
    # Now pick the good stars.  Limit the supplied aperture corrections and
    # airmass arrays to the same set of stars.  Ideally these should have
    # been part of the original SExtractorDetection design, but addind it
    # would mean reprocessing a lot of data at this stage...
    goodidx = where_good_stars(stars)
    if np.sum(goodidx) < 1:
        print "FATAL:  no good stars left!  Returning defaults."
        return 24.0, 5.0
    good_stars = np.array(stars)[goodidx]
    if airmasses is not None:
        print "including extinction term"
        airmasses = np.array(airmasses)[goodidx]
    else:
        print "ignoring extinction term"
    if apcors is not None:
        print "including aperture corrections"
        apcors = np.array(apcors)[goodidx]
    else:
        print "ignoring aperture corrections"
    # And only pick the good calibration stars too, while we're at it.
    goodattrs = [filter + ext for ext in ('', '_err', 'zpc', 'zpc_err')]
    good_calstars = np.array([cs for cs in calstars
        if np.all([getattr(cs, attr) not in dievals for attr in goodattrs])
        and getattr(cs, filter + '_err') < 0.2])
    best_calstars = np.array([cs for cs in good_calstars
        if cs._nodata not in [getattr(cs, attr) for attr in goodattrs]
        and getattr(cs, filter + 'zpc_err') < 0.2])
    # Match the calibration stars to the input star catalog.
    # Try to use the stars with color information first; require at least
    # 10 for a good solution (this is pretty arbitrary, will tune later).
    # NB:  type(matchidx != 0) == np.array, type(matchidx != None) == bool
    def where_ok(x):
        return np.array([i for i in range(len(x)) if x[i] != None])
    matchidx = np.array(match_coords(good_stars, best_calstars)[0])
    incl = where_ok(matchidx)
    if len(incl) < 10:
        # If that falls through, require 5 stars with no color information,
        matchidx = np.array(match_coords(good_stars, good_calstars)[0])
        incl = where_ok(matchidx)
        if len(incl) < 5:
            # If we don't even have that, we should get out of the kitchen.
            if verbose:
                print "Only", len(incl), "good cross-matches to the catalog!"
                print "Not enough for a good zeropoint; returning defaults."
            return 24.0, 5.0
        else:
            use_colors = False
            use_calstars = [good_calstars[i] for i in matchidx[incl]]
    else:
        use_colors = True
        use_calstars = [best_calstars[i] for i in matchidx[incl]]
    use_stars = good_stars[incl]
    if airmasses is not None:
        airmasses = airmasses[incl]
    if apcors is not None:
        apcors = apcors[incl]

    # Right!  So we made it this far.  Here we go with fitting.
    # Unpack all the magnitudes and their errors into numpy arrays.
    calmaglist = np.array([getattr(s, filter) for s in use_calstars])
    calmagerrs = np.array([getattr(s, filter+"_err") for s in use_calstars])
    calzpclist = np.array([getattr(s, filter+"zpc") for s in use_calstars])
    calzpcerrs = np.array([getattr(s, filter+"zpc_err") for s in use_calstars])
    maglist = -2.5*np.log10([s.flux[sexzpidx] for s in use_stars])
    if apcors is not None:
        maglist += apcors
    magerrs = -2.5*np.log10([1.0 - s.fluxerr[sexzpidx]/s.flux[sexzpidx]
                             for s in use_stars])
    zpmag = calmaglist - maglist
    zpmagerrs = np.sqrt(calmagerrs**2 + magerrs**2)

    # If we have colors, include them.  If not, don't.  In either case,
    # fit for a systematic dispersion such that the reduced chi-square ~ 1.
    # Let's try linear least squares since curve_fit is being fussy...
    sigma = zpmagerrs
    if use_colors and airmasses is not None:
        zpf = lambda c, X, zp, kc, kext: zp + kc*c + kext*X
        A = np.c_[np.ones(len(calzpclist)), calzpclist, airmasses]
    elif use_colors:
        zpf = lambda c, X, zp, kc: zp + kc*c
        A = np.c_[np.ones(len(calzpclist)), calzpclist]
    else:
        zpf = lambda c, X, zp: zp
        A = np.c_[np.ones(len(calzpclist))]
    popt, resid, rank, singvals = lstsq(A/sigma[:, None], zpmag/sigma)
    resids = zpf(calzpclist, airmasses, *popt) - zpmag
    zpstat = sigma.std()/np.sqrt(len(sigma) - 1)
    # Reiterate to make chisq_nu <= 1, if necessary
    sysvar = resids.var() - sigma.mean()**2
    print "resids.std = {0:.3f}, sigma.mean = {1:.3f}".format(
            resids.std(), sigma.mean())
    if sysvar <= 0:
        zpsys = 0.0
        chisq_nu = np.mean((resids/sigma)**2)
    else:
        zpsys = np.sqrt(sysvar)
        sigma = np.sqrt(sigma**2 + sysvar)
        popt, resid, rank, singvals = lstsq(A/sigma[:, None], zpmag/sigma)
        resids = zpf(calzpclist, airmasses, *popt) - zpmag
        chisq_nu = np.mean((resids/sigma)**2)
    zp, zp_err = popt[0], np.sqrt(zpstat**2 + zpsys**2)
    if verbose:
        print "Calibration (from {0} stars): ".format(len(incl)),
        print "zp = {0:.3f} +/- {1:.3f} (stat) +/- {2:.3f} (sys)".format(
                zp, zpstat, zpsys)
        if use_colors: 
            print "kc = {0:.3f}".format(popt[1])
        if use_colors and airmasses is not None:
            print "kext = {0:.3f}".format(popt[2])
        print "chisq_nu = {0:.3f}".format(chisq_nu)
    return zp, zp_err

def make_SM_cal(stars, filter, exptime=None, zp=None):
    """Generate calibration stars from a list of SExtractorDetections.

    Transfers a zeropoint (from calc_zp) to a list of SExtractorDetections
    and returns them as SkymapperCalibStars.  If no zeropoint exists, which
    may happen if we don't have any external catalog coverage of the area
    from which the SExtractorDetections are drawn, do the best we can with
    with guesstimated zeropoints; this will at least allow us to get reliable
    *differential* light curves with the right shape.

    stars:    a list of SExtractorDetections
    filter:   a one-character SkyMapper filter name
    exptime:  exposure time of stars in seconds
    seeing:   seeing of stars in pixels
    zp:       calibrated zeropoint, if it exists
    """
    # See whether we need to calculate a fake ZP
    magkw, magerrkw = filter, filter + '_err'
    if zp is None:
        if exptime is None:
            print "SMcal:  FAIL -- I need an exposure time or a zeropoint"
            return [ ]
        else:
            zp = Imager.fake_zp[magkw] + 2.5*np.log10(exptime/110.0)
    # Return a list of SkymapperCalibStars made from all stars that have
    # sensible instrumental magnitudes.  For the moment we're just using
    # Calculate instrumental magnitudes for all the stars.  Currently we're
    # aperture #2 (the 4-pixel aperture), and it's hard-coded.  :(
    # In the future this may become an aperture-corrected standard magnitude.
    good_stars = select_good_stars(stars, snrlim=5.0)
    return [SkymapperCalibStar(
            **{ 'ra': s.ra, 'dec': s.dec,
                magkw: -2.5*np.log10(s.flux[sexzpidx]) + zp,
                magerrkw: -2.5*np.log10(
                    1.0 - s.fluxerr[sexzpidx]/s.flux[sexzpidx]) })
            for s in good_stars]

def merge_cal(calstars, pcalstars, improve=False):
    """
    Merges two lists of SkymapperCalibStars, adding data from calstars into
    pcalstars.  If there's a star in calstars that doesn't exist in pcalstars,
    add it in.  Resolve conflicts via the "improve" keyword below.

    calstars:   list of SkyMapperCalibStars to merge into pcalstars
    pcalstars:  list of pre-existing SkyMapperCalibStars for this field
    improve:    a Boolean flag that tells the process how to merge objects
                when a magnitude for the same objects in the same filter
                exists in both calstars and pcalstars.
                If improve=False (default), pcalstars always wins conflicts.
                If improve=True, use the most precise measurement.
    """
    # First match coordinates between the two catalogs.
    matchidx, dr = match_coords(calstars, pcalstars)
    for i in range(len(matchidx)):
        # If no match, add the unmatched source to pcalstars.
        # Otherwise, merge attributes into the matching pcalstars item.
        if matchidx[i] == None:  pcalstars.append(calstars[i])
        else:  pcalstars[matchidx[i]].merge(calstars[i]) #, improve=True)
    return pcalstars
