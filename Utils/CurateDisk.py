#!/usr/bin/env python

"""
RS 2015/05/07:  Pipeline Disk Curation

This infrastructure is written to keep our disks relatively clean by
automatically identifying files that need to be deleted.  Eventually,
we should be able to run this as its own workflow, possibly during the
day as the flip side of STAP_ingest_from_markle.py --production
"""

import os
import re
import sys
import glob
import ephem
import pyfits
import numpy as np
import mydjango.followup.models as fu
from Utils.Constants import PipelinePath as CPP
from Utils.Constants import follow_types, Imager
from Utils.Photometry import calc_fluxlim_img, calc_apcor_gauss
from STAP import zeropoint


# Query the database for all active SNe -- we only want to do this once
_active_sne = [ ]
for type in follow_types:
    try:
        djtype = fu.TransientType.objects.get(type=type)
        _active_sne += [sn for sn in djtype.transient_set.all()]
    except fu.TransientType.DoesNotExist:
        continue
# Form tuples as a shorthand for which field and subfield each SN is in
_sn_locations = [(sn.field.id, sn.ccd) for sn in _active_sne]


def calc_maglim50(bkgsig, seeing, exptime, filter):
    """
    Quick limiting magnitude for Gaussian PSF
    """
    return (-2.5*np.log10(
            calc_fluxlim_img(bkgsig, siglim=7.0, CL=0.50)/exptime) +
            calc_apcor_gauss(6.0/seeing, np.inf) + Imager.zpinst[filter])

def is_old(hdr, debug=False):
    """
    Is this header from an image more than a month old?
    """
    date_obs = ephem.Date(hdr['DATE-OBS'].replace('T',' '))
    today = ephem.now()
    if debug:
        print "date_obs =", date_obs, "today =", today
        print "date_obs - today =", today - date_obs
    return (today - date_obs > 30)

def is_related_to_followup_target(hdr, debug=False):
    """
    Is this header from an image related to existing follow-up targets?
    """
    # If these keywords are missing, it's a flat field or something
    try:
        field = hdr['FIELD_ID']
        subfield = hdr['SUBFIELD']
    except KeyError:
        return False
    if debug:
        print "field, subfield =", field, subfield
        print "sn_locations =", _sn_locations
    return (field, subfield) in _sn_locations

def is_deep_enough_for_ref(hdr, debug=False):
    """
    Is this header from an image deep enough to serve as a good reference?
    """
    # Based on field 3201-10 (SN 2013hx data) ridge line for good 30-sec SUBs:
    #     MAGLIM50 > 20.65 - 0.5*(SEEING - 2.5")
    # Theoretical ridgeline fits the same data just fine
    try:
        bkgsig = hdr['BKGSIG']
        seeing = hdr['SEEING']
        filter = hdr['FILTNAME']
        zpmag  = hdr['ZPMAG']
    except KeyError:
        print "WARNING:  possible problem w/data; keeping file for testing"
        return True
    if 'MAGLIM50' in hdr:
        maglim = hdr['MAGLIM50']
    else:
        maglim = -2.5*np.log10(calc_fluxlim_img(bkgsig, CL=0.50)) + zpmag
        hdr.update('MAGLIM50', maglim, '95% upper limit flux from BKGSIG')
    exptime = 30.0
    nsbsig = np.sqrt(Imager.nsb_cps[filter] * exptime)
    ridgeline = calc_maglim50(nsbsig, seeing, exptime, filter) - 0.2
    if debug:
        print ("bkgsig = {0:.2f}, maglim = {1:.2f}, seeing = {2:.2f}".format(
               bkgsig, maglim, seeing))
        print ("filter = {0}, nsbsig = {1:.2f}, ridgeline = {2:.2f}".format(
               filter, nsbsig, ridgeline))
    return maglim > ridgeline

def is_good_new(hdr, debug=False):
    return (is_related_to_followup_target(hdr, debug=debug)
            or not is_old(hdr, debug=debug))

def is_good_ref(hdr, debug=False):
    return is_deep_enough_for_ref(hdr, debug=debug)

def is_good_sub(hdr, debug=False):
    return (is_related_to_followup_target(hdr, debug=debug)
            or not is_old(hdr, debug=debug))

def new_files(fname):
    """
    Which files to curate for NEW with the given FITS image filename
    """
    regex = "(.*Skymapper_.*).fits"
    match = re.search(regex, fname)
    if match is None:
        return [ ]
    rootfname = match.group(1)
    return [rootfname + ext for ext in ('.fits', '_mask.fits.gz')]

def ref_files(fname):
    """
    Which files to curate for REF with the given FITS image filename
    """
    regex = "(.*Skymapper_.*)_(wcs|mask).fits"
    match = re.search(regex, fname)
    if match is None:
        return [ ]
    rootfname = match.group(1)
    return [rootfname + ext for ext in
            ('_wcs.fits.gz', '_wcs.fits.stars', '_mask.fits.gz')]

def sub_files(fname):
    """
    Which files to curate for SUB with the given FITS image filename
    NB FRAGILE -- only works for files that exist on disk
    """
    # Have everything match the obs time
    regex = "(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})"
    match = re.search(regex, fname)
    if match is None:
        return [ ]
    obstag = match.group(1)
    subpath = os.path.dirname(fname)
    # Now grab all the associated .fits.gz files
    return glob.glob(os.path.join(subpath, '*{0}*.fits.gz'.format(obstag)))


def main():
    # Just search chip 17 on all the images (to make globs go faster)
    # and use that as an indicator for all 32 amps.
    new_glob = os.path.join(CPP.new, '*', '*', '*', '*.fits')
    ref_glob = os.path.join(CPP.ref, '*', '*', '*', '*_wcs.fits.gz')
    sub_glob = os.path.join(CPP.sub, '*', '*', '*', '*.fits.cands')
    for select_func, CPP_glob, select_files in [
            (is_good_new, new_glob, new_files),
            (is_good_ref, ref_glob, ref_files),
            (is_good_sub, sub_glob, sub_files),
    ]:
        for fname in glob.glob(CPP_glob):
            hdr = pyfits.getheader(fname)
            seeing, maglim = None, None
            keepme = select_func(hdr)
            if 'SEEING' in hdr:
                seeing = "{0:.2f}".format(hdr['SEEING']/2.0)
            if 'MAGLIM50' in hdr:
                maglim = "{0:.2f}".format(hdr['MAGLIM50'])
            print fname, seeing, maglim, keepme
            if not keepme:
                for delfn in select_files(fname):
                    cmd = "rm {0}".format(delfn)
                    print cmd
                    # os.unlink(delfn)
            sys.stdout.flush()

if __name__ == "__main__":
    main()
