#!/usr/bin/env python

"""
RS 2015/05/07:  Pipeline Disk Curation

This infrastructure is written to keep our disks relatively clean by
automatically identifying files that need to be deleted.  Eventually,
we should be able to run this as its own workflow, possibly during the
day as the flip side of STAP_ingest_from_markle.py --production

Design thoughts:

    Use pyfits to read headers
"""

import os
import nose
import ephem
import pyfits
from Utils.Constants import PipelinePath as CPP
from Utils.CurateDisk import is_good_new, is_good_ref, is_good_sub
from Utils.CurateDisk import new_files, ref_files, sub_files


# Make the test harness think it is 7 February 2014
ephem.now = lambda: ephem.Date("2014-02-07")


# Make some mock headers for testing. The implementation doesn't need real
# pyfits header objects right now, just something that acts like a hash with
# the right keywords in it.  This makes our testing more flexible too, since
# it isn't necessarily tied to specific outcomes for data, although the
# values below were originally derived from actual exposures.

_test_fits_mockkeys = [ 'DATE-OBS', 'FIELD_ID', 'SUBFIELD', 'FILTNAME',
                        'SEEING',   'BKGSIG',   'ZPMAG',    'MAGLIM50', ]
_test_fits_mockvals = [
    # Skymapper_893387295_00000_2014-01-31T22:05:09, from testdata/ref
    ('2014-01-31T11:06:10', 3201, 10, 'g', 9.02, 11.1, 25.77, 19.40),
    ('2014-01-31T11:06:10', 3201, 17, 'g', 9.09, 10.2, 25.76, 19.48),
    # Skymapper_893389330_00000_2013-09-11T23:54:45, from testdata/ref
    ('2013-09-11T13:56:20', 3201, 10, 'g', 4.86, 12.5, 27.97, 21.47),
    ('2013-09-11T13:56:20', 3201, 17, 'g', 4.92, 13.1, 28.03, 21.47),
    # Skymapper_893421120_00000_2013-10-02T00:32:42, from testdata/sub
    ('2013-10-01T14:34:14', 3201, 10, 'g', 6.16, 12.6, 27.60, 21.09),
    ('2013-10-01T14:34:14', 3201, 17, 'g', 6.34, 13.5, 27.64, 21.06),
    # Skymapper_893436573_00000_2013-11-16T02:28:04, from testdata/sub
    ('2013-11-15T15:29:07', 3201, 10, 'g', 5.21, 30.6, 26.93, 19.45),
    ('2013-11-15T15:29:07', 3201, 17, 'g', 5.09, 32.6, 27.00, 19.45),
]

# Construct a list of FITS-header-like dicts using the above data
_test_fits_mockhdrs = [
    dict([(key, val) for key, val in zip(_test_fits_mockkeys, valrow)])
    for valrow in _test_fits_mockvals
]


# Use of nose's test generator

def file_select_harness(func, fn, result):
    assert func(fn, debug=True) == result

def test_isgoodnew():
    # NEW:  less than 1 month old or related to existing SNe
    test_results = [1, 1, 1, 0, 1, 0, 1, 0]
    for hdr, result in zip(_test_fits_mockhdrs, test_results):
        yield file_select_harness, is_good_new, hdr, bool(result)

def test_isgoodref():
    # REF:  passes depth cuts (above seeing ridge line)
    test_results = [0, 0, 1, 1, 1, 1, 0, 0]
    for hdr, result in zip(_test_fits_mockhdrs, test_results):
        yield file_select_harness, is_good_ref, hdr, bool(result)

def test_isgoodsub():
    # SUB:  less than 1 month old or related to existing SNe
    test_results = [1, 1, 1, 0, 1, 0, 1, 0]
    for hdr, result in zip(_test_fits_mockhdrs, test_results):
        yield file_select_harness, is_good_sub, hdr, bool(result)


# Testing for which files are associated with what kind of objects.
# The following are collections of files with absolute pathnames associated
# with different kinds of pipeline entities (NEW, REF, SUB).  The first file
# in each collection is always the primary image file, while the others are
# associated files to be curated along with that main file..

_newfiles = [os.path.join(CPP.new, '3201', 'g', '32', fn) for fn in (
    'Skymapper_893387338_00000_2014-01-30T22:34:26_32.fits',
    'Skymapper_893387338_00000_2014-01-30T22:34:26_32_mask.fits.gz'
)]

_reffiles = [os.path.join(CPP.ref, '3201', 'g', '32', fn) for fn in (
    'Skymapper_893389330_00000_2013-09-11T23:54:45_32_wcs.fits.gz',
    'Skymapper_893389330_00000_2013-09-11T23:54:45_32_mask.fits.gz',
    'Skymapper_893389330_00000_2013-09-11T23:54:45_32_wcs.fits.stars'
)]

_subfiles = [os.path.join(CPP.sub, '3201', 'g', '32', fn) for fn in (
    'sub20150429_115440_3201_g_2013-10-20T01:55:33_32.fits.gz',
    'sub20150429_115440_3201_g_2013-10-20T01:55:33_32_mask.fits.gz',
    'Skymapper_893421199_00000_2013-10-20T01:55:33_32_wcs.fits.gz',
    'Skymapper_893421199_00000_2013-10-20T01:55:33_32_mask.fits.gz',
    'Skymapper_893389330_00000_2013-09-11T23:54:45_32_wcs_2013-10-20T01:55:33_remap.fits.gz',
    'Skymapper_893389330_00000_2013-09-11T23:54:45_32_wcs_2013-10-20T01:55:33_remap_mask.fits.gz'
)]

def file_assoc_harness(func, collect):
    assert set(func(collect[0])) == set(collect)

def test_file_collections():
    for func, collect in zip(
            (new_files, ref_files, sub_files),
            (_newfiles, _reffiles, _subfiles)):
        print "collect =", collect
        yield file_assoc_harness, func, collect
