#! /usr/bin/env python

import argparse
from Utils.DetectionMerger import DetectionMerger
from STAP_comm import print_cmd_line

def xref(mergedfname, candsfnames=[ ], lcfnames=[ ], skybot_cachefname=None, sdss_cachefname=None,
         replace=False):
    """Cross-identify detections in different filters

    Cross-identify detections in the one or more .fits.cands files with
    pre-existing detections (if any) in mergedfname.  Write out new Real
    detections in mergedfname with their sources identified.
    If N_coinc (default 2) or more detections are found at the same RA and DEC,
    give that candidate a name and store it in a 'NAMExxxx' keyword in the
    primary FITS header, with 'xxxx' being the internal integer ID used inside
    the table to identify candidates.
    Finally, write out light curves for all NAMEd vettable candidates.

    candfiles: one or more RealBogus .fits.cands files
    mergedfname: a FITS table of SkymapperDetections corresponding to
        the Real detections in multiple filters in that sky area
    N_coinc: number of times a detection must be re-detected to count as
             a source we want to bother vetting
    """

    print_cmd_line("STAP_xref.py", mergedfname, candsfnames=candsfnames,
                   lcfnames=lcfnames, skybot_cachefname=skybot_cachefname,sdss_cachefname=sdss_cachefname,
                   replace=replace)

    # First merge all available detections into the xref file.
    merger = DetectionMerger(mergedfname, verbose=True)
    for fn in candsfnames:  merger.update_predetections(fn, replace=replace)
    merger.register_xids(write_lc=(lcfnames == None))
    
    # The "vetbot" part:  Run some common vetting tasks.
    merger.vetbot(skybot_cachefname=skybot_cachefname,sdss_cachefname=sdss_cachefname)

    # FY - run vetbot first, then save the results
    merger.write_fits_file(mergedfname)

    # If there were any light curve files given, compile light curves.
    merger.compile_lightcurves(lcfnames)


    # Cross-check asteroid completeness
    # merger.compile_roidxchcke(candsfnames + lcfnames)

    # Return the list of what we detected.  This includes all pre-existing
    # SkyMapper transients -- we'll want to let the user know that these
    # now have new light curve data.
    return merger.name_registry

def main():
    parser = argparse.ArgumentParser(description="Detection merging!")
    parser.add_argument("mergedfname",
                         help="FITS table filename(s) for multi-file cross-IDs")
    parser.add_argument("--candsfnames", nargs="*",
                         help='SUB candidates FITS table filename(s)')
    parser.add_argument("--lcfnames", nargs="*",
                         help='FITS table filename(s) containing object photometry')
    parser.add_argument("--skybot_cachefname",
                         help='ASCII table filename(s) containing Skybot query')
    parser.add_argument("--sdss_cachefname",
                         help='ASCII table filename(s) containing SDSS VizieR query')
    parser.add_argument("--replace", action="store_true", default=False,
                         help='Clobber old detections in survey history file(s)')
    args = parser.parse_args()

    xref(args.mergedfname, args.candsfnames, args.lcfnames, args.skybot_cachefname, args.sdss_cachefname,
         args.replace)

if __name__ == "__main__":
    main()
