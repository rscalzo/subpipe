#!/usr/bin/env python

# ============================================================================
# RS 2012/02/15:  Creating Multi-Filter Detection Files
# ----------------------------------------------------------------------------
# This code takes as input several *.cands (RealBogus) files as input,
# cross-references them and stores the result in another file.
# ============================================================================

import re
import os
import sys
import datetime
import argparse
import pyfits
import numpy as np
import ephem
from .Record import FITSRecord, KeywordRecord
from .Catalog import CatalogEntry, SExtractorDetection, \
        SkybotObject, APASSObject, match_coords, match_coords_v2,SDSSObject
from .RealBogus import RealBogus
from .TrackableException import TrackableException
from .Photometry import calc_fluxlim_img, sexzpidx
from .Constants import get_obsseq, FilenamesXsient, Filenames, SMfield_lookup


# What do I want this code to be able to do?
# -- Search for spatially coincident detections among various .cand files
#    in multiple filters and on multiple days
# -- Dump a final list of coincidences among filters for each (field, CCD).
# -- For objects which have more than N detections, register them in some
#    central SQL database.
# -- Compute a crude (differential) historical light curve for each source
#    in the SQL database.
# -- Take care of any Python plotting necessary.
# In short, this could be the engine for the historical light curve builder
# which I plan to parallelize and operate in the normal pipeline.


# ============================================================================
#     Definition of class DetectionMerger, where most of the action is
# ============================================================================


class SkymapperDetection(FITSRecord, KeywordRecord, CatalogEntry):
    """
    Represents a FITS record used in the SkyMapper multi-filter catalogs.
    Used for things like light curves, etc.
    """

    _fits_fields = \
    [
        [ "id",      "XFILT_ID",     "1J", "",      "I4",    0    ],
        [ "runtag",  "RUNTAG",      "15A", "",     "A15",    ""   ],
        [ "filter",  "FILTER",       "5A", "",      "A5",    ""   ],
        [ "field",   "SM_FIELD",     "1I", "",      "I4",    ""   ],
        [ "obsseq",  "OBSSEQ",      "19A", "",     "A19",    ""   ],
        [ "subfield", "SUBFIELD",    "1I", "",      "I2",    ""   ],
        [ "ccd",     "SM_CCD",       "1I", "",      "I2",    ""   ],
        [ "x",       "X_IMAGE",      "1E", "pixel", "F7.2",  0.0  ],
        [ "y",       "Y_IMAGE",      "1E", "pixel", "F7.2",  0.0  ],
        [ "jd",      "JD_OBS",       "1D", "day",   "F11.3", 0.0  ],
        [ "ra",      "ALPHA_J2000",  "1D", "deg",   "F10.6", None ],
        [ "dec",     "DELTA_J2000",  "1D", "deg",   "F10.6", None ],
        [ "exptime", "EXPTIME",      "1E", "sec",   "F8.3",  None ],
        [ "flux",    "FLUX_AP4",     "1E", "count", "G10.5", 0.0  ],
        [ "fluxerr", "FLUX_AP4_ERR", "1E", "count", "G10.5", 0.0  ],
        [ "flag",    "FLAGS",        "1I", "",      "I3",    0    ],
        [ "rbscore", "RBSCORE",      "1E", "",      "F5.3",  0.0  ],
        [ "apcor",   "APCORR04",     "1E", "mag",   "F5.3",  0.00 ],
        [ "apcerr",  "APCERR04",     "1E", "mag",   "F5.3", 99.99 ],
        [ "zp",      "ZPMAG",        "1E", "mag",   "F5.3",  0.00 ],
        [ "zperr",   "ZPMAGERR",     "1E", "mag",   "F5.3", 99.99 ],
        [ "mag",     "MAG_CAL",      "1E", "mag",   "F5.3",  0.00 ],
        [ "magerr",  "MAG_CAL_ERR",  "1E", "mag",   "F5.3", 99.99 ],
    ]
    _keyword_fields = [ (f[0],f[5]) for f in _fits_fields ]

    def is_upper_limit(self):
        return (self.fluxerr not in [None, np.nan, np.inf]
                and self.fluxerr < 0.0)

    def _init_from_catobj(self, cat, **kwargs):
        # First initialize from keywords
        zpi = sexzpidx
        if isinstance(cat, RealBogus):
            self._init_from_keywords(
                x    = cat.xsub,      y=cat.ysub,
                ra   = cat.ra,        dec=cat.dec,
                flux = cat.f4sub,     fluxerr=cat.df4sub,
                flag = cat.flagsub,   rbscore=cat.rbscore, **kwargs)
        elif isinstance(cat, SExtractorDetection):
            self._init_from_keywords(
                x    = cat.x,         y=cat.y,
                ra   = cat.ra,        dec=cat.dec,
                flux = cat.flux[zpi], fluxerr=cat.fluxerr[zpi],
                flag = cat.flag,      rbscore=100, **kwargs)
        # Selected keywords have meanings:  apply zeropoint.
        # Treat upper limits differently from standard photometry;
        # "fluxerr" and "magerr", if negative, are confidence levels.
        if self.is_upper_limit():
            self.magerr = self.fluxerr
        # Normal photometry points will have positive flux.
        elif self.flux > 0.0 and self.flux/self.fluxerr > 2.0:
            self.mag = -2.5*np.log10(self.flux) + self.zp
            self.magerr = np.sqrt(
                (1.0857*(self.fluxerr/self.flux))**2 + self.zperr**2)
        # If the measured flux is negative, it must be because SExtractor
        # detected some clearly-bogus-thing there, so set to junk values.
        else:
            self.mag, self.magerr = 10.0, -1.0

    def __init__(self, *args, **kwargs):
        """
        Try to init in this order:  from a RealBogus object (1 argument);
        from a FITS table row (2 arguments); keywords as available.
        """
        if len(args) == 1:   self._init_from_catobj(*args,**kwargs)
        if len(args) == 2:   self._init_from_fitsrow(*args)
        if len(kwargs) > 0:  self._init_from_keywords(**kwargs)


class SkymapperTransient(CatalogEntry, KeywordRecord):
    """A small wrapper class to help with cross-IDs."""

    _keyword_fields = [ [ "id",        0 ], [ "name",     "" ],
                        [ "ra",      0.0 ], [ "dec",     0.0 ],
                        [ "field",     0 ], [ "subfield",  0 ],
                        [ "nlcpts",    0 ], [ "last-jd", 0.0 ], ]

    def __init__(self, **kwargs):
        self._init_from_keywords(**kwargs)
        super(SkymapperTransient, self).__init__(self.ra, self.dec)

    def __str__(self):
        return "{0} (id {1}) in field {2}-{3} with {4} obs".format(
                self.name, self.id, self.field, self.subfield, self.nlcpts)


class DetectionMerger(object):
    """Handles cross-IDs of separate detections of a source."""

    version_tag = "v1.0.0a"

    def __init__ (self, predetect_file=None, coinc=2, verbose=False):
        """Initialize the DetectionMerger
        
        predetect_file:  filename of existing SkymapperDetections
        coinc:  number of predetections needed to register a new transient
        """
        self.coinc = coinc
        self.verbose = verbose
        if predetect_file != None and os.path.exists(predetect_file):
            # Read the candidates themselves in from the file.
            self.predetections, self.header = \
                SkymapperDetection.read_fits_file(predetect_file)
            # Sort out the upper limits from the actual predetections.
            self.upperlimits = [s for s in self.predetections if s.id < 0]
            self.predetections = [s for s in self.predetections if s.id >= 0]
            # The names of registered candidates (up to 1e+4 per file)
            # are stored in the header as keywords of the form "NAMExxxx".
            # These are mirrored in an internal "name registry" which
            # looks up the internal ID associated with a given source name.
            self.name_registry = { }
            for (key, value) in self.header.items():
                namematch = re.match("NAME(\d{4})", key)
                if namematch:
                    field = self.predetections[0].field
                    subfield = self.predetections[0].subfield
                    id = int(namematch.group(1))
                    self.name_registry[id] = { "name": value, 'id': id }
        else:
            self.predetections = [ ]
            self.upperlimits = [ ]
            self.name_registry = { }
            self.header = None

    def write_fits_file(self, predetect_file):
        SkymapperDetection.write_fits_file(predetect_file,
                self.predetections + self.upperlimits, header=self.header)

    @staticmethod
    def candidate_name(ra, dec):
        """
        Standardized naming convention for SkyMapper candidates.
        Use (RA,DEC) converted to sexagesimal to the nearest arcsec.
        """
        # RS 2014/02/05:  Fiddled with it to get all the zeroes right.
        # ephem doesn't give the most standardized output by default.  :/
        coo = ephem.Equatorial(ra*ephem.pi/180.0, abs(dec*ephem.pi/180.0))
        ra_str  = re.sub('[:.]','',str(coo.ra))
        dec_str = re.sub('[:.]','',str(coo.dec))
        if ra < 150.0:       ra_str = "0" + ra_str
        if abs(dec) < 10.0:  dec_str = "0" + dec_str
        if dec >= 0.0:  dec_str = "+" + dec_str
        if dec < 0.0:   dec_str = "-" + dec_str
        return "SMTJ{0}{1}".format(ra_str,dec_str)

    @staticmethod
    def convert_detections_file(candsfname, lightcurve=False):
        """
        Given a RealBogus file, produces a list of SkymapperDetections,
        plus an upper limit computed from the header information.
        """
        # Check the header to see what kind it is
        newheader = pyfits.getheader(candsfname, 1)
        if 'CLASSNAM' in newheader:
            classname = newheader['CLASSNAM']
        else:
            raise TrackableException("FITS file incompatible with "
                                     "this version of DetectionMerger")

        if classname == 'RealBogus':
            newdetections, newheader = RealBogus.read_fits_file(candsfname)
            realdetections = [rb for rb in newdetections
                              if rb.rbscore > RealBogus.rbscore_thresh]
            if not lightcurve:  newdetections = realdetections
        elif classname == 'SExtractorDetection':
            newdetections, newheader = \
                    SExtractorDetection.read_fits_file(candsfname)

        # RS 2014/07/02:  Apparently we're tracking the obsseq now.  Fixing a
        # bug by filling this up front, instead of in an if/then below where
        # it may be missed.
        obsseq = get_obsseq(candsfname)

        # Extract some key information from the FITS header.  Not all of
        # these quantities are currently being filled, so we'll need to start
        # filling them at least with example values (probably in STAP_SEx.py)
        # and then getting the real keyword names from Stefan later.
        field, filter, ccd, exptime, jd = 0, 0, 0, 30.0, 2451545.0
        zp, zperr = 24.0, 5.0
        subfield = ccd
        runtag = ""
        bkgsig = newheader["BKGSIG"]
        # RS 2014/11/10:  HACK to work in catalog space
        if runtag not in newheader:
            newheader.update('RUNTAG', 'catspace', '')
        try:
            kw_missing = [kw for kw in
                          ('RUNTAG',  'FILTNAME', 'FIELD_ID', 'SUBFIELD',
                           'IMAGEID', 'EXPTIME',  'ZPMAG',    'ZPMAGERR')
                          if kw not in newheader]
            if len(kw_missing) > 0:
                raise KeyError("{0} missing light curve keywords {1}".format
                               (candsfname, str(kw_missing)))
            runtag   = newheader["RUNTAG"]
            filter   = newheader["FILTNAME"]
            field    = newheader["FIELD_ID"]
            subfield = newheader["SUBFIELD"]
            ccd      = newheader["IMAGEID"]
            exptime  = newheader["EXPTIME"]
            zp       = newheader["ZPMAG"]
            zperr    = newheader["ZPMAGERR"]
            date = re.findall('\d+', newheader["DATE-OBS"])
            jd = ephem.Date(tuple([int(i) for i in date])) + 2415020
            converted_detections = \
                [SkymapperDetection(nd, subfield=subfield, obsseq=obsseq,
                                jd=jd, field=field, filter=filter, ccd=ccd,
                                zp=zp, zperr=zperr, exptime=exptime,
                                runtag=runtag) for nd in newdetections]
        except Exception as e:  # KeyError as e:
            print "DetectionMerger.convert_detections_file:  ", e
            converted_detections = [ ]
            
        # Record a limiting flux for this subtraction, if we haven't yet
        # done so (check for existing records within 5 seconds).  Empirically,
        # I find the 50% completeness flux is around 6.0 (background) sigma,
        # so the 95% CL upper limit is ~4.5 sigma in the 4-pixel aperture.
        # This is a rough limit for screening purposes; for science-grade
        # upper limits I recommend going back to the original pixels.
        # Present the upper limit as a SkymapperDetection itself.
        flux_ulim = calc_fluxlim_img(bkgsig, aperture=6.0, CL=0.95)
        mag_lim = -2.5*np.log10(flux_ulim) + zp
        ulimit = SkymapperDetection(id=-1, runtag=runtag, filter=filter,obsseq=obsseq,
                    field=field, subfield=subfield, ccd=ccd, x=0.0, y=0.0,
                    jd=jd, ra=0.0, dec=0.0, exptime=exptime,
                    flux=flux_ulim, fluxerr=-0.95, flag=0, rbscore=0.0,
                    mag=mag_lim, magerr= -0.95, zp=zp, zperr=zperr)

        # Return the results!
        return converted_detections, ulimit
    
    def update_predetections(self, candsfname, replace=False):
        """
        Merges a list of RealBoguses into a list of SkymapperDetections
        for this {field,CCD} combination.
        """
        # Read in the new detections and cache the most recent ones.
        # This will be helpful for vetting purposes (see below).
        newdetections, ulimit = DetectionMerger.convert_detections_file(
                candsfname, lightcurve=False)
        if len(newdetections) < 1:  return
        self.newdetections = newdetections

        # See which new detections match which pre-detections
        # (one index for each pre, each pointing to a matched new)
        if self.verbose:
            print "DetectionMerger.update_predetections: ",
            print "Matching", len(newdetections), "detections in", candsfname
        backmatchidx, dist = match_coords(
                self.predetections, newdetections, tol=1.0)
        # Compile a list of lists of pre-matches for each new detection
        # RS 2013/10/04:  For clarity, backmatchidx is a list with length
        # len(self.predetections) containing the indices in newdetections
        # of new detections matched to each element of self.predetections.
        # Then fwdmatchidx is a list of lists with length len(newdetections),
        # each location containing possibly multiple predetections for each
        # element of newdetections.
        fwdmatchidx = [[] for nd in newdetections]
        for i in range(len(backmatchidx)):
            if backmatchidx[i] != None:
                fwdmatchidx[backmatchidx[i]].append(i)
        # Now go through and assign pre-detected IDs to new detections.
        running_idx = len(set([obj.id for obj in self.predetections]))
        for i in range(len(fwdmatchidx)):
            redundant = False
            sd = newdetections[i]
            # If we found no predetection(s) for this new detection, assign
            # it to a new source.
            if len(fwdmatchidx[i]) == 0:
                if self.verbose:
                    print "-- No pre-existing objects at location of new "\
                          "detection", i,"; adding new xref ID", running_idx
                sd.id = running_idx
                running_idx += 1
            # If there are pre-detections, do some bookkeeping...
            else:
                # For the new detection, use the source ID of an arbitrary
                # matched pre-detection; they should all be the same.
                # (We might want to check whether this is really true!)
                idset = set([self.predetections[j].id for j in fwdmatchidx[i]])
                if len(idset) == 0:
                    print "FATAL:  no forward matches for new detection", i
                    print "We should never get here!  Why are we even here?!"
                    raise TrackableException("Bad pre-detection cross-match")
                elif len(idset) > 1:
                    print "WARNING:  ambiguous detection --",
                    print "new detection", i, "matches IDs", idset
                # Check whether any pre-detections have the same observation
                # time as this one, to within 5 seconds.  This is relevant in
                # case we're re-running subtractions for this area; we don't
                # want to keep re-finding the same stuff.
                sd.id = self.predetections[fwdmatchidx[i][0]].id
                for j in fwdmatchidx[i]:
                    if abs(sd.jd - self.predetections[j].jd) < 5*ephem.second:
                        if self.verbose:
                            print "-- Found old detection at location of "\
                                  "new detection", i, ";",
                        redundant = True
                        if replace:
                            if self.verbose:  print "replacing with new detection"
                            self.predetections[j] = sd
                        else:
                            if self.verbose:  print "ignoring new detection"
                        break
            if not redundant:
                self.predetections.append(sd)

    def compile_lightcurves(self, candsfnamelist, sublc=True):
        """
        Compiles historical light curves, given a list of filenames containing
        RealBogus objects, i.e., the whole survey history for this {field,CCD}.
        For now, we're just using the internal registry of candidates, since
        this is for vetting.  For science we'll have something better.
        """

        self.predetections, self.upperlimits, self.boundslist = [ ], [ ], [ ]
        for candsfname in set(candsfnamelist):
            # Grab the minra, mindec, maxra, maxdec
            hdr = pyfits.getheader(candsfname)
            if 'MINRA' in hdr:
                minra, mindec = hdr['MINRA'], hdr['MINDEC']
                maxra, maxdec = hdr['MAXRA'], hdr['MAXDEC']
                coomin = ephem.Equatorial(minra, mindec)
                coomax = ephem.Equatorial(maxra, maxdec)
            else:
                raise TrackableException("FITS file incompatible with "
                                         "this version of DetectionMerger")
            # Grab the actual detections
            newdetections, ulimit = DetectionMerger.convert_detections_file(
                    candsfname, lightcurve=True)
            # Stop and make sure they have all the necessary info
            if len(newdetections) < 1:  continue
            self.predetections += newdetections
            # This is a bit of a HACK, but it should work as long as we have
            # access to all the light curve files.  What would be much better
            # is to make an actual SkymapperFluxLimit class.
            rad2deg = 180/ephem.pi
            self.upperlimits.append(ulimit)
            self.boundslist.append((coomin.ra*rad2deg, coomin.dec*rad2deg,
                                    coomax.ra*rad2deg, coomax.dec*rad2deg))

        # Turn the name registry hash entries into CatalogEntries for matching
        sources = [SkymapperTransient(**kw)
                   for kw in self.name_registry.values()]

        # Since we already know where the sources are, all we have to do is
        # match them -- we don't need to iterate.
        if self.verbose:
            print "DetectionMerger.compile_lightcurves: ",
            print "Matching {0} detections to {1} pre-detected events".format(
                len(self.predetections), len(sources))
        backmatchidx, dist = match_coords(
                self.predetections, sources, tol=1.0)
        # Compile detections corresponding to each source.
        # RS 2013/10/04:  For clarity, backmatchidx is a list with length
        # len(self.predetections) containing the indices in sources 
        # of transient events matched to each element of self.predetections.
        # Then fwdmatchidx is a list of lists with length len(sources), with
        # each location containing possibly multiple predetections for each
        # element of sources.  NB that if a source has been registered,
        # it must have a light curve, so all elements of self.lightcurves
        # should have at least one light curve point after this!
        self.lightcurves = [[] for s in sources]
        for i in range(len(backmatchidx)):
            if backmatchidx[i] != None:
                self.lightcurves[backmatchidx[i]].append(self.predetections[i])
        # self.lightcurves now contains light curves including all Real and
        # Bogus detections of each source.  Now add in the non-detections.
        for i in range(len(sources)):
            src, lc = sources[i], self.lightcurves[i]
            # RS 2014/02/06:  Sometimes we get a weird ValueError indicating
            # that the light curve for an existing candidate is empty.
            # This can happen during re-runs when a candidate that exists in
            # the registry isn't in the part of the survey history we're
            # currently examining.  What a mess.
            if len(lc) < 1:
                print "-- WARNING:  ", sources[i].name,
                print "not matched ANYWHERE in survey history -- skipping"
                del self.name_registry[src.id]
                sys.stdout.flush()
                continue
            # Go through the upper limits and see which ones correspond to
            # exposures which *don't* have a corresponding positive detection.
            # Also insist that an exposure be within the image bounds!
            for ul, bounds in zip(self.upperlimits, self.boundslist):
                ul_mindt = np.abs([ul.jd - sd.jd for sd in lc]).min()
                minra, mindec, maxra, maxdec = bounds
                if (ul_mindt > 5*ephem.second
                        and src.ra > minra and src.dec > mindec
                        and src.ra < maxra and src.dec < maxdec):
                    ul.x, ul.y, ul.ra, ul.dec = 0.0, 0.0, src.ra, src.dec
                    lc.append(ul)
            # Sort in time order and establish uniqueness.
            lc = sorted(lc, key=lambda sd: sd.jd)
            lc = [lc[ii] for ii in range(len(lc))
                  if ii == 0 or abs(lc[ii].jd - lc[ii-1].jd) > 5*ephem.second]
            # update name_registry with all observations
            islim = [sd.fluxerr <=0 for sd in lc]
            self.name_registry[src.id].update(
                    { "jd": [sd.jd for sd in lc], "islim": islim })
            lcheader = pyfits.PrimaryHDU().header
            lcheader.update("CANDNAME", src.name,
                            "Unique name of candidate")
            lcheader.update("FIELD_ID", src.field,
                            "SkyMapper main survey ID for this (RA,DEC)")
            lcheader.update("SUBFIELD", src.subfield,
                            "Subfield (SkyMapper CCD for ROTSKYPA=0)")
            lcfnames = FilenamesXsient(src.name, src.field, src.subfield)
            if sublc:
                destfn = lcfnames.basesublcfname
            else:
                destfn = lcfnames.basenewlcfname
            SkymapperDetection.write_fits_file(destfn, lc, header=lcheader)

    def register_xids(self, write_lc=True):
        """
        Goes through the multi-detection list to see if we got anything fun.
        """

        idset = set([obj.id for obj in self.predetections])
        self.lightcurves = [[ ] for id in idset]
        for sd in self.predetections:
            self.lightcurves[sd.id].append(sd)

        # Did we detect anything more than N times?  Consider all filters all
        # together for now, we'll consider specific filters later if need be.
        for id in range(len(self.lightcurves)):
            lc = self.lightcurves[id]
            if len(lc) >= self.coinc or id in self.name_registry:
                # Create a name for the thing if it doesn't already have one.
                ra  = np.mean([sd.ra  for sd in lc if sd.fluxerr > 0 and
                               sd.rbscore > RealBogus.rbscore_thresh])
                dec = np.mean([sd.dec for sd in lc if sd.fluxerr > 0 and 
                               sd.rbscore > RealBogus.rbscore_thresh])
                field, subfield = lc[0].field, lc[0].subfield
                if id in self.name_registry:
                    name = self.name_registry[id]["name"]
                else:
                    name = self.candidate_name(ra, dec)
                    self.name_registry[id] = { "name": name, "autotype": "?" }
                    if self.verbose:
                        print "Adding new transient", name, "with id", id

                # Add other information to the registry; this will be passed
                # back through the layers eventually to be written to disk.
                islim = [sd.fluxerr <=0 for sd in lc]
                self.name_registry[id].update({ "id": id, "ra": ra, "dec": dec,
                    "field": field, "subfield": subfield,
                    "jd": [sd.jd for sd in lc], "islim": islim })
                # Also write the new object ID to the xref FITS header,
                # with an ID matching the order in which it was discovered.
                if self.header == None:
                    self.header = pyfits.PrimaryHDU().header
                namekw = "NAME{0:04d}".format(int(id))
                self.header.update(namekw, name,
                                   "Name of source with index {0}".format(id))

                # Create a light curve file for the thing, and store it in a
                # directory the name of which we can construct on demand.
                # In practice, when updating light curve information we will
                # simply clobber the old file with the new (complete) data.
                if not write_lc:  continue
                lcheader = pyfits.PrimaryHDU().header
                lcheader.update("CANDNAME", name,
                                "Unique name of candidate")
                lcheader.update("FIELD_ID", field,
                                "SkyMapper main survey ID for this (RA,DEC)")
                lcheader.update("SUBFIELD", subfield,
                                "Subfield (SkyMapper CCD for ROTSKYPA=0)")
                lcfnames = FilenamesXsient(name, field, subfield)
                SkymapperDetection.write_fits_file(
                        lcfnames.basesublcfname, lc, header=lcheader)

    def vetbot(self, skybot_cachefname=None, sdss_cachefname=None):
        """Runs some common vetting tasks so humans don't have to

        This function does the following for each candidate:
        -- Cross-matches against obvious catalogs such as 2MASS, APASS DR6,
           and Skybot to look for asteroids, variable stars, host galaxies
           or anything else interesting near the candidate position
        -- Calculates the candidate's proper motion in case it might be a roid
        -- Looks at the candidate's historical light curve to see if it is
           rising, falling, periodic, or anything else interesting
        It then notes everything it found in the journal for that candidate.
        """

        def format_comment(name, type_str, comment_str):
            # 09. August 2012 09:04PM: user rscalzo thinks it is Junk.
            # Comment:  This candidate sucks!  What is the deal with it?!
            datestr = datetime.datetime.now().strftime("%d. %B %Y %I:%M %p")
            return "#\n{0}:  VetBot {1} thinks {2} is a {3}.\nComment: {4}\n"\
                    .format(datestr, self.version_tag,
                            name, type_str, comment_str)

        def match_objtype(sources,objects,tol,objtype,objsource,tag,update=True):
            matchidx, dist = match_coords(sources, objects, tol)
            for i in range(len(matchidx)):
                if matchidx[i] != None:
                    if objects[matchidx[i]].__dict__.has_key('name'):
                        objname=objects[matchidx[i]].name
                    else:
                        objname=objects[matchidx[i]]
                    if update or not self.name_registry[sources[i].id].has_key("autotype"):
                        comment = format_comment(sources[i].name, tag,
                                                 "{0:.1f} arcsec from {1} {2} {3}".format
                                                 (dist[i], objsource, objtype, objname))
                        print comment
                        self.name_registry[sources[i].id].update(
                            { "autotype": tag,
                              "match": objects[matchidx[i]],
                              "comment":  comment })
            print "VetBot:  matched {0}/{1} {2} with {3} candidates"\
                .format(len([m for m in matchidx if m != None]),
                        len(objects), objtype, len(sources))
            

        # Turn the name registry hash entries into CatalogEntries for matching
        # (increasingly I'm thinking they should just be in that form already)
        sources = [SkymapperTransient(**kw)
                   for kw in self.name_registry.values()]
        if len(sources) == 0:  return

        # RS 2014/02/01:  Annoyingly, a couple of weeks of screwed-up
        # astrometry in our headers has leaked bogus candidates into our
        # survey history.  The right thing to do would be to regenerate it
        # from the beginning, but that's too hard a problem for me to cope
        # with at 2 in the morning.  So just use the field RA, DEC.
        chcke_field_id = set([obj.field for obj in self.predetections])
        if len(chcke_field_id) < 1:
            print "DetectionMerger.vetbot():  WHOA PROBLEM --",
            print "no predetections!  Why are you vetting nothing?!"
            return
        elif len(chcke_field_id) > 1:
            print "DetectionMerger.vetbot():  WHOA PROBLEM --",
            print "predetections come from multiple fields:", chcke_field_id
            return
        chcke_field = SMfield_lookup(id=list(chcke_field_id)[0])

        # Get ready to match stars to the catalog
        ralist, declist = [s.ra for s in sources], [s.dec for s in sources]
        ra_min, dec_min = np.min(ralist), np.min(declist)
        ra_max, dec_max = np.max(ralist), np.max(declist)
        ra_mid, dec_mid = (ra_max + ra_min)/2, (dec_max + dec_min)/2
        cD = np.cos(np.pi*0.5*(dec_max - dec_min)/180.0)
        R = 1.5*np.max([cD*(ra_max - ra_min)/2, (dec_max - dec_min)/2])
        # RS 2014/02/01:  I used to think the above code was clever until it
        # started querying 1500 deg^2 areas in the APASS catalog.  If that
        # threatens to happen, just use the field RA and DEC; it's inefficient
        # but not horribly so.
        if R > 3.5:
            ra_mid, dec_mid, R = chcke_field.ra, chcke_field.dec, 3.5
            print "WARNING:  something's fishy about this field"

        # generalized catalog matching:
        catalogs = \
        [
            # class tol type comment
            # [ SkybotObject, 1.0, "Roid", "Skybot asteroid" ],
            # [ APASSObject, 1.0, "Star", "APASS DR6 star" ],
            # [ TwoMassPSCObject, 1.0, "Star", "2MASS PSC star" ],
            # [ TwoMassXSCObject, 1.0, "Cand", "2MASS XSC galaxy" ],
        ]

        # First!  Pull up the cached SkyBot info and cross-match.
        if skybot_cachefname is not None:
            if not os.path.exists(skybot_cachefname):
                print "VetBot:  can't find Skybot cache", skybot_cachefname
            else:
                print "VetBot:  Reading Skybot cache", skybot_cachefname
                roids = SkybotObject.read_ascii_file(skybot_cachefname)
                match_objtype(sources,roids,1.0,"Roid","Skybot","Roid")

        # Next:  APASS DR6, from Simon's local thing.
        print "VetBot:  Querying APASS DR6 for stars"
        apass_stars = APASSObject.pgquery(ra_mid, dec_mid, R)
        match_objtype(sources,apass_stars,1.0,"Star","APASS DR6","Star")
            
        # Next: Check SDSS from VizieR.
        # Check seperately for star, qso and galaxy
        if sdss_cachefname is not None:
            if not os.path.exists(sdss_cachefname):
                print "VetBot:  can't find SDSS cache", sdss_cachefname
            else:
                print "VetBot:  Reading SDSS cache", sdss_cachefname
                sdss = SDSSObject.read_ascii_file(sdss_cachefname)
                if len(sdss)==0:
                    print "VetBot: no SDSS sources available"
                    return
                #tag spec star
                #not phot star, because:
                #1. phot star includes many QSO
                #2. phot star/qso (point source) has low ~1% false rate in
                #spec sample, but that is definitely not representitive of the whole sample
                sdss_star = [obj for obj in sdss if 'STAR' in obj.sptype] 
                #but, too many stars remain in candidates, so try to eliminate these
                sdss_psc =[obj for obj in sdss if (obj.rmag>12 and obj.rmag <19 and len(obj.sptype)==0 and obj.phtype==6)] 
                
                #tag spec qso
                sdss_qso = [obj for obj in sdss if 'QSO' in obj.sptype]
                
                #tag all possible galaxy, false rate ~3% based on spec sample
                sdss_xsc = [obj for obj in sdss if 'GALAXY' in obj.sptype or
                            (len(obj.sptype)==0 and obj.phtype==3)]
                
                match_objtype(sources,sdss_star,1.0,"Star","SDSS DR9","Star")
                match_objtype(sources,sdss_qso,1.0,"QSO","SDSS DR9","QSO")
                match_objtype(sources,sdss_xsc,5.0,"Galaxy","SDSS DR9","Cand",update=False)
                match_objtype(sources,sdss_psc,1.0,"PointSource","SDSS DR9","Star",update=False)

 
    def compile_roidxchcke(self, candsfnamelist):
        """Cross-checks RealBogus completeness/efficiency against asteroids

        Does what it says.  More details soon.
        """
        for candsfname in set(candsfnamelist):
            # print "Processing", candsfname
            # Grab the minra, mindec, maxra, maxdec
            hdr = pyfits.getheader(candsfname)
            if 'MINRA' in hdr:
                minra, mindec = hdr['MINRA'], hdr['MINDEC']
                maxra, maxdec = hdr['MAXRA'], hdr['MAXDEC']
                coomin = ephem.Equatorial(minra, mindec)
                coomax = ephem.Equatorial(maxra, maxdec)
            else:
                raise TrackableException("FITS file incompatible with "
                                         "this version of DetectionMerger")
            # Grab the actual detections
            sources, ulimit = DetectionMerger.convert_detections_file(
                    candsfname, lightcurve=True)
            if len(sources) < 1:  continue

            # Convert boundaries to decimal degrees
            rad2deg = 180/ephem.pi
            minra, mindec = coomin.ra*rad2deg, coomin.dec*rad2deg
            maxra, maxdec = coomax.ra*rad2deg, coomax.dec*rad2deg

            # Retrieve SkyBot cache and match against detections
            sub = Filenames(candsfname)
            skybot_cachefname = sub.skybot_cachefname
            if not os.path.exists(skybot_cachefname):
                print "RoidXchcke:  can't find Skybot cache", skybot_cachefname
                continue
            roids = SkybotObject.read_ascii_file(skybot_cachefname)
            matchidx, dist = match_coords(roids, sources, tol=1.0)
            for i in range(len(matchidx)):
                if not (minra < roids[i].ra < maxra
                    and mindec < roids[i].dec < maxdec):
                    continue
                # roid:  name, V magnitude;
                print "RoidXchcke: ",
                print "{0:20s}".format(roids[i].name.replace(" ","_")),
                print "{0:5.2f}".format(roids[i].mV),
                # image:  band, zeropoint, seeing, background;
                mjdobs = ephem.Date(sub.dateobs) + 15020.0
                print "{0:.3f} {1} {2:5.2f} {3:5.2f} {4:5.1f}".format(
                      mjdobs, sub.filter, sub.zpmag, sub.seeing, sub.bkgsig),
                mag_lim = -2.5*np.log10(4.5*3.5*sub.bkgsig) + sub.zpmag
                print "{0:5.2f}".format(mag_lim),
                # candidate:  flux, fluxerr, apsig, rbscore
                if matchidx[i] != None:
                    mycand = sources[matchidx[i]]
                    flux = float(mycand.flux)
                    fluxerr = float(mycand.fluxerr)
                    if fluxerr > 0.0:
                        apsig = flux/fluxerr
                    else:
                        apsig = 0.0
                    rbscore = float(mycand.rbscore)
                else:
                    flux, fluxerr, apsig, rbscore = 0.0, 0.0, 0.0, 0.0
                print "{0:6.1f} {1:5.1f} {2:5.1f} {3:5.1f}".format(
                        flux, fluxerr, apsig, rbscore)
                sys.stdout.flush()
