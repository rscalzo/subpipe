#!/usr/bin/env python

# ============================================================================
# RS 2012/02/13:  Implementation of CatalogEntry class hierarchy
# ----------------------------------------------------------------------------
# This code defines a class hierarchy for objects we get from, or put into,
# catalogs, i.e., things which have RAs and DECs, and possibly magnitudes.
# The "as a function of time" thing will probably be handled by collections
# of these catalog entries, so for example LightCurves or RockOrbits might
# be iterable collections of CatalogEntry objects.
# ============================================================================

import re
import shlex
import subprocess as sp
import ephem
import pyfits
import numpy as np
import datetime
from Record import AsciiRecord, KeywordRecord, FITSRecord, PGRecord, BadInit


class CatalogEntry(object):
    """
    This represents anything that's guaranteed to have J2000 RA and DEC.
    It's agnostic about how these fields are eventually filled, although
    the "init" method gives you a hook to put in simple values.
    """

    def __init__(self, ra, dec):
        """Initialize with RA, DEC in decimal degrees; handle 0 = 360 wrap."""
        try:  self.ra, self.dec = float(ra) % 360, float(dec)
        except:
            raise BadInit ("Bad CatalogEntry input:  ra = {0}, dec = {1}"
                           .format(ra,dec))


def match_coords(list1, list2, tol = 1.5):
    """
    This just goes through and checks which members (ra1,dec1) are close to
    which other members (ra2,dec2), within matching tolerance tol in arcsec.
    list1 and list2 are nominally lists of CatalogEntry instances, but could
    be lists of anything with 'ra' and 'dec' attributes in decimal degrees.
    """
    tol1, tol2 = tol/3600.0, (tol/3600.0)**2
    jmin, match_idx, best_dr2 = 0, [None]*len(list1), [tol2]*len(list1)
    # This next bit is rather tricky because we have to sort the lists,
    # but also remember where we put each item.  Here i, j = sorted indices,
    # while ci[i], ci[j] = corresponding actual indices in original lists.
    ci = np.array([c.ra for c in list1]).argsort()
    cj = np.array([c.ra for c in list2]).argsort()
    for i in range(len(list1)):
        CD = np.cos(list1[ci[i]].dec*np.pi/180.0)
        ramin, ramax = list1[ci[i]].ra - tol1/CD, list1[ci[i]].ra + tol1/CD
        decmin, decmax = list1[ci[i]].dec - tol1, list1[ci[i]].dec + tol1
        # Inch along in list2 until we find the part that matches list1 in RA
        while jmin < len(list2) and list2[cj[jmin]].ra < ramin:  jmin += 1
        if jmin == len(list2):  break
        # Now go through all the RA matches and check the RA+DEC distance.
        j = jmin
        while j < len(list2) and list2[cj[j]].ra < ramax:
            dr2 = ((CD*(list1[ci[i]].ra-list2[cj[j]].ra))**2 +
                       (list1[ci[i]].dec-list2[cj[j]].dec)**2)
            if dr2 < best_dr2[ci[i]]:
                match_idx[ci[i]], best_dr2[ci[i]] = cj[j], dr2
            j += 1
    # Save and return the index in list2 of the best-matched object to each
    # item in list1, and the distances between corresponding best matches.
    best_dr = [np.sqrt(dr2)*3600 for dr2 in best_dr2]
    return match_idx, best_dr


def match_coords_v2(list1, list2, tol = 1.5):
    """
    A new version of match_coords above, more complex but slightly faster.
    list1 and list2 are nominally lists of CatalogEntry instances, but could
    be lists of anything with 'ra' and 'dec' attributes in decimal degrees.
    """
    tol1, tol2 = tol/3600.0, (tol/3600.0)**2
    jmin, match_idx, best_dr2 = 0, [None]*len(list1), [tol2]*len(list1)
    # Unpack the RAs and DECs of everything into numpy arrays.
    list1ra = np.array([c.ra for c in list1])
    list2ra = np.array([c.ra for c in list2])
    list1dec = np.array([c.dec for c in list1])
    list2dec = np.array([c.dec for c in list2])
    # This next bit is rather tricky because we have to sort the lists,
    # but also remember where we put each item.  Here i, j = sorted indices,
    # while ci[i], ci[j] = corresponding actual indices in original lists.
    # Introduce transformed lists to economize on indexing.
    ci, cj = list1ra.argsort(), list2ra.argsort()
    slist1ra, slist1dec = list1ra[ci], list1dec[ci]
    slist2ra, slist2dec = list2ra[cj], list2dec[cj]
    # Since there's no sense in searching outside the bounds of the arrays,
    # extract the indices of only those elements within the rectangular
    # overlap area (we'll mostly be dealing with equatorial rectangles).
    # Reindex the lists to the originals.
    R0, R1 = max(min(list1ra), min(list2ra)), min(max(list1ra), max(list2ra))
    D0, D1 = max(min(list1dec), min(list2dec)), min(max(list1dec), max(list2dec))
    ci = ci[np.all([slist1ra >= R0, slist1dec >= D0,
                    slist1ra <= R1, slist1dec <= D1], axis=0)]
    cj = cj[np.all([slist2ra >= R0, slist2dec >= D0,
                    slist2ra <= R1, slist2dec <= D1], axis=0)]
    slist1ra, slist1dec = list1ra[ci], list1dec[ci]
    slist2ra, slist2dec = list2ra[cj], list2dec[cj]
    # Finally, start going through the lists again.
    for i in range(len(slist1ra)):
        decmin, decmax = slist1dec[i] - tol1, slist1dec[i] + tol1
        CD = np.cos(0.5*(decmin+decmax)*np.pi/180.0)
        ramin, ramax = slist1ra[i] - tol1/CD, slist1ra[i] + tol1/CD
        # Inch along in list2 until we find the part that matches list1 in RA
        while jmin < len(slist2ra) and slist2ra[jmin] < ramin:  jmin += 1
        # No point going past the end of the list
        if jmin == len(slist2ra):  break
        # Now go through all the RA matches and check the RA+DEC distance.
        j = jmin
        while j < len(slist2ra) and slist2ra[j] < ramax:
            # Check in the box before finding the angular distance
            if slist2dec[j] > decmin and slist2dec[j] < decmax:
                dr2 = ((CD*(slist1ra[i]-slist2ra[j]))**2 +
                           (slist1dec[i]-slist2dec[j])**2)
                if dr2 < best_dr2[ci[i]]:
                    match_idx[ci[i]], best_dr2[ci[i]] = cj[j], dr2
            j += 1
    # Save and return the index in list2 of the best-matched object to each
    # item in list1, and the distances between corresponding best matches.
    best_dr = [np.sqrt(dr2)*3600 for dr2 in best_dr2]
    return match_idx, best_dr


# ============================================================================
#                Classes for retrieval of online catalogs
# ============================================================================


class DS9Region(CatalogEntry):
    """
    On thinking about it, a DS9 region file is itself a kind of catalog.
    The parsing for DS9 region files is weird and variable, so I doubt
    any of the above generic classes will really do the trick; it has to
    be parsed in stages.
    """

    _mod_fields = ["coord","color","width","point","text","angle"]

    def _init_from_region_str(self, ra=None, dec=None, string=None):
        # Baseline default values; all we need to initialize a Region
        for mod in self._mod_fields:  setattr(self, mod, None)
        self.coord, self.shape, self.size = "j2000", "circle", 5.0/3600

        # If we're given a region string, try to parse it
        if string != None:
            try:
                # First separate the basic parameters out
                mymatch = re.match\
                    ("([a-zA-Z ]*)[(]?([\d\.\s,-]*)[)]?\s*#(.*)", string)
                shape, pars, mods = mymatch.groups()
                self.shape = shape.strip()
                pars = re.findall("([\w\.-]+)", pars)
                mods = dict(re.findall("(\w*)\s*=\s*(\w*)", mods))
                # Based on the shape, figure out what the pars mean
                super(DS9Region,self).__init__(pars[0], pars[1])
                point_subtype = re.match("(\w+)\s*point",self.shape)
                if self.shape == "circle":
                    self.size = float(pars[2])
                elif self.shape == "box":
                    self.size, self.angle = float(pars[3]), float(pars[4])
                elif self.shape == "text":
                    self.text = re.match("[{ ]?(.*)[} ]?",pars[2])
                    if self.text != None:  self.text = self.text.group(1)
                elif point_subtype != None:
                    self.shape == "point"
                    self.point == point_subtype.group(1)
                else:
                    raise BadInit("Unsupported ds9 region type:  {0}"
                                  .format(self.shape), string=string)
                # Finally, parse the mods
                for attr in self._mod_fields:
                    if mods.has_key(attr):    setattr(self, attr, mods[attr])
                    elif hasattr(self,attr):  pass
                    else:                     setattr(self, attr, None)
            except:
                raise BadInit("Trouble parsing ds9 region string:\n{0}"
                              .format(string), string=string)

        # Else, if RA and DEC are given, initialize from there
        elif ra != None and dec != None:
            super(DS9Region,self).__init__(ra,dec)

        # If no initialization given, fail epically
        else:
            raise BadInit ("No initialization information given to {0}"
                           .format(self.__class__.__name__))

    def region_str(self):
        # Now do all the above, except in reverse
        regstr = "{0}({1:.5f},{2:.5f}".format(self.shape,self.ra,self.dec)
        if self.shape == "circle":
            regstr += ",{0:.5f})".format(self.size)
        elif self.shape == "box":
            regstr += ",{0:.5f},{0:.5f},{1:.0f})".format(self.size,self.angle)
        elif self.shape == "point" or self.shape == "text":
            regstr += ")"
        else:
            raise BadInit ("I never should have gotten here in class {0}!"
                           .format(self.__class__.__name__))
        modlist = [m for m in self._mod_fields if getattr(self,m) != None]
        if len(modlist) > 0:
            regstr += " # " + ",".join(["{0}={1}".format(m,getattr(self,m))
                                       for m in modlist])
        return regstr

    def __init__(self, ra=None, dec=None, string=None):
        # Wraps region_line so it can be overloaded in derived classes
        self.__init__(ra=ra, dec=dec, string=string)

    def __str__(self):
        # Wraps region_line so it can be overloaded in derived classes
        return self.region_str()


class SkybotObject(AsciiRecord, KeywordRecord, CatalogEntry):
    """
    Represents a SkyBot object as read in from an ASCII file, in the format
    retrieved by wwwget.
    """

    _ascii_separator = " | "
    _ascii_fields = \
    [
        [ "id",      0, lambda x:      int(x), '{0:6d}',        0 ],
        [ "name",    1, lambda x:   x.strip(), '{0:20s}',      '' ],
        [ "ra_hour", 2, lambda x:    float(x), '{0:10.7f}',  None ],
        [ "dec",     3, lambda x:    float(x), '{0:9.5f}',   None ],
        [ "class",   4, lambda x:   x.strip(), '{0:8s}',     None ],
        [ "mV",      5, lambda x:    float(x), '{0:4.1f}',    0.0 ],
        [ "sigpos",  6, lambda x:    float(x), '{0:4.1f}',    0.0 ],
        [ "mu_tot",  7, lambda x:    float(x), '{0:5.1f}',    0.0 ],
        [ "mu_ra",   8, lambda x:    float(x), '{0:5.1f}',    0.0 ],
        [ "mu_dec",  9, lambda x:    float(x), '{0:5.1f}',    0.0 ],
        [ "Dg",     10, lambda x:    float(x), '{0:14.10f}', None ],
        [ "Dh",     11, lambda x:    float(x), '{0:14.10f}', None ],
    ]
    _keyword_fields = [ (f[0],f[4]) for f in _ascii_fields ]

    def __str__(self):
        ra = ephem.hours(self.ra*ephem.pi/180.0)
        dec = ephem.degrees(self.dec*ephem.pi/180.0)
        return "{0} ({1}, {2})".format(self.name, ra, dec)

    def __init__(self, *args, **kwargs):
        """
        Initialize from ASCII if possible, overriding with existing keywords.
        """
        if len(args)   > 0:
            self._init_from_ascii(*args)
            # The ASCII files store ra in hours, so convert to degrees.
            # RA as a keyword init argument should still be degrees.
            self.ra = 15.0*self.ra_hour
        if len(kwargs) > 0:  self._init_from_keywords(**kwargs)

    @staticmethod
    def webquery(ra=None, dec=None, datetime=None, cachefname=None,
                 rm=120, verbose=False):
        """
        Retrieves asteroids from Skybot in a FOV-shaped square around the
        central (RA,DEC) of a given SkyMapper exposure.  This will allow us
        to tag known rocks automatically in our data.
        RA and DEC must be in decimal degrees, datetime a datetime.datetime.
        """

        # If the user hasn't specified a full RA, DEC, and time, complain.
        if (ra is None or dec is None or datetime is None):
            print "SkybotObject.webquery():  I need to know where to look!"
            return [ ]
    
        # Call wwwget to query Skybot website and retrieve data as ASCII.
        if verbose:
            print "Querying Skybot database with RA = {0:.5f}, DEC = {1:.5f}"\
              .format(ra, dec)
        url = "http://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php"
        params = \
        [
            "-ep={0}".format(datetime),     # 2011-08-25T14:15:31
            "-ra={0:.5f}".format(ra),       # 300.00000 [degrees]
            "-dec={0:.5f}".format(dec),     # -35.00000 [degrees]
            "-rm={0:f}".format(rm),         # 120 [arcminutes]
            "-mime=text",
            "-output=basic",
            "-filter=1",
            "-loc=260",                     # Siding Spring
            "-from=skymapper"
        ]
        cmd = "wwwget '{0}?{1}'".format(url,"&".join(params))
        if verbose:  print cmd

        # Run the command and retrieve the results
        proc = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.STDOUT)
        (ascii_query, stderrstr) = proc.communicate()
        objlist = [ ]
        for line in ascii_query.split('\n'):
            try:  objlist.append(SkybotObject(line.strip()))
            except BadInit:  continue
        if cachefname != None:
            if verbose:
               print "Writing results of query to {0}".format(cachefname)
            SkybotObject.write_ascii_file(cachefname, objlist)
        return objlist


class VizieRObject(AsciiRecord, KeywordRecord, CatalogEntry):
    """
    Represents anything you can get through VizieR.  This is an abstract
    class as it does not implement _keyword_fields[] or _ascii_fields.
    It does implement an abstract "webquery" method, which queries VizieR
    from a virtual _vizcat to be implemented in derived classes.
    For formats of individual catalogs, including column names along with
    a basic query form, see
        http://vizier.u-strasbg.fr/viz-bin/VizieR
    which should be useful for deriving new wrapper classes.
    """

    _vizcat = None
    _vizpars = ""

    def __init__(self, *args, **kwargs):
        # This will initialize from ASCII if possible, and will override
        # the ASCII fields with keywords if they are provided.
        if len(args)   > 0:  self._init_from_ascii(*args)
        if len(kwargs) > 0:  self._init_from_keywords(**kwargs)

    @classmethod
    def webquery(cls, ra=None, dec=None, cachefname=None, rm=120, verbose=False):
        """
        Retrieves objects from VizieR in a cone around the central
        (RA,DEC) of a given SkyMapper exposure.  From this we can get crappy
        zeropoints until the interns get that stellar locus fit going.
        """

        # If the user hasn't specified a full RA and DEC, complain.
        # RA and DEC must be in decimal degrees.
        if (ra is None or dec is None):
            print "{0}.webquery():".format(cls.__name__)
            print "I need to know where to look!"
            return [ ]
    
        # Use vizquery to query VizieR and retrieve the data as ASCII.
        if verbose:
            print "Querying {0} via VizieR with RA = {1:.5f}, DEC = {2:.5f}"\
              .format(cls._vizcat, ra, dec)
        cmd = "vizquery -get -mime='|' -source='{0}' -out.max=100000 "\
              "-c={1:.5f}{2:+.5f},rm={3:.1f},eq=j2000 {4}".format(
                      cls._vizcat, ra, dec, rm, cls._vizpars)
        if verbose:  print cmd

        # Run the command and retrieve the results
        proc = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.STDOUT)
        (ascii_query, stderrstr) = proc.communicate()
        if cachefname != None:
            if verbose:
               print "Writing results of query to {0}".format(cachefname)
            with open(cachefname, "wt") as cache:
                cache.write(ascii_query)
        objlist = [ ]
        for line in ascii_query.split('\n'):
            # There's a lot of crap at the top of what vizquery returns.
            # For most catalogs, lines should have at least some numbers
            # and have no *lowercase* letters.
            # if re.search('[a-z]+', line) or not re.search('\d+', line):
            # RS 2013/09/20:  Requiring the first 3 characters of each line
            # to be digits.  Vet each new catalog to see if it's an issue.
            if not re.search('^\d{3}', line):
                continue
            try:  objlist.append(cls(line.strip()))
            except BadInit:  continue
        return objlist


class USNOB1_Object(VizieRObject):
    """
    Represents a USNO-B 1.0 catalog object as read in from an ASCII file,
    in the format retrieved by vizquery.
    """

    _vizcat = "USNO-B1.0"
    _vizpars = "-out=RAJ2000,DEJ2000,USNO-B1.0,e_RAJ2000,e_DEJ2000,"\
               "Epoch,pmRA,pmDE,Ndet,B1mag,R1mag,B2mag,R2mag,Imag,"\
               "B1s/g,R1s/g,B2s/g,R2s/g,Is/g"
    _ascii_separator = "|"
    _ascii_fields = \
    [
        [ "ra",      0, lambda x:        float(x), '{0:9.5f}', None ],
        [ "dec",     1, lambda x:        float(x), '{0:9.5f}', None ],
        [ "name",    2, lambda x:       x.strip(),  '{0:12s}',   '' ],
        [ "ra_err",  3, lambda x: 2.8e-7*float(x), '{0:5.1e}', None ],
        [ "dec_err", 4, lambda x: 2.8e-7*float(x), '{0:5.1e}', None ],
        [ "epoch",   5, lambda x:        float(x), '{0:6.1f}', None ],
        [ "pmRA",    6, lambda x:        float(x), '{0:6.1f}',  0.0 ],
        [ "pmDEC",   7, lambda x:        float(x), '{0:6.1f}',  0.0 ],
        [ "Ndet",    8, lambda x:          int(x),    '{0:d}',    0 ],
        [ "B1mag",   9, lambda x:        float(x), '{0:5.2f}',  0.0 ],
        [ "R1mag",  10, lambda x:        float(x), '{0:5.2f}',  0.0 ],
        [ "B2mag",  11, lambda x:        float(x), '{0:5.2f}',  0.0 ],
        [ "R2mag",  12, lambda x:        float(x), '{0:5.2f}',  0.0 ],
        [ "Imag",   13, lambda x:        float(x), '{0:5.2f}',  0.0 ],
        [ "B1sg",   14, lambda x:          int(x), '{0:2d}',      0 ],
        [ "R1sg",   15, lambda x:          int(x), '{0:2d}',      0 ],
        [ "B2sg",   16, lambda x:          int(x), '{0:2d}',      0 ],
        [ "R2sg",   17, lambda x:          int(x), '{0:2d}',      0 ],
        [ "Isg",    18, lambda x:          int(x), '{0:2d}',      0 ],
    ]
    _keyword_fields = [ (f[0],f[4]) for f in _ascii_fields ]

class SDSSObject(VizieRObject):
    """
    Represents a SDSS Photometry Catalog DR9 catalog object as read in from an ASCII file,
    in the format retrieved by vizquery.
    """

    _vizcat = "SDSS9"
    _vizpars = "-out=RAJ2000,DEJ2000,SDSS9,e_RAJ2000,e_DEJ2000,"\
               "rmag,cl,spCl,zsp,e_zsp,f_zsp "\
               "-out.max=1000000 mode=1|2 rmag=<=21.0 "
    _ascii_separator = "|"
    _ascii_fields = \
    [
        [ "ra",      0, lambda x:        float(x), '{0:9.5f}', None ],
        [ "dec",     1, lambda x:        float(x), '{0:9.5f}', None ],
        [ "name",    2, lambda x:       x.strip(),  '{0:20s}',   '' ],
        [ "ra_err",  3, lambda x: 2.8e-7*float(x), '{0:5.1e}', None ],
        [ "dec_err", 4, lambda x: 2.8e-7*float(x), '{0:5.1e}', None ],
        [ "rmag",    5, lambda x:        float(x), '{0:5.2f}',  0.0 ],
        [ "phtype",  6, lambda x:          int(x),   '{0:2d}',    0 ],
        [ "sptype",  7, lambda x:       x.strip(),  '{0:10s}',   '' ],
        [ "zsp",     8, lambda x:        float(x), '{0:8.5f}',  -1. ],
        [ "zsp_err", 9, lambda x:        float(x), '{0:7.5f}',  0.0 ],
        [ "zsp_flag",10, lambda x:          int(x),   '{0:3d}',    0 ],
   ]
    _keyword_fields = [ (f[0],f[4]) for f in _ascii_fields ]

class TwoMassPSCObject(VizieRObject):
    """
    Represents a 2MASS Point Source Catalog object as read in from
    an ASCII file, in the format retrieved by vizquery.
    """

    _vizcat = "2MASS-PSC"
    _ascii_separator = "|"
    _ascii_fields = \
    [
        [ "ra",       0, lambda x: float(x), '{0:9.5f}', None ],
        [ "dec",      1, lambda x: float(x), '{0:9.5f}', None ],
        [ "name",     2, lambda x: x.strip(), '{0:12s}',   '' ],
        [ "Jmag",     3, lambda x: float(x), '{0:5.2f}', None ],
        [ "Jmag_err", 4, lambda x: float(x), '{0:5.2f}', None ],
        [ "Hmag",     5, lambda x: float(x), '{0:5.2f}', None ],
        [ "Hmag_err", 6, lambda x: float(x), '{0:5.2f}', None ],
        [ "Kmag",     7, lambda x: float(x), '{0:5.2f}', None ],
        [ "Kmag_err", 8, lambda x: float(x), '{0:6.2f}',    0 ],
        [ "Qflg",     9, lambda x: x.strip(),  '{0:3s}',  0.0 ],
        [ "Rflg",    10, lambda x: x.strip(),  '{0:3s}',  0.0 ],
        [ "Bflg",    11, lambda x: x.strip(),  '{0:3s}',  0.0 ],
        [ "Cflg",    12, lambda x: x.strip(),  '{0:3s}',  0.0 ],
        [ "Xflg",    13, lambda x: x.strip(),  '{0:3s}',  0.0 ],
        [ "Aflg",    14, lambda x: x.strip(),  '{0:3s}',  0.0 ],
    ]
    _keyword_fields = [ (f[0],f[4]) for f in _ascii_fields ]


class TwoMassXSCObject(VizieRObject):
    """
    Represents a 2MASS Extended Source Catalog object as read in from
    an ASCII file, in the format retrieved by vizquery.
    """

    _vizcat = "2MASX"
    _ascii_separator = "|"
    _ascii_fields = \
    [
        [ "name",      0, str,   '{0:12s}',    '' ],
        [ "ra",        2, float, '{0:9.5f}', None ],
        [ "dec",       3, float, '{0:9.5f}', None ],
        [ "Kba",       4, float, '{0:9.5f}', None ],
        [ "rfe",       5, float, '{0:9.5f}', None ],
        [ "Jmag",      6, float, '{0:5.2f}', None ],
        [ "Jmag_err",  7, float, '{0:5.2f}', None ],
        [ "Hmag",      8, float, '{0:5.2f}', None ],
        [ "Hmag_err",  9, float, '{0:5.2f}', None ],
        [ "Kmag",     10, float, '{0:5.2f}', None ],
        [ "Kmag_err", 11, float, '{0:6.2f}',    0 ],
        [ "xid",      12, int,     '{0:3s}',  0.0 ],
    ]
    _keyword_fields = [ (f[0],f[4]) for f in _ascii_fields ]


class PesstoObject(AsciiRecord, KeywordRecord, CatalogEntry):
    """
    Represents a PESSTO object as read in from an ASCII file, in the format
    retrieved by wwwget.
    """

    _ascii_separator = "|"
    _ascii_fields = \
    [
        [ "name",     0, str.strip, '{0:20s}',   '' ],
        [ "ra_str",   1, str.strip, '{0:11s}', None ],
        [ "dec_str",  2, str.strip, '{0:11s}', None ],
        [ "discdate", 3, str.strip, '{0:8s}',  None ],
        [ "source",   4, str.strip, '{0:8s}',    '' ],
        [ "mag",      5, float,     '{0:4.1f}', 0.0 ],
        [ "z",        6, float,     '{0:5.3f}', 0.0 ],
        [ "type",     7, str.strip, '{0:s}',     '' ],
    ]
    _keyword_fields = [ (f[0],f[4]) for f in _ascii_fields ]
    _keyword_fields += [ ('atel', None) ]

    def __str__(self):
        return "{0} ({1}, {2})".format(self.name, self.ra_str, self.dec_str)

    def __init__(self, *args, **kwargs):
        """
        Initialize from ASCII if possible, overriding with existing keywords.
        """
        if len(args) > 0:
            self._init_from_ascii(*args)
            # The ASCII files store ra and dec in sexigesimal; convert.
            # RA and DEC as keyword init arguments should still be degrees.
            self.ra = 180.0/ephem.pi * ephem.hours(self.ra_str)
            self.dec = 180.0/ephem.pi * ephem.degrees(self.dec_str)
        if len(kwargs) > 0:  self._init_from_keywords(**kwargs)

    @staticmethod
    def webquery(atel, cachefname=None, verbose=False):
        """
        Retrieves an ATel from the Web and parses it, returning a list of all
        PESSTO objects in the ATel.
        """

        # If the user hasn't given us a sensible ATel number, complain.
        if not isinstance(atel, int) or atel < 0:
            print "PesstoObject.webquery(): bad ATel number"
            return [ ]
        # Call wwwget to query ATel website and retrieve data as ASCII.
        if verbose:
            print "Querying web for ATel", atel
        url = "http://www.astronomerstelegram.org/?read={0}".format(atel)
        cmd = "wwwget '{0}'".format(url)
        if verbose:  print cmd

        # Run the command and retrieve the results
        proc = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.STDOUT)
        (ascii_query, stderrstr) = proc.communicate()
        objlist = [ ]
        for line in ascii_query.split('\n'):
            try:  objlist.append(PesstoObject(line.strip(), atel=atel))
            except BadInit:  continue
        if cachefname != None:
            if verbose:
               print "Writing results of query to {0}".format(cachefname)
            PesstoObject.write_ascii_file(cachefname, objlist)
        return objlist


class CrtsObject(PesstoObject):
    """
    Represents a CRTS object as read in from a CSV file.
    """

    _ascii_separator = None
    _ascii_fields = \
    [
        [ "ra",   3, float, '{0:10.6f}', 0.0 ],
        [ "dec",  4, float, '{0:10.6f}', 0.0 ],
        [ "mag",  5, float, '{0:6.3f}',  0.0 ],
        [ "jd",   6, float, '{0:14.6f}', 0.0 ],
    ]
    _keyword_fields = [ (f[0],f[4]) for f in _ascii_fields ]

    def __init__(self, *args, **kwargs):
        """
        Initialize from ASCII if possible, overriding with existing keywords.
        """
        if len(args) > 0:  self._init_from_ascii(*args)
        if len(kwargs) > 0:  self._init_from_keywords(**kwargs)
        discdate = ephem.Date(self.jd - 2415020.0)
        self.ra_str = str(ephem.hours(self.ra * ephem.pi/180))
        if self.ra < 150.0:
            self.ra_str = '0' + self.ra_str
        self.dec_str = str(ephem.degrees(self.dec * ephem.pi/180))
        if abs(self.dec) < 10.0:
            self.dec_str = self.dec_str[0] + '0' + self.dec_str[1:]
        if self.dec > 0:
            self.dec_str = '+' + self.dec_str
        self.date_str = discdate.datetime().strftime("%d %b %Y")
        self.name = "CSS{}:{}{}".format(
            discdate.datetime().strftime("%y%m%d"),
            self.ra_str.replace(':','')[:6],
            self.dec_str.replace(':','')[:7])

    @staticmethod
    def webquery(cachefname=None, verbose=False):
        """
        Retrieves the CRTS list from the Web and parses it, returning a list
        of all CRTS (probable) supernovae.
        """

        # Call wwwget to query CRTS website and retrieve data as ASCII.
        if verbose:
            print "Querying CRTS Supernova webpage"
        url = "http://nesssi.cacr.caltech.edu/catalina/Results.outall2xSNfind"
        cmd = "wwwget '{0}'".format(url)
        if verbose:  print cmd

        # Run the command and retrieve the results
        proc = sp.Popen(shlex.split(cmd), stdout=sp.PIPE, stderr=sp.STDOUT)
        (ascii_query, stderrstr) = proc.communicate()
        objlist = [ ]
        for line in ascii_query.split('\n'):
            try:  objlist.append(CrtsObject(line.strip()))
            except BadInit:  continue
        if cachefname != None:
            if verbose:
               print "Writing results of query to {0}".format(cachefname)
            PesstoObject.write_ascii_file(cachefname, objlist)
        return objlist


class LSQObject(PesstoObject):
    """
    Represents a LaSilla-QUEST object as read in from a tab-delimited file.
    """

    _ascii_separator = '\t'
    _ascii_fields = \
    [
        [ "name",     0, str.strip, '{:<8}',   ' ' ],
        [ "discdate", 1, str.strip, '{}',      ' ' ],
        [ "lastobs",  2, str.strip, '{}',      ' ' ],
        [ "type",     3, str.strip, '{}',      ' ' ],
        [ "cattype",  4, str.strip, '{}',      ' ' ],
        [ "spectype", 5, str.strip, '{}',      ' ' ],
        [ "discmag",  6, float,     '{:6.3f}', 0.0 ],
        [ "lastmag",  7, float,     '{:6.3f}', 0.0 ],
        [ "host_z",   8, float,     '{:6.3f}', 0.0 ],
        [ "sn_z",     9, float,     '{:6.3f}', 0.0 ],
        [ "ra_str",  12, str.strip, '{:11}',   0.0 ],
        [ "dec_str", 13, str.strip, '{:11}',   0.0 ],
    ]
    _keyword_fields = [ (f[0],f[4]) for f in _ascii_fields ]

    def __init__(self, *args, **kwargs):
        """
        Initialize from ASCII if possible, overriding with existing keywords.
        """
        if len(args) > 0:  self._init_from_ascii(*args)
        if len(kwargs) > 0:  self._init_from_keywords(**kwargs)
        self.ra = ephem.hours(self.ra_str)*180/ephem.pi
        self.dec = ephem.degrees(self.dec_str)*180/ephem.pi
        if self.sn_z is None:  self.sn_z = 0.0
        if self.host_z is None:  self.host_z = 0.0


# ============================================================================
#        Classes for pipeline-related catalogs internal to SkyMapper
# ============================================================================


class SExtractorDetection(FITSRecord, CatalogEntry):
    """
    Represents a SExtractor detection as read in from a binary FITS file.
    """

    _fits_fields = \
    [
        [ "id",      "NUMBER",       "1J", "",      "I4",    0    ],
        [ "x",       "X_IMAGE",      "1E", "pixel", "F7.2",  0.0  ],
        [ "y",       "Y_IMAGE",      "1E", "pixel", "F7.2",  0.0  ],
        [ "a",       "A_IMAGE",      "1E", "pixel", "F5.2",  0.0  ],
        [ "b",       "B_IMAGE",      "1E", "pixel", "F5.2",  0.0  ],
        [ "e",       "ELONGATION",   "1E", "",      "F5.3",  0.0  ],
        [ "theta",   "THETA_IMAGE",  "1E", "deg",   "F5.1",  0.0  ],
        [ "ra",      "ALPHA_J2000",  "1D", "deg",   "F10.6", None ],
        [ "dec",     "DELTA_J2000",  "1D", "deg",   "F10.6", None ],
        [ "fwhm",    "FWHM_IMAGE",   "1E", "pixel", "F5.2",  0.0  ],
        [ "flux",    "FLUX_APER",    "7E", "count", "G10.5", 0.0  ],
        [ "fluxerr", "FLUXERR_APER", "7E", "count", "G10.5", 0.0  ],
        [ "flag",    "FLAGS",        "1I", "",      "I3",    0    ],
        [ "star",    "CLASS_STAR",   "1E", "",      "F5.3",  0.0  ],
        [ "quikmag", "MAG_BEST",     "1E", "mag",   "F5.3",  0.0  ],
    ]


class SExtractorWindowedDetection(FITSRecord, CatalogEntry):
    """
    Same as SExtractorDetection, but uses windowed measurement values.
    """

    _fits_fields = \
    [
        [ "id",      "NUMBER",         "1J", "",      "I4",    0    ],
        [ "x",       "XWIN_IMAGE",     "1E", "pixel", "F7.2",  0.0  ],
        [ "y",       "YWIN_IMAGE",     "1E", "pixel", "F7.2",  0.0  ],
        [ "a",       "AWIN_IMAGE",     "1E", "pixel", "F5.2",  0.0  ],
        [ "b",       "BWIN_IMAGE",     "1E", "pixel", "F5.2",  0.0  ],
        [ "e",       "ELONGATION",     "1E", "",      "F5.3",  0.0  ],
        [ "theta",   "THETAWIN_IMAGE", "1E", "deg",   "F5.1",  0.0  ],
        [ "ra",      "ALPHA_J2000",    "1D", "deg",   "F10.6", None ],
        [ "dec",     "DELTA_J2000",    "1D", "deg",   "F10.6", None ],
        [ "fwhm",    "FWHMWIN_IMAGE",  "1E", "pixel", "F5.2",  0.0  ],
        [ "flux",    "FLUX_APER",      "7E", "count", "G10.5", 0.0  ],
        [ "fluxerr", "FLUXERR_APER",   "7E", "count", "G10.5", 0.0  ],
        [ "flag",    "FLAGS",          "1I", "",      "I3",    0    ],
        [ "star",    "CLASS_STAR",     "1E", "",      "F5.3",  0.0  ],
        [ "quikmag", "MAG_BEST",       "1E", "mag",   "F5.3",  0.0  ],
    ]


class SkymapperCalibStar(FITSRecord, KeywordRecord, CatalogEntry):
    """
    Represents a list of sources we might use to get photometric zeropoints
    for our fields.  Will belong in the REF cache ultimately.
    """

    # Generate filter attributes + keywords below, e.g., "u_err", "UMAGERR".
    # _nodata is an effectively infinite error value to gracefully de-weight
    # missing data in calibration algorithms. _zpcolors shows color combos
    # for calculating zeropoints, e.g., g' = zp_g + C_g*(g - r).
    _nodata = 99.99
    _filters = ['u', 'v', 'g', 'r', 'i', 'z']
    _zpcolors = [('u', 'g'), ('v', 'g'), ('g', 'r'),
                 ('g', 'r'), ('r', 'i'), ('i', 'z')]
    _fits_fields = \
    [
        [ "ra",  "ALPHA_J2000", "1D", "deg", "F10.6", None ],
        [ "dec", "DELTA_J2000", "1D", "deg", "F10.6", None ],
    ]
    for f in _filters:
        magattr, magkw = f, f.upper() + "MAG"
        errattr, errkw = f + "_err", f.upper() + "MAGERR"
        _fits_fields += [ [ magattr, magkw, "1E", "mag", "F5.3", 15.00 ],
                          [ errattr, errkw, "1E", "mag", "F5.3", _nodata ], ]

    _keyword_fields = [ (f[0], f[5]) for f in _fits_fields ]

    def __init__(self, *args, **kwargs):
        """
        Try to init in this order:  from a FITS table row (2 arguments);
        then override with keywords as available.
        """
        if len(args) == 2:   self._init_from_fitsrow(*args)
        self._init_from_keywords(**kwargs)
        self._calc_zpcolors()

    def _calc_zpcolors(self):
        """
        Once all the magnitudes are in there, form the colors which will
        be needed to achieve a color solution when calibrating.
        """
        for f, zpc in zip(self._filters, self._zpcolors):
            m0, merr0 = getattr(self, zpc[0]), getattr(self, zpc[0] + '_err')
            m1, merr1 = getattr(self, zpc[1]), getattr(self, zpc[1] + '_err')
            if merr0 < 9.99 and merr1 < 9.99:
                setattr(self, f + 'zpc', m0 - m1)
                setattr(self, f + 'zpc_err', np.sqrt(merr0**2 + merr1**2))
            else:
                setattr(self, f + 'zpc', self._nodata)
                setattr(self, f + 'zpc_err', self._nodata)

    def merge(self, fill, improve=False):
        """
        Merge in magnitudes from another SkymapperCalibStar instance to fill
        empty fields in this instance.  If the "improve" keyword is set,
        replace more accurate magnitudes in self with the fill values.
        """
        for mag in self._filters:
            magerr, col, colerr = [mag + s for s in '_err', 'zpc', 'zpc_err']
            smag, smagerr = getattr(self, mag), getattr(self, magerr)
            fmag, fmagerr = getattr(fill, mag), getattr(fill, magerr)
            if smagerr > 9.99 or (improve and smagerr > fmagerr):
                setattr(self, mag, fmag)
                setattr(self, magerr, fmagerr)
        # Recalculate colors after we're done
        self._calc_zpcolors()


class APASSObject(PGRecord, SkymapperCalibStar):
    """
    Represents a record from the APASS DR6 hosted via the Postgres catalog
    Simon has recently set up on martini (2012/06/20).  Uses PGrecord.
    """

    _nodata = SkymapperCalibStar._nodata
    _dbname = 'catalogs'
    _dbhost = 'martini'
    _dbuser = 'untrusted'
    _dbpw = 'trustno1'
    _sqlquery = "SELECT * FROM apass_dr6 WHERE " \
                "q3c_radial_query(ra, decl, %s, %s, %s)"

    _sql_fields = \
    [
        [ "ra", "ra", None ], [ "dec",   "decl", None ],
        [ "g",  "g", 15.00 ], [ "g_err", "e_g", _nodata ],
        [ "r",  "r", 15.00 ], [ "r_err", "e_r", _nodata ],
        [ "i",  "i", 15.00 ], [ "i_err", "e_i", _nodata ],
    ]

    def __init__(self, *args, **kwargs):
        """
        Try to init in this order:  from a FITS table row (2 arguments);
        then override with SQL and keywords.
        """
        if len(args) == 2:
            self._init_from_fitsrow(*args)
        if len(kwargs) > 0:
            self._init_from_sqlrow(**kwargs)
            self._init_from_keywords(**kwargs)
        self._calc_zpcolors()
        # RS 2014/06/30:  Don't try this at home, kids, but to get something
        # like a reasonable zeropoint for v, I'm just going to zeropoint the
        # v filter using the g magnitude, and use a g-r color term to make up
        # the difference.  This isn't really going to be right, and I would
        # never advocate this e.g. for metallicity science.  But for SN work
        # it might not be completely horrible.  Let's see what happens.
        self.v, self.v_err = self.g, self.g_err
        self.vzpc, self.vzpc_err = self.gzpc, self.gzpc_err
