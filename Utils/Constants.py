#!/usr/bin/env python

# ============================================================================
# RS 2012/02/24:  Global Constants module
# ----------------------------------------------------------------------------
# This code is a central namespace for any variable that gets referred to in
# more than one part of the pipeline and can't really be made strictly local.
# Good examples include SkyMapper detector characteristics, structures of
# directory trees where files are stored, and so forth.
# ============================================================================

import re
import os
import time
import datetime
import subprocess
import ephem
import pyfits
import numpy as np
import copy
from Utils.Record import AsciiRecord
from Utils.Catalog import CatalogEntry
from Utils.TrackableException import TrackableException

def svn_version():
    cmd = "svnversion -n " + os.environ['SUBPIPEHOME']
    proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    return proc.communicate()[0]

# Is this a production install?
is_production_install = bool(int(os.environ["SUBPIPE_PRODUCTION_INSTALL"]))

# ----------------------------------------------------------------------------
#          Imager constants like gain, well depth, pixel scale etc.
# ----------------------------------------------------------------------------


class Imager(object):
    """Namespace wrapper for global detector information"""
    gain = 1.3                      # CCD gain (e-/ADU)
    saturate_adu = 45000            # Saturation level (ADU)
    read_noise = 10                 # Read noise (e-)
    rpsf = 15                       # PSF cutoff radius (pix) used in hotpants
    pixel_scale = 0.497             # Pixel scale in arcsec/pix
    flipx = 1                       # Parity flip for x w/rt (RA,DEC)
    flipy = -1                      # Parity flip for y w/rt (RA,DEC)

    # RS 2012/07/12:  Rotational mapping of CCDs to sky-invariant "subfields".
    # The subfield equals the CCD number for ROTSKYPA = 0, and is given as
    # follows for ROTSKYPA = 180:  subfield = ccd_map[ccd], ccd in [1..32].
    # This will help us with subtraction bookkeeping.  See Fang's email of
    # 21 July 2012 for details.
    ccd_map = [0, 29, 30, 31, 32, 25, 26, 27, 28,
                  21, 22, 23, 24, 17, 18, 19, 20,
                  13, 14, 15, 16,  9, 10, 11, 12,
                   5,  6,  7,  8,  1,  2,  3,  4]

    # Zeropoints in SkyMapper filters assuming typical SSO sky brightness
    # at airmass = 1, seeing = 1.5", and sky-limited 1-second exposure.
    # For guesstimating calibration w/ no catalog matches in a given field.
    # From SkyMapper pages at http://rsaa.anu.edu.au/observatories/
    filters = ['u', 'v', 'g', 'r', 'i', 'z']
    fake_zp = { 'u': 21.5, 'v': 21.3, 'g': 21.9,
                'r': 21.6, 'i': 21.0, 'z': 20.6, }

    # Instrumental zeropoints predicted from LPNHE model and validated on
    # real SkyMapper data; these are for AB magnitudes with fluxes in e-/s,
    # projected to airmass 1.
    zpinst = { 'u': 21.6, 'v': 21.8, 'g': 24.1,
               'r': 23.9, 'i': 23.1, 'z': 22.4, }

    # Night sky background levels at Siding Spring as measured by LPNHE,
    # in magnitudes per square arcsec
    nsb_mag = { 'u': 25.0, 'v': 22.5, 'g': 21.5,
                'r': 20.7, 'i': 19.3, 'z': 25.0 }

    # Night sky background levels at Siding Spring as measured by LPNHE,
    # in counts per pixel per second
    nsb_cps = dict([(k, (0.25*10**(-0.4*(nsb_mag[k] - zpinst[k]))))
                    for (k, v) in nsb_mag.items()])

    # RS 2013/06/11:  Characteristics of SkyMapper more relevant to
    # scheduling calculations.
    overhead = 20                   # Observation overhead (seconds)
    min_airmass = 1.0               # Minimum airmass at which to observe
    max_airmass = 1.4               # Maximum airmass at which to observe
    priority_thresh = 50.0          # Priority threshold at which to observe

    min_alt = np.arcsin(1.0/max_airmass)*180/np.pi
    max_alt = np.arcsin(1.0/min_airmass)*180/np.pi


def sso_sitegen():
    """Generates an ephem.Observer() for Siding Spring Observatory
    
    All this really does is to set some of the ephem.Observer() properties
    to give sensible answers, because copy.copy() doesn't respect them.
    Not even copy.deepcopy() works!
    """

    site = ephem.Observer()             # Siding Spring Observatory:
    site.lon = '149.066086'             #    longitude (decimal degrees)
    site.lat = '-31.277039'             #    longitude (decimal degrees)
    site.horizon = str(int(Imager.min_alt))
    return site

try:
    sso = sso_sitegen()
except:
    pass

follow_types = ['Cand', 'CandFaint','SN', 'Ia', 'Ibc', 'II', 'IIn', 'IIP', 'IIL', 'SLSN']

FLAG_BADREF=80
FLAG_OFFDISK=99
MAX_REFELONG=1.3
MAX_REFSEEING=10 #pixel
MIN_REFDEPTH={'g':20.,'v':15.,'r':20.,'i':20.,'u':15.,'z':15.}  #check for different filters

# ----------------------------------------------------------------------------
#                     Pathname and filename conventions
# ----------------------------------------------------------------------------


class PipelinePath(object):
    """Namespace wrapper for pipeline directory structure"""
    home = os.environ["SUBPIPEHOME"]
#    ext = os.environ["SUBPIPEEXT"]
    bin = home + "/bin"                     # Binaries of useful executables
    etc = home + "/etc"                     # Global configuration files
    data = home + "/data"                   # Pipeline staging + output
    raw = data + "/raw"                     # Raw SkyMapper mosaics staging
    cal = data + "/cal"                     # Calibration data staging
    masks = cal + "/masks"                  # Bad pixel masks
    flats = cal + "/flats"                  # Bias-subtracted flats by CCD
    newflats = cal + "/newflats"            # Superflats by CCD ready for use
    new = data + "/new"                     # NEW images to be subtracted
    new_glob_str = new + "/*/*/*/*.fits"
    ref = data + "/ref"                     # REF images for subtraction
    ref_glob_str = ref + "/*/*/*/*.fits"
    sub = data + "/sub"                     # Finished SUBs
    sub_glob_str = sub + "/*/*/*/*.fits"
    cand = data + "/cand_info"              # Info on detected candidates
    xref = cand + "/xref"                   # Cross-referenced detections
    lcdata = cand + "/lcdata"               # Light curve data
    skybot = cand + "/skybot"               # Skybot cache data
    sdss = cand + "/sdss"                   # SDSS cache data
    django = home + "/mydjango"             # Source path for django website
    django_static = django + "/staticfiles" # Source path for django statics
    django_data = data + "/mydjango"        # django website related data
    django_db = django_data + "/db/test.db" # SQLite3 database for django
    django_media = django_data + "/media"   # Files uploaded by django users
    thumbpath = django_media + "/thumbs"    # Thumbnails
    finderpath = django_media + "/fcharts"  # Finding charts
    scratch = os.environ["SUBSCRATCH"]      # Scratch space (/ramdisk)
    scratchetc = scratch + "/etc"
    # Configuration files for various external programs
    SEx_files = [ "default.conv","STAP_SEx.sex", "STAP_SEx.nnw", "STAP_SEx.par",
                  "gauss_1.5_3x3.conv", "gauss_2.0_3x3.conv",
                  "gauss_2.0_5x5.conv", "gauss_2.5_5x5.conv",
                  "gauss_3.0_5x5.conv", "gauss_3.0_7x7.conv",
                  "gauss_4.0_7x7.conv", "gauss_5.0_9x9.conv" ]
    SWarp_config = [ "STAP_SWarp.swarp", "STAP_flats.swarp" ]
    random_pickle = [ "randomforest.pkl" ]
    etc_files = SEx_files + SWarp_config + random_pickle
    scratchetc_files = [ "{0}/{1}".format(scratchetc,fn) for fn in etc_files ]
    # Accounting of what's been recently observed (for scheduler)
    lastobs_state_fn = etc + "/skymapper_lastobs.txt"
    # List of known SNe for VetBot to auto-tag
    known_sne = etc + "/snradec.txt"
    # Pause lock filename for temporarily halting job submission during a run
    pauselockfname = etc + "/pause.switch"
    # SExtractor settings related to Fang's flux subtraction?
    fsatconf = "fsat.sex"
    fsatpar = "fsat.param"

# Miscellaneous filename conventions

def sextractor_output(fname):
    """SExtractor detection FITS table name from FITS image name"""
    return fname.replace(".fits",".fits.stars")

def ds9_region_file(fname):
    """ds9 region file name from FITS image name"""
    return fname.replace(".fits",".fits.reg")

def get_runtag(start_time):
    """Return unique runtag corresponding to the current time"""
    return time.strftime("%Y%m%d_%H%M%S", time.localtime(start_time))

def datetime2str(mydt):
    """Our favorite string format for datetime.datetime instances"""
    return mydt.strftime("%Y-%m-%dT%H:%M:%S")

def str2datetime(mystr):
    """Our favorite string format for datetime.datetime instances"""
    return datetime.datetime.strptime(mystr, "%Y-%m-%dT%H:%M:%S")

def candidate_name(ra, dec):
    """Standard naming convention for SkyMapper candidates"""
    # Use (RA,DEC) converted to sexagesimal to the nearest arcsec.
    coo = ephem.Equatorial(ra*ephem.pi/180.0, dec*ephem.pi/180.0)
    ra_str  = re.sub('[:.]','',str(coo.ra))[-1]
    dec_str = re.sub('[:.]','',str(coo.dec))[-1]
    if dec > 0: dec_str = "+" + dec_str
    return "SKY_J{0}{1}".format(ra_str,dec_str)

def get_obsseq(filename):
    match = re.search("(\d+-\d+-\d+T\d+:\d+:\d+)", filename)
    if match != None:  return match.group(1)
    else: return None

def get_obsid(filename):
    match = re.search("(\d+)_\d+-\d+-\d+T\d+:\d+:\d+", filename)
    if match != None:  return match.group(1)
    else: return ''
    
def get_basename(filename):
    match = re.search("(\S+\d+-\d+-\d+T\d+:\d+:\d+)",filename)
    if match != None:  return match.group(1)
    else: return None
   
def fits_exists(filename):
    if os.path.exists(filename):
        return filename
    elif os.path.exists(filename+'.gz'):
        return filename+'.gz'
    else:
        return None

class Filenames(object):
    """A class to do storage bookkeeping associated with SkyMapper images.

    It's basically a namespace to store filename conventions, but it also
    reads related metadata from the headers.  This class is a base class,
    but it can be used on anything with a header to provide some baseline
    global information on any image accessible from the whole pipeline.
    """

    def __init__(self, filename, runtag=""):
        """Initializes SkyMapper filenames according to conventions."""
        with pyfits.open(filename) as hdulist:
            hdr = hdulist[0].header
            # Keywords which will always be present
            self.filter = hdr['FILTNAME']
            self.imgtype = hdr['IMAGETYP']
            self.airmass = hdr['AIRMASS']
            self.exptime = hdr['EXPTIME']
            self.dateobs = str2datetime(hdr['DATE-OBS'])
            self.obstype = hdr['IMAGETYP']
            self.basefname = os.path.basename(filename)
            # Keywords which may not be present:
            # SkyMapper field ID
            if 'FIELD_ID' in hdr:  self.field = int(hdr['FIELD_ID'])
            else:  self.field = 0
            # Physical CCD number
            if 'IMAGEID' in hdr:  self.ccd = int(hdr['IMAGEID'])
            else:  self.ccd = 0
            # Virtual CCD number, or subfield number
            if 'SUBFIELD' in hdr:  self.subfield = int(hdr['SUBFIELD'])
            else:  self.subfield = self.ccd
            # Seeing median and normalized interquartile width
            if 'SEEING' in hdr:  self.seeing = float(hdr['SEEING'])
            else:  self.seeing = 99.9
            if 'SEEWID' in hdr:  self.seewid = float(hdr['SEEWID'])
            else:  self.seewid = 99.9
            # Elongation median and normalized interquartile width
            if 'ELONG' in hdr:  self.elong = float(hdr['ELONG'])
            else:  self.elong = 99.9
            if 'ELONGWID' in hdr:  self.elongwid = float(hdr['ELONGWID'])
            else:  self.elongwid = 99.9
            # Zeropoint
            if 'ZPMAG' in hdr:  self.zpmag = float(hdr['ZPMAG'])
            else:  self.zpmag = 0.0
            # Standard deviation of background pixels
            if 'BKGSIG' in hdr:  self.bkgsig = float(hdr['BKGSIG'])
            else:  self.bkgsig = 0.0
            if 'BKGMED' in hdr:  self.bkgmed = float(hdr['BKGMED'])
            else:  self.bkgmed = 0.0
            if 'MAGLIM50' in hdr:  self.maglim50 = float(hdr['MAGLIM50'])
            else:  self.maglim50 = 0.0
            if 'MAGLIM95' in hdr:  self.maglim95 = float(hdr['MAGLIM95'])
            else:  self.maglim95 = 0.0
            if 'MINRA' in hdr:
                minra, mindec = hdr['MINRA'], hdr['MINDEC']
                maxra, maxdec = hdr['MAXRA'], hdr['MAXDEC']
                coomin = ephem.Equatorial(minra, mindec)
                coomax = ephem.Equatorial(maxra, maxdec)
                self.minra,self.mindec=coomin.ra,coomin.dec
                self.maxra,self.maxdec=coomax.ra,coomax.dec
            else:
                self.minra,self.mindec=0.,0.
                self.maxra,self.maxdec=0.,0.
            if 'PROGRAM' in hdr: self.program = hdr['PROGRAM']
            else: 
                try: 
                    self.program = decode_obs_id(filename.split('/')[-1].split('_')[1])[0]
                except:
                    self.program = 'unknown'
        # Fill dummy "id" field so we can use this for job descriptions
        self.id = None
        self.runtag = runtag
        # "obsseq" (the date and time at which the OB was sent to TAROS),
        # isn't yet available as a keyword, so grep the filename.
        self.obsseq = get_obsseq (self.basefname)
        self.obsid  = get_obsid (self.basefname)
        # RS 2012/08/30:  Add Skybot cache name for this pointing.
        self.skybot_cachefname = "{0}/{1}/skybot_{2}.txt".format(
                PipelinePath.skybot, self.field, datetime2str(self.dateobs))
        self.sdss_cachefname = "{0}/{1}/sdss_{1}.txt".format(
                PipelinePath.sdss, self.field)
        # FY 2015/08/03: Add skymapper observation type
        


class FilenamesFlat(Filenames):
    """Filenames subclass for SkyMapper superflats made from twilight flats"""

    def __init__(self, filename, runtag=""):
        super(FilenamesFlat, self).__init__(filename, runtag=runtag)
        self.basedir = PipelinePath.flats
        self.relative_dir = "{0}/{1:02d}".format(self.filter, self.ccd)
        self.absolute_dir = "{0}/{1}".format(self.basedir, self.relative_dir)
        #FY - this should be for testing only
        #if not os.path.exists("{0}/{1}".format(self.absolute_dir,filename)):
        #    self.basedir = PipelinePath.new
        #    self.absolute_dir = "{0}/{1}".format(self.basedir, self.relative_dir)
        self.flatfname = "flat_{0}_{1}_{2}_{3}.fits".format(
                self.obsid, self.obsseq, self.filter, self.ccd)
        self.wtmapfname = "weightmap_{0}_{1}_{2}_{3}.fits".format(
                self.obsid, self.obsseq, self.filter, self.ccd)
        self.log_file_name = self.flatfname.replace(".fits",".log")
        self.absfname = "{0}/{1}".format(self.absolute_dir, self.basefname)


class FilenamesNewflat(Filenames):
    """Filenames subclass for processed SkyMapper twilight flat images"""

    kernfpars = [(6, 0.5), (4, 1.0), (2, 2.0)]
    def __init__(self, filename, runtag=""):
        super(FilenamesNewflat, self).__init__(filename, runtag=runtag)
        self.basedir = PipelinePath.newflats
        self.relative_dir = "{0}/{1:02d}".format(self.filter, self.ccd)
        self.absolute_dir = "{0}/{1}".format(self.basedir, self.relative_dir)
        self.absfname = "{0}/{1}".format(self.absolute_dir, self.basefname)


class FilenamesNew(Filenames):
    """Filenames subclass for SkyMapper science images"""

    @staticmethod
    def relative_dir (field,filter,subfield):
        return "{0:04d}/{1}/{2:02d}".format(field, filter, subfield)
    
    @staticmethod
    def absolute_dir (field, filter, subfield):
        return "{0}/{1}".format(PipelinePath.new, FilenamesNew.relative_dir(field,filter,subfield))

    def __init__(self, filename, runtag=""):
        super(FilenamesNew, self).__init__(filename, runtag=runtag)
        self.basedir = PipelinePath.new
        self.relative_dir = FilenamesNew.relative_dir(
            self.field, self.filter, self.subfield)
        self.absolute_dir = FilenamesNew.absolute_dir(
            self.field, self.filter, self.subfield)
        self.absfname = "{0}/{1}".format(self.absolute_dir, self.basefname)


class FilenamesRef(Filenames):
    """Filenames subclass for SkyMapper galaxy reference images"""

    # FY - 2015/10/01: adding obsid, because obsseq is not unique.
    subnameregex = "ref(\d+_\d+)_([^_]+)_([uvgriz])_(\d+)_([^_]+)_(\d+)"
    subnamegroup = ["runtag", "field", "filter", "obsid","obsseq", "subfield"]
    subnameregex0 = "ref(\d+_\d+)_([^_]+)_([uvgriz])_([^_]+)_(\d+)"
    subnamegroup0 = ["runtag", "field", "filter", "obsseq", "subfield"]

    def __init__(self, filename, runtag="",subid=None):
        if subid:
            # If instead we're given a job ID
            match = re.match(self.subnameregex, subid)
            if not match:
                match = re.match(self.subnameregex0, subid)
                if match:
                    self.basefname = ''
                    self.runtag = match.group(1)
                    self.field = int(match.group(2))
                    self.filter = match.group(3)
                    self.obsid = None
                    self.obsseq = match.group(4)
                    self.subfield = int(match.group(5))
                else:
                    raise TrackableException("Bad 'subid' input to FilenamesRef")
            else:
                self.basefname = ''
                self.runtag = match.group(1)
                self.field = int(match.group(2))
                self.filter = match.group(3)
                self.obsid = match.group(4)
                self.obsseq = match.group(5)
                self.subfield = int(match.group(6))
        else:
            super(FilenamesRef, self).__init__(filename, runtag=runtag)
        self.basedir = PipelinePath.ref
        self.relative_dir = "{0:04d}/{1}/{2:02d}".format(
                self.field, self.filter, self.subfield)
        self.absolute_dir = "{0}/{1}".format(self.basedir, self.relative_dir)
        self.absfname = "{0}/{1}".format(self.absolute_dir, self.basefname)
        self.id = "ref{0}_{1:04d}_{2}_{3}_{4}_{5:02d}".format(self.runtag,
                self.field, self.filter, self.obsid, self.obsseq, self.subfield)
        self.workdir = "{0}/{1}".format(PipelinePath.scratch, self.id)
        self.log_file_name = self.id + ".log"


class FilenamesSub(Filenames):
    """Filenames subclass for SkyMapper subtraction images"""

    # FY - 2015/10/01: adding obsid, because obsseq is not unique.
    subnameregex = "sub(\d+_\d+)_([^_]+)_([uvgriz])_(\d+)_([^_]+)_(\d+)"
    subnamegroup = ["runtag", "field", "filter", "obsid", "obsseq", "subfield"]
    subnameregex0 = "sub(\d+_\d+)_([^_]+)_([uvgriz])_([^_]+)_(\d+)"
    subnamegroup0 = ["runtag", "field", "filter", "obsseq", "subfield"]

    @staticmethod
    def subid_from_fields(runtag, field, filter, obsseq, subfield,obsid=None):
       """For convenience in guessing subids"""
       if obsid:
           if obsid!='None':
               return "sub{0}_{1:04d}_{2}_{3}_{4}_{5:02d}".format(
                   runtag, field, filter, obsid, obsseq, subfield)
       return "sub{0}_{1:04d}_{2}_{3}_{4:02d}".format(
           runtag, field, filter, obsseq, subfield)
       
    def __init__(self, new=None, runtag=None, subid=None):
        """NB:  you can initialize SUB info from a Filenames instance"""
        if new:
            # Form a new Subtraction instance from NEW and possibly runtag.
            if isinstance(new, Filenames):
                # If this is a Filenames instance, just copy some stuff over
                self.basefname = new.basefname
                self.field = new.field
                self.filter = new.filter
                self.ccd = new.ccd
                self.subfield = new.subfield
                self.obsid = new.obsid
                self.obsseq = new.obsseq
                if hasattr(new, 'imgtype'):   self.imgtype = new.imgtype
                if hasattr(new, 'airmass'):   self.airmass = new.airmass
                if hasattr(new, 'exptime'):   self.exptime = new.exptime
                if hasattr(new, 'dateobs'):   self.dateobs = new.dateobs
                if hasattr(new, 'runtag'):    self.runtag = new.runtag
                if hasattr(new, 'seeing'):    self.seeing = new.seeing
                if hasattr(new, 'seewid'):    self.seewid = new.seewid
                if hasattr(new, 'elong'):     self.elong = new.elong
                if hasattr(new, 'elongwid'):  self.elongwid = new.elongwid
                if hasattr(new, 'zpmag'):     self.zpmag = new.zpmag
                if hasattr(new, 'bkgmed'):    self.bkgmed = new.bkgmed
                if hasattr(new, 'bkgsig'):    self.bkgsig = new.bkgsig
                if hasattr(new, 'maglim50'):  self.maglim50=new.maglim50
                if hasattr(new, 'maglim95'):  self.maglim95=new.maglim95
                if hasattr(new, 'skybot_cachefname'):
                    self.skybot_cachefname = new.skybot_cachefname
                if hasattr(new, 'sdss_cachefname'):
                    self.sdss_cachefname = new.sdss_cachefname
                if hasattr(new, 'program'):
                    self.program = new.program
            elif isinstance(new, str):
                # If it's just a string filename, see if it's a SUB filename;
                # if so, reroute to "subid" since filename = subid + ".fits".
                # Yes, this amounts to parsing the filename, which is evil...
                # but at least the filename convention is here in this class.
                match = re.match(self.subnameregex, new)
                match0= re.match(self.subnameregex0, new)
                if match or match0:
                    self.__init__(subid=new.replace('.fits',''))
                    return
                # Otherwise, call again and construct Filenames instance
                else:
                    self.__init__(new=Filenames(new), runtag=runtag)
                    return
            else:
                raise TrackableException("Bad 'new' input to FilenamesSub")
            # Make sure we've got the runtag!
            if runtag:  self.runtag = runtag
            else:  raise TrackableException("Can't find runtag in input")
            self.id = FilenamesSub.subid_from_fields(self.runtag,
                    self.field, self.filter, self.obsseq, self.subfield,self.obsid)
        elif subid:
            # If instead we're given a subtraction ID, parse what we can.
            # The regex below is based on the self.id format seen above.
            # We won't be able to get everything but we may not need it all.
            # This means we need to watch which FilenamesSub instances have
            # the airmass, exptime, etc. attributes defined.
            match = re.search(self.subnameregex, subid)
            if not match:
                match=re.search(self.subnameregex0, subid)
                if match:
                    self.basefname = ''
                    self.runtag = match.group(1)
                    self.field = int(match.group(2))
                    self.filter = match.group(3)
                    self.obsid  = None
                    self.obsseq = match.group(4)
                    self.subfield = int(match.group(5))
                else:
                    raise TrackableException("Bad 'subid' input to FilenamesSub")
            else:
                self.basefname = ''
                self.runtag = match.group(1)
                self.field = int(match.group(2))
                self.filter = match.group(3)
                self.obsid  = match.group(4)
                self.obsseq = match.group(5)
                self.subfield = int(match.group(6))
            self.ccd = 0
            self.id = subid
            self.airmass = 1.0
            self.exptime = 0.0
            self.dateobs = None
            self.seeing, self.seewid = 99.9, 99.9
            self.elong, self.elongwid = 99.9, 99.9
            self.program=''
        else:
            raise TrackableException("Missing parameters in FilenamesSub")

        # Obviously all Subs are IMAGETYPE='object'
        self.obstype = 'object'

        # All right!  Now that we've got all THAT stuff initialized, go ahead
        # and construct the pathnames of various SUB-related files.
        self.basedir = PipelinePath.sub
        self.relative_dir = "{0:04d}/{1}/{2:02d}".format(
                self.field, self.filter, self.subfield)
        self.absolute_dir = "{0}/{1}".format(self.basedir, self.relative_dir)
        self.workdir = "{0}/{1}".format(PipelinePath.scratch,self.id)
        self.sub_fits_name = self.id + ".fits"
        self.sub_stars_name = sextractor_output(self.sub_fits_name)
        self.sub_cands_name = self.id + ".fits.cands"
        self.hotpants_kernel_reg = self.id + ".kern.reg"
        self.sub_reg_name = ds9_region_file(self.sub_fits_name)
        self.log_file_name = self.id + ".log"


class FilenamesXref(Filenames):
    """Filenames subclass for cross-referenced detections FITS tables"""

    def __init__(self, field, subfield):
        self.field = field
        self.subfield = subfield
        self.basedir = PipelinePath.xref
        self.relative_dir = "{0:04d}/{1:02d}".format(field, subfield)
        self.absolute_dir = "{0}/{1}".format(self.basedir, self.relative_dir)
        self.basefname = "{0:04d}_{1}_xref.fits".format(field, subfield)
        self.absfname = "{0}/{1}".format(self.absolute_dir, self.basefname)
        self.calfname = "{0:04d}_{1}_xcal.fits".format(field, subfield)
        self.abscalfname = "{0}/{1}".format(self.absolute_dir, self.calfname)


class ThumbNames(object):
    """Names of thumbnails"""
    
    @staticmethod
    def thumb_dir(objname):
        return PipelinePath.thumbpath +'/'+objname

    def __init__(self, objname, subid, find = False):
        if isinstance(objname, FilenamesXsient):
            name=objname.name
        else:
            name=objname
        self.base = name + "_" + subid
        self.new = self.base + "_new.png"
        self.ref = self.base + "_ref.png"
        self.diff = self.base + "_diff.png"
        self.absolute_dir = ThumbNames.thumb_dir(name)
        self.exists = None
        if find:
            self.exists = True
            if not os.path.exists(self.absolute_dir + '/' + self.diff):
                self.exists = False
            if not os.path.exists(self.absolute_dir + '/' + self.new):
                self.exists = False
            if not os.path.exists(self.absolute_dir + '/' + self.ref):
                self.exists = False
               

class FilenamesXsient(object):
    """Separate filenames subclass for transient data and meta-data"""

    def __init__(self, name, field, subfield):
        self.name = name
        self.field = field
        self.subfield = subfield
        self.basedir = PipelinePath.lcdata
        self.relative_dir = "{0:04d}/{1:02d}".format(int(field), int(subfield))
        self.absolute_dir = "{0}/{1}".format(self.basedir, self.relative_dir)
        self.basesublcfname = "{0}_sub_lc.fits".format(self.name)
        self.basesublcglob = "*_sublc.fits".format(self.name)
        self.basenewlcfname = "{0}_new_lc.fits".format(self.name)
        self.basenewlcglob = "*_newlc.fits".format(self.name)
        self.abssublcfname = "{0}/{1}".format(
                self.absolute_dir, self.basesublcfname)
        self.absnewlcfname = "{0}/{1}".format(
                self.absolute_dir, self.basenewlcfname)
        self.basehistfname = "{0}_history.txt".format(self.name)
        self.abshistfname = "{0}/{1}".format(
                self.absolute_dir, self.basehistfname)
        self.thumb_dir = ThumbNames.thumb_dir(self.name)
        self.finderfname = "{0}/{1}_fchart.png".format(
                PipelinePath.finderpath, self.name)
        self.basexrefname = "{0:04d}_{1}_xref.fits".format(int(field), int(subfield))
        self.absxrefname = "{0}/{1}/{2}".format(PipelinePath.xref,self.relative_dir,self.basexrefname)

def Find_Image_Set(subid,localonly=False):
    """find new, ref and sub"""
    try:
        subnames = FilenamesSub(subid)
    except:
        print "Couldn't construct FilenamesSub instance for", subid
        return None
    subfname = subnames.sub_fits_name
    islocal = True
    subfile = fits_exists(subfname)
    if not subfile:
        if localonly:
            return None
        islocal = False
        subfile = fits_exists(subnames.absolute_dir + "/" + subfname)
        if not subfile:
            print "Missing sub image:",subnames.absolute_dir,subfname
            return None
    
    subhead = pyfits.getheader(subfile)
    newfname = subhead["TARGET"]
    reffname = subhead["TEMPLATE"]
    if islocal:
        newfile = fits_exists(newfname)
        reffile = fits_exists(reffname)
    else:
        newfile = fits_exists(subnames.absolute_dir + "/" + newfname)
        if not newfile:
            newfile = fits_exists(FilenamesNew.absolute_dir(
                    subnames.field, subnames.filter, subnames.subfield)+ "/" + newfname)
                                  
        reffile = fits_exists(subnames.absolute_dir + "/" + reffname)
        if not newfile:
            print "Missing new image:",newfname
            return None
        if not reffile:
            print "Missing ref image:",reffname
            return None
    
    imnames={"newname":newfile,"refname":reffile,"subname":subfile}
    return imnames
            

# ----------------------------------------------------------------------------
#              List of SkyMapper fields (mirrored in database)
# ----------------------------------------------------------------------------


def _sex2ra(ra):  return ephem.hours(ra)*180/ephem.pi
def _sex2dec(dec):  return ephem.degrees(dec)*180/ephem.pi


class SkymapperFieldCenter(AsciiRecord, CatalogEntry):
    """Represents a SkyMapper field center
    
    For text format, see $SUBPIPEHOME/etc/skymapper_field_centers.txt.

    RS 2013/06/11:  Incorporated code from my scheduling scripts, in an
    attempt to integrate and reuse some of the messier code I've written.
    """

    _ascii_separator = None
    _ascii_fields = \
        [ [ "id",   0, int,      '{0:6d}',   0    ],
          [ "ra",   1, _sex2ra,  '{0:9.5f}', None ],
          [ "dec",  2, _sex2dec, '{0:9.5f}', None ],
          [ "ebmv", 3, float,    '{0:5.3f}', None ],
          [ "exp",  4, str,      '{0}',      None ], ]

    def __str__(self):
        return "{0:4d} [{1:9.5f},{2:9.5f}]".format(self.id, self.ra, self.dec)

    def __init__(self, *args, **kwargs):
        """ScheduleField constructor"""
        # Initialize from an object
        if 'dict' in kwargs:
            self.__dict__ = kwargs['dict'].__dict__.copy()
        else:
            # First initialize as a SkymapperFieldCenter
            super(SkymapperFieldCenter, self).__init__(*args, **kwargs)
            # Set cadence and date last observed to default values
        if not hasattr(self,'cadence'): self.set_cadence()
        if not hasattr(self,'last_observed'): self.last_observed = 0.0
        if not hasattr(self,'pmult'): self.pmult = 1.0
        if hasattr(self,'exp') and (not hasattr(self,'filters') or not hasattr(self,'exptime')):
            # Then set the exposure time
            filters, exptime = [ ], { }
            for filtinfo in self.exp.split(","):
                filt, xpt = filtinfo.split(":")
                if not filt in filters:filters.append(filt)
                if filt in exptime:
                    exptime[filt].append(float(xpt))
                else:
                    exptime[filt]=[float(xpt)]
            self.filters = filters
            self.exptime = exptime
        # How long will this observation take, including overhead?
        if not hasattr(self,'obs_length'):
            self.obs_length = sum([et + Imager.overhead
                                   for f in self.filters for et in self.exptime[f]]) * ephem.second
 
    def init_ephem_body(self):
        """Sets a pyephem body at the coordinates of the field center
        
        This pre-loads use of the field center as a body with pyephem, e.g.,
        to find rise and set times.
        """
        # RA and DEC will have been retrieved in decimal degrees.
        tmpra = ephem.hours(self.ra*ephem.pi/180.0)
        tmpdec = ephem.degrees(self.dec*ephem.pi/180.0)
        coo = ephem.Equatorial(tmpra, tmpdec, epoch=2000)
        self.body = ephem.readdb("{id},f,{ra},{dec},0.0,2000.0".format(
                                 id=self.id, ra=coo.ra, dec=coo.dec))

    def pickle_safe(self):
        """Returns a picklable version of self

        Since pyephem bodies can't be pickled, this function returns a
        version of the field center with the pyephem body set to None.
        (It can easily be re-initialized when the object is reloaded.)
        """
        tmpself = copy.copy(self)
        tmpself.body = None
        return tmpself

    def set_cadence(self, cadence=4.0, cwidth=1.0):
        """Set the cadence with which this field should be observed"""
        self.cadence, self.cwidth = cadence, cwidth

    def set_obs_date(self, date):
        """Set the date on which the field was last observed."""
        self.last_observed = 1.0*date

    def calc_priority(self, date, verbose=False, min_alt=Imager.min_alt):
        """Assigns a scheduling priority to a field

        Priority is based on the time of observing and the parameters of the
        survey (cadence etc.).  0 = lowest, 100 = highest.

        RS 2015/07/13:  Remove vibrations and reference-cadence constraints,
        since these aren't really relevant when running w/multiple queues.
        Also, since the cost function is discrete like a decision tree,
        change fancy exp() factors to simple if-thens to reflect this.
        """

        # Baseline priority is to observe it (100).  Everything else whittles
        # away to *keep* us from observing bad targets.
        # RS 2015/07/13:  Add "pmult" factor to account for queue-dependent
        # things SkymapperFieldCenter doesn't necessarily need to know about.
        self.priority = 100*self.pmult
        # annoying reassigns because ephem.Observer() is impervious to copy
        site = ephem.Observer()
        site.lat, site.lon, site.horizon = sso.lat, sso.lon, sso.horizon
        site.date = date
        # ok carry on
        self.body.compute(site)
        self.az, self.alt = self.body.az, self.body.alt
        dt = site.date - self.last_observed
        if verbose:
            print "# >>> priority = {:.1f}".format(self.priority),
        # Is the field above the horizon (or a reasonable airmass)?
        if np.degrees(self.alt) < min_alt:
            self.priority = 0.0
            if verbose:
                print "# >> priority = {0:.1f} (airmass, alt:{1:.1f}) ".format(self.priority,np.degrees(self.alt)),
        # Are we too early?
        if dt < self.cadence:
            self.priority = 0.0
            if verbose:
                print "# >> priority = {0:.1f} (cadence too small, dt:{1:.1f}) ".format(self.priority,dt),
        # Bump priority up if the field is badly overdue
        if dt > 1.5*self.cadence and dt < 30.0: #but not if it's way too overdue...
            self.priority *= 2.0
            if verbose: 
                print "# >> {:.1f} (deadman) ".format(self.priority)
        # Account for any extraneous priority multipliers
        if verbose:
            print "# >> {:.1f}".format(self.priority)
        # Return final value
        return self.priority


SMfield_list = SkymapperFieldCenter.read_ascii_file(
        PipelinePath.etc + "/skymapper_field_centers.txt")

def SMfield_lookup(ra=None, dec=None, id=None, maxdist=0.5, maxra=0.3):
    """Find SkyMapper field corresponding to (ra, dec) in decimal degrees.
    
    RS 2012/05/25:  A mirror of its namesake in mydjango.jobstats.STAP_API.
    I put it here so that Fang could get started while I continued to work
    on the Django stuff.  It does something similar, but reads the fields
    from the file in $SUBPIPEHOME/etc/ and returns a SkymapperFieldCenter.

    RS 2013/10/30:  Made "maxdist" a keyword argument, since this code is
    used in the pipeline to match NEWs to REFs (default case), but also in
    vetting code to find a given (RA, DEC) in a field, which we might want
    to be less restrictive if a candidate is near the field's edge.

    RS 2014/01/10:  Did same with maxra.
    """

    # If a field ID is provided, just use that to look it up.
    if id is not None:
        checklist = [smfield for smfield in SMfield_list if smfield.id == id]
        if len(checklist) == 0:  return None
        return checklist[0]
    if ra is None or dec is None:
        print "SMfield_lookup:  Need either field ID or field (RA,DEC)!"
        return None

    # Calculate the arc distance to all field centers; return the closest.
    # def hpos(ra, dec):  return ephem.degrees(str(ra)), ephem.degrees(str(dec))
    def hpos(ra, dec):
        if isinstance(ra, str):  ephra = ephem.hours(ra)
        else:  ephra = ephem.degrees(str(ra))
        ephdec = ephem.degrees(str(dec))
        return ephra, ephdec
    arcdist = np.array([ ephem.separation(hpos(ra, dec), hpos(f.ra, f.dec))
                         * 180/ephem.pi for f in SMfield_list ])
    mindist=arcdist[arcdist.argmin()]
    if mindist > maxdist:
        return SMfield_list[0]
    elif mindist >maxra:
        fid=arcdist.argmin()
        radist=ephem.separation(hpos(ra, SMfield_list[fid].dec), hpos(SMfield_list[fid].ra, SMfield_list[fid].dec))*180/ephem.pi
        if radist <maxra:
            return SMfield_list[fid]
        #offset is too large for match to be meaningful
        return SMfield_list[0]
    else:
        return SMfield_list[arcdist.argmin()]


# ----------------------------------------------------------------------------
#                         Scheduler log file format
# ----------------------------------------------------------------------------


class SchedulerLogEntry(CatalogEntry):
    """Represents a single line in the scheduler log

    This is a record to hold entries in the scheduler log, which describe
    the status of an observation which was requested during the night.
    """

    def __init__(self, logline):
        """Parses a line from the scheduler log using regular expressions"""

        # default for if something fails
        self.obsdate = None
        self.obsseq = None
        self.program = None
        self.obsid = None
        self.imgid = None
        self.filter = None
        self.rotskypa = None
        self.ra = None
        self.dec = None
        self.ra_str, self.dec_str = "", ""

        # here's the line with what the scheduler told TAROS to do
        mm = re.search("(\S+ \S+) .*Ob: (3PS|BAD); id (\d+); .* filter (.); rot (\d+);"
                       ".*RA/dec. (\S+) / (\S+) \((\S+) / (\S+)\)", logline)
        if mm is None:  return
    
        # otherwise, actually start parsing things
        # WARNING:  ephem 3.7.4.1 doesn't like hyphens in dates
        strdate = mm.group(1)[:19].replace('-','/')
        self.obsdate = ephem.Date(strdate)
        self.obsseq = strdate.replace(' ','T')
        self.program = mm.group(2)
        self.obsid = int(mm.group(3))
        self.imgid = "{0}_{1}".format(self.obsid, self.obsseq[:-3])
        self.filter = mm.group(4).strip()
        self.rotskypa = float(mm.group(5))
        self.ra = ephem.hours(mm.group(6))
        self.dec = ephem.degrees(mm.group(7))
        ra_str, dec_str = str(self.ra), str(self.dec)
        self.ra *= 180/ephem.pi
        self.dec *= 180/ephem.pi
        # force format of RA and DEC to conform to normal FITS headers
        # horribly clunky, but it works
        ra_groups = [float(f) for f in ra_str.split(":")]
        ra  = "{0:02.0f}:{1:02.0f}:{2:04.1f}".format(*ra_groups)
        dec_groups = [float(f) for f in dec_str.split(":")]
        if self.dec > 0:
            dec_groups = ['+'] + dec_groups
        else:
            dec_groups = ['-'] + dec_groups
            dec_groups[1] = abs(dec_groups[1])
        dec = "{0}{1:02.0f}:{2:02.0f}:{3:02.0f}".format(*dec_groups)
        self.ra_str, self.dec_str = ra_str, dec_str


# ----------------------------------------------------------------------------
#                         Skymapper observation id 
# ----------------------------------------------------------------------------

ID_SEQ_NO_BITS   =  9 
ID_OB_SEQ_BITS   =  5 
ID_DR_SER_BITS   =  1 
ID_FIELD_ID_BITS = 13 
ID_OB_TYPE_BITS  =  4 
ID_LOW_BITS = ( ID_SEQ_NO_BITS + ID_OB_SEQ_BITS + ID_DR_SER_BITS +
                ID_FIELD_ID_BITS )

TID_SEQ_NO_BITS    =  9 
TID_OBSID_BITS     = 12 
TID_SURVEY_ID_BITS =  7 
TID_OB_TYPE_BITS   =  4

OBSERVATION_TYPE = {"MS":1,
                    "TFF":2,
                    "3PS":3,
                    "STD":4,
                    "5SS":5,
                    "BAD":6,
                    "BIAS":7}

def get_ob_name(ob_type):
    ob_name=''
    for key, value in OBSERVATION_TYPE.iteritems():
        if value==ob_type:
            ob_name=key
            break
    return ob_name

def decode_obs_id(encoded_ids):
    encoded_ids=int(encoded_ids)
    if  (( encoded_ids >> ID_LOW_BITS ) & ( 2 ** ID_OB_TYPE_BITS - 1 )==OBSERVATION_TYPE["3PS"]  
         or (encoded_ids >> ID_LOW_BITS ) & ( 2 ** ID_OB_TYPE_BITS - 1 )==OBSERVATION_TYPE["BAD"]):   
        seq_no_3PS = encoded_ids & ( 2 ** TID_SEQ_NO_BITS   - 1 )   # Get sequence number
        encoded_ids >>= TID_SEQ_NO_BITS    
        obsid = encoded_ids & ( 2 ** TID_OBSID_BITS   - 1 )   # Get obsid number
        encoded_ids >>= TID_OBSID_BITS
        survey_id = encoded_ids & ( 2 ** TID_SURVEY_ID_BITS - 1 )   # Get survey number
        encoded_ids >>= TID_SURVEY_ID_BITS 
        ob_type  = encoded_ids & ( 2 ** TID_OB_TYPE_BITS  - 1 )   # Get observation type
        return get_ob_name(ob_type), survey_id, obsid, seq_no_3PS
    else:
        encoded_ids >>= ID_SEQ_NO_BITS                                # Scrap observation id sequence number
        ob_seq_no = encoded_ids & ( 2 ** ID_OB_SEQ_BITS   - 1 )   # Get sequence (filter) number
        encoded_ids >>= ID_OB_SEQ_BITS         
        dr_ser_no = encoded_ids & ( 2 ** ID_DR_SER_BITS   - 1 )   # Get data release or series number
        encoded_ids >>= ID_DR_SER_BITS        
        field_id  = encoded_ids & ( 2 ** ID_FIELD_ID_BITS - 1 )   # Get field id
        encoded_ids >>= ID_FIELD_ID_BITS         
        ob_type   = encoded_ids & ( 2 ** ID_OB_TYPE_BITS  - 1 )   # Get observation type
        return get_ob_name(ob_type), field_id, dr_ser_no, ob_seq_no

# ----------------------------------------------------------------------------
