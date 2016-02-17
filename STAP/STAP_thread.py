import os
import sys
from Utils import Pipeline as Pipeline
from Utils import Constants
import STAP
import glob
import subprocess

class STAP_thread(Pipeline.Workflow):
    """Runs a standard production-style subtraction thread.

    RS:  Will write proper comments later.
    """

    def _setup_workflow(self, newnames, wrefnames, runtag):

        # These things all belong to stuff derived from Constants
        newimg     = newnames.absfname
        wrefimg    = wrefnames.absfname
        subnames   = Constants.FilenamesSub(newimg, runtag)
        wnewimg    = os.path.basename(newimg.replace(".fits","_wcs.fits"))
        wrefimg_zipped = wrefimg
        wrefimg=wrefimg.replace('.gz','')
        rwrefimg   = os.path.basename(wrefimg.replace(".fits",'_'+subnames.obsseq+"_remap.fits"))
        subimg     = subnames.sub_fits_name
        wnewstars  = Constants.sextractor_output(wnewimg)
        rwrefstars = Constants.sextractor_output(rwrefimg)
        substars   = subnames.sub_stars_name
        subcands   = subnames.sub_cands_name
        subreg     = subnames.sub_reg_name
        subkern    = subnames.hotpants_kernel_reg
        logfile    = subnames.log_file_name
        skybotfn   = subnames.skybot_cachefname
        sdssfn     = subnames.sdss_cachefname
        xrefnames  = Constants.FilenamesXref(subnames.field, subnames.subfield)
        xrefcat    = xrefnames.absfname
        calstars   = xrefnames.abscalfname
        configjunk = Constants.PipelinePath.scratchetc_files
        newmask    = newimg.replace('.fits','_mask.fits')
        refmask    = wrefimg.replace('_wcs.fits','_mask.fits')
        submask    = subimg.replace('.fits','_mask.fits')
        newmask_zipped=newmask+'.gz'
        refmask_zipped=refmask+'.gz'

        # Save some of these so the cleanup phase can refer to them
        self.newnames = newnames
        self.refnames = wrefnames
        self.subnames = subnames
        self.wnewimg = wnewimg
        self.rwrefimg = rwrefimg
        self.subimg = subimg

        # Get base filenames of first stage inputs; these will be in workdir
        lnewimg = os.path.basename(newimg)
        lwrefimg = os.path.basename(wrefimg)
        lxrefcat = os.path.basename(xrefcat)
        lcalstars = os.path.basename(calstars)
        lsubstars = os.path.basename(substars)
        lwrefimg_zipped=os.path.basename(wrefimg_zipped)
        lnewmask   = os.path.basename(newmask)
        lrefmask   = os.path.basename(refmask)
        lnewmask_zipped=lnewmask+'.gz'
        lrefmask_zipped=lrefmask+'.gz'
        rwrefmask  = rwrefimg.replace('.fits','_mask.fits')

        # Define _copy_over and _copy_back; these will all be absolute paths
        self._copy_over = [newimg, wrefimg_zipped, xrefcat, calstars,
                           skybotfn, sdssfn, 
                           newmask_zipped, refmask_zipped] + configjunk
        self._copy_back = [subnames.absolute_dir + "/" + fn for fn in
                           (wnewstars, rwrefstars, substars, subcands,
                            subreg, subkern, logfile)] + [xrefcat, calstars]

        # RS 2012/08/29:  Don't forget previous light curve info
        lcfnames = [ ]
        for filter in Constants.Imager.filters:
            histnewnames = Constants.FilenamesNew(newimg)
            histnewnames.filter = filter
            histsubnames = Constants.FilenamesSub(histnewnames, runtag)
            srcpath = histsubnames.absolute_dir
            fnames = [os.path.basename(fn)
                      for fn in glob.glob(srcpath + "/sub*.fits.stars")]
            #         for fn in glob.glob(srcpath + "/*.fits.cands")]
            lcfnames += fnames
            self._copy_over += [srcpath + "/" + fn for fn in fnames]
        lcfnames.append(lsubstars)
        #glob_xsient= Constants.FilenamesXsient(
        #        "*", subnames.field, subnames.ccd)
        #self._copy_over += [fn for fn in glob.glob(glob_xsient.abshistfname)]

        # List of pipeline stages to execute
        self._stagelist = \
        [
            ["unzip_REF", STAP.unzip, [lwrefimg_zipped]],
            ["unzip_NEWMASK", STAP.unzip, [lnewmask_zipped]],
            ["unzip_REFMASK", STAP.unzip, [lrefmask_zipped]],
            ["WCS",      STAP.WCS, [lnewimg, wnewimg]],
            ["remap",    STAP.remap, [wnewimg, lwrefimg, rwrefimg], 
                { "module": "swarp", "weightmap": lrefmask,
                  "weightmap_out":rwrefmask }],
            ["maskcombine", STAP.maskcombine, [[lnewmask,rwrefmask]],
                { "outname":submask, "combinetype":"AND" }],            
            ["SEx_NEW",  STAP.SEx, [wnewimg], 
                { "do_seeing": True, "do_apcor": True }],
            ["SEx_REF",  STAP.SEx, [rwrefimg], 
                { "do_seeing": True, "do_apcor": True }],
            ["zp_NEW",   STAP.zeropoint, [wnewstars], 
                { "calstarsfname": lcalstars, "imgfname": wnewimg }],
            ["hotpants", STAP.hotpants, [wnewimg, rwrefimg, subimg], 
                { "sexstamps" : False }],
            ["SEx_SUB",  STAP.SEx, [subimg], 
                { "SEx_switches": "-WEIGHT_TYPE MAP_WEIGHT "
                                  "-WEIGHT_THRESH 0.0001 -WEIGHT_IMAGE "
                                  + submask }],
            ["classify", STAP.classify, [wnewimg, rwrefimg, subimg]],
            ["xref",     STAP.xref, [lxrefcat],
        	{ "candsfnames": [subcands], "lcfnames": lcfnames,
               	"skybot_cachefname": skybotfn, "sdss_cachefname": sdssfn,
		"replace": True }],
            ["makethumbs",   STAP.makethumbs,[[],lxrefcat,]],
        ]

    def _cleanup_workflow(self):

        # These things all belong to stuff derived from Constants
        subnames = self.subnames
        wnewimg  = self.wnewimg
        rwrefimg = self.rwrefimg
        subimg   = self.subimg

        # Convert all the FITS images to integer format in prep for copying,
        # and try to compress to gzip.  If the compression stage fails then
        # just copy back the original uncompressed file (for debugging).
        for fn in (wnewimg, rwrefimg, subimg):
            if not os.path.exists(fn):  continue
            absfn = subnames.absolute_dir + "/" + fn
            if fn == subimg or fn == rwrefimg:
                fnmask = fn.replace('.fits','_mask.fits')
            else:
                fnmask = fn.replace('_wcs.fits','_mask.fits')
            if os.path.exists(fnmask):mask=fnmask
            else: mask=None
            try:
                print "Reintegerizing and compressing", fn
                STAP.roundint(fn, mask=mask,maskflag=0,gzip=True)
                self._copy_back += [absfn + ".gz"]
            except:
                self._copy_back += [absfn]
            if os.path.exists(fnmask):
                absfn = subnames.absolute_dir + "/" + fnmask
                subprocess.call(["gzip",fnmask])
                self._copy_back += [absfn + ".gz"]

        # Did we find any transients?  If so, copy light curves etc. back.
        xrefout = [o["results"] for o in self.output["stages"]
                   if o["name"] == "xref" and "results" in o]
        if len(xrefout) == 1:
            print "Found", len(xrefout[0]), "transients"
            for id, xdict in xrefout[0].items():
                xsientnames = Constants.FilenamesXsient(
                        xdict["name"], xdict["field"], xdict["subfield"])
                print "Copying back", xsientnames.abssublcfname
                self._copy_back.append(xsientnames.abssublcfname)
                #print "Copying back", xsientnames.abshistfname
                #self._copy_back.append(xsientnames.abshistfname)
                #find thumbnails
                thumbs=glob.glob("*{0}*.png".format(xsientnames.name))
                for thumb in thumbs:
                    self._copy_back.append("{0}/{1}".format(xsientnames.thumb_dir,thumb))


class STAP_thread_runref(Pipeline.Workflow):
    """ Runs a production thread designed to initiate a reference cache

    RS:  Will write proper comments later.
    """

    def _setup_workflow(self, newimg):

        # These things all belong to stuff derived from Constants
        newnames = Constants.FilenamesRef(newimg)
        xrefnames  = Constants.FilenamesXref(newnames.field, newnames.subfield)
        calstars   = xrefnames.abscalfname
        configjunk = Constants.PipelinePath.scratchetc_files
        newmask    = newimg.replace('.fits','_mask.fits')
        newmask_zipped=newmask + '.gz'

        # Get base filenames of first stage inputs; these will be in workdir
        lnewimg   = os.path.basename(newimg)
        lcalstars = os.path.basename(calstars)
        wnewimg   = lnewimg.replace(".fits","_wcs.fits")
        wnewstars = wnewimg.replace(".fits",".fits.stars")
        logfile   = self.logfile #newnames.log_file_name
        lnewmask   = os.path.basename(newmask)
        lnewmask_zipped=lnewmask+'.gz'
 
        # Save some of these for later in case we need them
        self.copyback_dir = newnames.absolute_dir
        self.wnewimg = wnewimg
        self.wnewstars = wnewstars

        # Define _copy_over and _copy_back; these will all be absolute paths
        self._copy_over = [newimg, calstars, newmask_zipped] + configjunk
        self._copy_back = [self.copyback_dir + "/" + fn for fn in
                           (wnewstars, logfile)] + [calstars]

        # List of pipeline stages to execute
        self._stagelist = \
        [
            ["unzip_NEWMASK", STAP.unzip, [lnewmask_zipped]],
            ["WCS", STAP.WCS, [lnewimg, wnewimg]],
            ["SEx", STAP.SEx, [wnewimg], { "rerun": True }],
            ["zp",  STAP.zeropoint, [wnewstars],
                { "calstarsfname": lcalstars, "imgfname": wnewimg } ],
        ]

    def _cleanup_workflow(self):

        # These things all belong to stuff derived from Constants
        wnewimg = self.wnewimg
        wnewstars = self.wnewstars

        # Convert all the FITS images to integer format in prep for copying,
        # and try to compress to gzip.  If the compression stage fails then
        # just copy back the original uncompressed file (for debugging).
        absfn = self.copyback_dir + "/" + wnewimg
        try:
            print "Reintegerizing and compressing", wnewimg
            STAP.roundint(wnewimg, gzip=True)
            self._copy_back += [absfn + ".gz"]
        except:
            self._copy_back += [absfn]
                        
        fnmask = wnewimg.replace('_wcs.fits','_mask.fits')
        if os.path.exists(fnmask):
            absfn = self.copyback_dir + "/" + fnmask
            subprocess.call(["gzip",fnmask])
            self._copy_back += [absfn + ".gz"]


class STAP_thread_runflats(Pipeline.Workflow):
    """ Runs a production thread designed to create flat fields & weight maps

    RS:  Will write proper comments later.
    """

    def _setup_workflow(self, imglist):

        # These things all belong to stuff derived from Constants
        imgnamelist = [Constants.FilenamesNewflat(fn) for fn in imglist]
        configjunk = Constants.PipelinePath.scratchetc_files

        # Give the superflat the obsseq of the *oldest* constituent flat name.
        # This will make it easy to tell whether it has gotten stale.
        imgnamelist.sort(key=lambda x: x.obsseq)
        ffnames = Constants.FilenamesFlat(imgnamelist[0].absfname)
        ffpixname = ffnames.flatfname
        wtmapname = ffnames.wtmapfname
        logfile = ffnames.log_file_name
        limglist = [img.basefname for img in imgnamelist]

        #find badpix mask
        maskpath=Constants.PipelinePath.masks+'/%s/%02d'%(ffnames.filter,ffnames.ccd)
        maskfiles=glob.glob(maskpath+'/badpix*.fits')
        if len(maskfiles)>0:
            maskfiles.sort()
            maskname=maskfiles[-1]
            domask=True
        else:
            domask=False
        
        # Define _copy_over and _copy_back; these will all be absolute paths
        self._copy_over = imglist + configjunk
        if domask: self._copy_over += [maskname]
        self._copy_back = [ffnames.absolute_dir + "/" + fn
                           for fn in [ffpixname]]  #, logfile

        # List of pipeline stages to execute
        if domask:
            self._stagelist = \
                [
                ["flatcombine", STAP.flatcombine,
                 [limglist, ffpixname], { "mask": os.path.basename(maskname)}],
                ]
        else:
            self._stagelist = \
                [
                ["flatcombine", STAP.flatcombine,
                 [limglist, ffpixname]],
                ]
    

class STAP_quickmag(Pipeline.Workflow):
    """Runs a quick crappy subtraction against DSS to get a magnitude.

    RS:  Will write proper comments later.
    """

    def _setup_workflow(self, newimg, ra, dec):

        # These things all belong to stuff derived from Constants
        staging = Constants.PipelinePath.data + "/quickmag"
        lnewimg = os.path.basename(newimg)
        wnewimg = lnewimg.replace(".fits","_wcs.fits")
        lrefimg = "dssimg.fits"
        wrefimg = lrefimg.replace(".fits","_wcs.fits")
        rwrefimg = wrefimg.replace(".fits","_remap.fits")
        subimg = "subimg.fits"
        lnewstars = lnewimg + ".stars"
        wrefstars = wrefimg + ".stars"
        substars = subimg + ".stars"
        configjunk = Constants.PipelinePath.scratchetc_files

        # Figure out which USNO catalog we need to match against
        meta = Constants.Filenames(newimg)
        if meta.filter in ['g']:
            dsscat, usnomag = "poss2ukstu_blue", "B1mag"
        elif meta.filter in ['r', 'i']:
            dsscat, usnomag = "poss2ukstu_red", "R1mag"

        # Define _copy_over and _copy_back; these will all be absolute paths
        self._copy_over = [newimg] + configjunk
        self._copy_back = [staging + "/" + fn for fn in (wnewimg, rwrefimg)]

        # List of pipeline stages to execute
        self._stagelist = \
        [
            ["SEx_NEW",  STAP.SEx, [lnewimg], { "rerun": False }],
            ["WCS_new",  STAP.WCS, [lnewimg, wnewimg], { "astronet": True }],
            ["getDSS",   STAP.getDSS, [ra, dec, lrefimg], { "survey": dsscat }],
            ["WCS_ref",  STAP.WCS, [lrefimg, wrefimg], { "astronet": True }],
            ["remap",    STAP.remap, [wnewimg, wrefimg, rwrefimg], { "module": "wcsremap" }],
            ["SEx_REF" , STAP.SEx, [rwrefimg], { "rerun": False }],
            ["hotpants", STAP.hotpants, [lnewimg, rwrefimg, subimg], { "sexstamps": True }],
            ["SEx_SUB",  STAP.SEx, [subimg]],
            ["getmag",   STAP.getmag, [lnewstars, substars, ra, dec], { "usnomag": usnomag }]
        ]
