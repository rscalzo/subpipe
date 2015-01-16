#!/usr/bin/env python

# ============================================================================
# RS 2014/02/18:  SkyMapper Subtraction Pipeline Quick-Fix
# ----------------------------------------------------------------------------
# This script is meant to re-run only the zeropoint step based on existing
# main pipeline output (to avoid having to re-run stacks of subtractions).
# It takes as input FITS files from a pre-existing pipeline run or runs
# and re-builds the survey history from them.
# ============================================================================


# Modules to include

import time
import sys
import os
import re
import glob
import ephem
import pyfits
import argparse
import datetime
import subprocess
import numpy as np
from STAP.STAP_comm import psmemchk
from STAP.STAP_tools.headerprop import headerprop
from Utils import Constants
from Utils import Pipeline as Pipeline
from mydjango.jobstats.STAP_API import register_job, register_run
# def register_run(*args, **kwargs):  pass
# def register_job(*args, **kwargs):  pass


# ----------------------------------------------------------------------------
#                            The main routine
# ----------------------------------------------------------------------------


def main():
    """Top level ("master process") of the pipeline"""

    # Get command-line arguments
    today_utc = datetime.datetime.utcnow().strftime("%Y-%m-%d")
    parser = argparse.ArgumentParser(description=
            "Correct headers of incoming images based on scheduler logs")
    parser.add_argument('--utdate', default=None,
            help="path to images (default: %(default)ss)'")
    parser.add_argument('--nproc', type=int, default=32,
            help="number of processors to use (default: 32)")
    args = parser.parse_args()
    print "Using pyfits version", pyfits.__version__

    # Initialize a bunch of important directories that may not exist
    etc = Constants.PipelinePath.etc
    scratch = Constants.PipelinePath.scratch
    scratchetc = Constants.PipelinePath.scratchetc
    subprocess.call("mkdir -p {0}".format(scratch), shell=True)
    subprocess.call("rm -rf {0}/sub*".format(scratch), shell=True)
    subprocess.call("rm -rf {0}/Sky*".format(scratch), shell=True)
    subprocess.call("rm -rf {0}".format(scratchetc), shell=True)
    subprocess.call("cp -r {0} {1}".format(etc,scratchetc), shell=True)

    # Pause lock file name
    pauselockfname = Constants.PipelinePath.pauselockfname

    # Sort out glob strings for *local* dates for today UT.
    if args.utdate in [None, 'all']:
        today = ephem.now()
    else:
        today = ephem.Date(args.utdate.replace("-","/"))
    tomorrow = ephem.Date(today + 1.0)
    # today, tomorrow = ephem.Date(ephem.now() - 1.0), ephem.now()
    today_str = today.datetime().strftime("%Y-%m-%d")
    if args.utdate == 'all':
        ngstr = ['', 'grzbdlx', 'grzbdlx']
    else:
        ngstr = [today.datetime().strftime("%Y-%m-%dT1"),
                 today.datetime().strftime("%Y-%m-%dT2"),
                 tomorrow.datetime().strftime("%Y-%m-%dT0")]
    print "ngstr =", ngstr

    # Set up the queue for pipeline processes
    Q = Pipeline.LowMemQueue(name="Main Pipeline Queue", processors=32)

    start_run_time = time.time()        # start time for run, use as run ID
    register_run(start_run_time)        # register run in Django
    os.chdir(scratch)                   # scratch should be our working path

    donefnames = [ ]                    # names of files already processed
    fieldccdlock = { }                  # Field+CCD combinations in progress
    loop_delay = 30                     # Minimum duration of loop in seconds
    pause_delay = 600                   # Length of pause cycle in seconds

    # Get unique identifier for this run
    runtag = Constants.get_runtag(start_run_time)

    def locktag(thing):
        return "{0}_{1:02d}".format(thing.field, thing.subfield)

    # The inputs will be the .cands files.
    CPP_base_glob = Constants.PipelinePath.sub + "/*/*/*/*"
    freshfnames = [f for f in
                   glob.glob(CPP_base_glob + ngstr[0] + "*.fits.cands") +
                   glob.glob(CPP_base_glob + ngstr[1] + "*.fits.cands") +
                   glob.glob(CPP_base_glob + ngstr[2] + "*.fits.cands")
                   if f not in donefnames]
    # Which will need to be run in chronological order.
    obsseq = []
    for ffname in freshfnames:
        obsseq.append(os.path.basename(ffname).split('_')[-2])
    sortf = np.array(obsseq).argsort()
    freshfnames = list(np.array(freshfnames)[sortf]) # [:1]

    # The main event loop
    while 1:

        start_loop_time = time.time()
        print "subpipe_master.py:  running main loop at local time ",
        print time.asctime(time.localtime(start_loop_time))

        # Which ones haven't been done yet?
        freshfnames = [f for f in freshfnames if f not in donefnames]
        # Check memory use
        memtot, pmemtot = psmemchk('python subpipe_master.py')
        print "Running {0:d}/{1:d} SUBs now; {2:d} NEW images waiting to run;"\
                .format(len(fieldccdlock),Q.N_processors,len(freshfnames)),
        print "memory use {0:.1f} MB ({1:.1f}%)".format(memtot, pmemtot)
        sys.stdout.flush()

        # If we're done, quit!
        if len(freshfnames) + len(fieldccdlock) == 0:
            print "Ran out of stuff to do; exiting normally."
            return

        # Start as many new jobs as we have empty slots
        for subcandfname in freshfnames:

            # Check to see if we've paused.  If we have, break out of this
            # loop but continue to retrieve results for finished jobs.
            if os.path.exists(pauselockfname):
                print "subpipe_master.py:  Someone hit the pause switch"
                print "Not starting any new jobs, but collecting old ones."
                sys.stdout.flush()
                break

            # Get the metadata in two stages:  first use a regex parse to get
            # the original sub's runtag, then the rest from the FITS header.
            subid = os.path.basename(subcandfname).replace(".fits.cands","")
            sub = Constants.FilenamesSub(subid)
            subruntag = sub.runtag
            # RS 2014/02/21:  Some issues with corrupt FITS files...
            try:
                sub = Constants.FilenamesSub(subcandfname, runtag=subruntag)
            except IOError:
                print "FITS error -- can't even set up {0}!  Skipping.".format(
                        subid)
                donefnames.append(subcandfname)
                continue

            # If there's no room left, don't submit any more
            if len(fieldccdlock) >= Q.N_processors:  break

            # Don't process the same field+CCD combination at the same time,
            # it'll screw up candidate cross-referencing.
            fieldccdtag = locktag(sub)
            if fieldccdtag in fieldccdlock:  continue

            # Set up the workflow
            workflow = STAP_thread(name=sub.id+"_rezpcal",
                                   workdir=sub.workdir,
                                   logfile=sub.id+"_rezpcal.log")
            Qargs = (workflow, subcandfname, subruntag)

            # Try to submit the job to the queue.
            status = Q.start_job(*Qargs)
            if status == 0:
                print "Submitted job {0} to queue".format(workflow.name)
                fieldccdlock[fieldccdtag] = sub
                donefnames.append(subcandfname)
            elif status == 1:
                print "Couldn't submit job ", workflow.name,
                print "; will try again later"
            elif status == -1:
                print "Fatal error setting up for job", workflow.name
                donefnames.append(subcandfname)
            sys.stdout.flush()

        # Wait for a bit
        time.sleep(loop_delay/2)

        # Deal with output of terminated processes
        while not Q.empty():
            wf = Q.finish_job()
            job_out = wf.output
            # print "Retrieved output from job", job_out["name"]

            # Figure out what happened to the job
            print job_out["name"],
            if "exception" in job_out["stages"][-1]:
                print "failed in stage", job_out["stages"][-1]["name"],
                print "from an uncaught exception:"
                print job_out["stages"][-1]["exception"]
            elif "results" in job_out["stages"][-1]:
                print "has successfully completed"
            else:
                print "probably got killed before finishing"

            # Register the job along with its result, then release the lock
            # on this field+CCD combination.
            for fieldccdtag in fieldccdlock:
                sub = fieldccdlock[fieldccdtag]
                if re.match(sub.id, job_out["name"]):
                    register_job(sub, job_out)
                    del fieldccdlock[fieldccdtag]
                    break

        end_loop_time = time.time()

        # If we're still paused, hang out for a few minutes
        if os.path.exists(pauselockfname):
            pause_wait = pause_delay - (end_loop_time - start_loop_time)
            if pause_wait > 0 and Q.empty():
                print "Waiting", pause_delay, "seconds while paused",
                print "before checking", pauselockfname
                sys.stdout.flush()
                time.sleep(pause_wait)

        # If necessary, hang out till the end of this delay cycle
        sys.stdout.flush()
        if len(fieldccdlock) < Q.N_processors and len(freshfnames) > 0:
            continue
        loop_wait = loop_delay-(end_loop_time-start_loop_time)
        if loop_wait > 0:  time.sleep(loop_wait)


# ----------------------------------------------------------------------------
#                     Definition of the workflow class
# ----------------------------------------------------------------------------


import os
from Utils import Pipeline as Pipeline
from Utils import Constants
import STAP
import glob


class STAP_thread(Pipeline.Workflow):
    """Runs the thread to re-run the xref step."""

    def fixhdrs(self, newstars, substars, subcands, subimg):
        """Fixes headers of files on ramdisk
        
        This locates mission-critical header information and propagates it
        between files.  Maybe all this junk should be in SQL.  Argh.
        """
        # First what STAP.SEx would normally do
        # WARNING:  This won't work if we failed in step classify.
        # We need to make STAP.classify or STAP.makethumbs less brittle.
        sexkw = ['FILTNAME', 'FIELD_ID', 'SUBFIELD',
                 'IMAGEID',  'EXPTIME',  'DATE-OBS', 'RUNTAG']
        if os.path.exists(subcands):
            print "Propagating from", subcands, "to", substars
            headerprop(subcands, substars, sexkw)
        # Then pass on the new zeropoint information
        zpckw = ['ZPMAG',    'ZPMAGERR', 'ZPCALSRC', 'ZPCALTYP',
                 'FAINTLIM', 'FLXLIM50', 'FLXLIM95', 'MAGLIM50', 'MAGLIM95']
        for fn in (substars, subimg, subcands):
            if os.path.exists(fn):
                print "Propagating from", newstars, "to", fn
                headerprop(newstars, fn, zpckw)

    def _setup_workflow(self, subcands, runtag):
        """Set up lots of annoying bookkeeping for the run"""

        # These things all belong to stuff derived from Constants
        # Find_Image_Set() appears only to work when the images themselves
        # exist, and when we're getting things from the SUB header.
        # So we have to work around that.
        subnames   = Constants.FilenamesSub(subcands, runtag=runtag)
        subhdr     = pyfits.getheader(subcands)
        subpath    = subnames.absolute_dir + '/'
        subimggz   = subpath + subnames.sub_fits_name + '.gz'
        newimggz   = subpath + subhdr['TARGET'] + '.gz'
        refimggz   = subpath + subhdr['TEMPLATE'] + '.gz'
        xrefnames  = Constants.FilenamesXref(subnames.field, subnames.ccd)
        xrefcat    = xrefnames.absfname
        calstars   = xrefnames.abscalfname
        configjunk = Constants.PipelinePath.scratchetc_files
        logfile    = subnames.absolute_dir + '/' + self.logfile

        # Don't forget the .stars files
        subimg     = subimggz.replace('.gz','')
        newimg     = newimggz.replace('.gz','')
        refimg     = refimggz.replace('.gz','')
        substars   = Constants.sextractor_output(subimg)
        newstars   = Constants.sextractor_output(newimg)
        refstars   = Constants.sextractor_output(refimg)

        # Get base filenames of first stage inputs; these will be in workdir
        # Man this is a lot of bookkeeping.
        lsubimggz  = os.path.basename(subimggz)
        lnewimggz  = os.path.basename(newimggz)
        lrefimggz  = os.path.basename(refimggz)
        lsubstars  = os.path.basename(substars)
        lnewstars  = os.path.basename(newstars)
        lrefstars  = os.path.basename(refstars)
        lsubcands  = os.path.basename(subcands)
        lxrefcat   = os.path.basename(xrefcat)
        lcalstars  = os.path.basename(calstars)
        lsubimg    = lsubimggz.replace('.gz','')
        lnewimg    = lnewimggz.replace('.gz','')
        lrefimg    = lrefimggz.replace('.gz','')

        # Save some of these so the cleanup phase can refer to them
        self.subnames = subnames
        self.substars = substars
        self.newstars = newstars
        self.refstars = refstars

        # This part is made more complicated by the fact that we may have
        # blown away some of the old images, and we should do what we can
        # even if we don't have those around.  But we're guaranteed to have
        # the survey history.
        self._stagelist = [ ]
        self._copy_back = [calstars, substars, newstars, refstars, subcands]
        self._copy_over = self._copy_back + configjunk
        self._copy_back.append(logfile)

        # If the image files exist, copy them over and arrange to reduce.
        # This is kind of shitty code, it's annoying.
        round_stages = [ ]
        if os.path.exists(subimggz):
            self._copy_over.append(subimggz)
            self._copy_back.append(subimggz)
            self._stagelist.append(["unzip_SUB", STAP.unzip, [lsubimggz]])
            round_stages.append(
                ["round_SUB", STAP.roundint, [lsubimg], { 'gzip': True }])
        if os.path.exists(refimggz):
            self._copy_over.append(refimggz)
            self._copy_back.append(refimggz)
            self._stagelist.append(["unzip_REF", STAP.unzip, [lrefimggz]])
            round_stages.append(
                ["round_REF", STAP.roundint, [lrefimg], { 'gzip': True }])
        if os.path.exists(newimggz):
            self._copy_over.append(newimggz)
            self._copy_back.append(newimggz)
            self._stagelist.append(["unzip_NEW", STAP.unzip, [lnewimggz]])
            self._stagelist.append(["SEx_NEW", STAP.SEx, [lnewimg],
                                    { 'do_seeing': True, 'do_apcor': True }])
            round_stages.append(
                ["round_NEW", STAP.roundint, [lnewimg], { 'gzip': True }])
            zpkw = { 'calstarsfname': lcalstars, 'imgfname': lnewimg }
        else:
            zpkw = { 'calstarsfname': lcalstars }

        # Rezeropointing stage
        self._stagelist.append(["zp_NEW", STAP.zeropoint, [lnewstars], zpkw])

        # And finally, the cleanup
        self._stagelist.append(["fixhdrs", self.fixhdrs,
                               [lnewstars, lsubstars, lsubcands, lsubimg]])
        self._stagelist += round_stages

    def _cleanup_workflow(self):
        """If all goes well, we should have nothing to clean up"""
        pass


# ----------------------------------------------------------------------------
#                        This is a callable script
# ----------------------------------------------------------------------------


if __name__ == "__main__":
    main()
