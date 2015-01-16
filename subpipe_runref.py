#!/usr/bin/env python

# ============================================================================
# RS 2012/02/27:  SkyMapper Subtraction Pipeline v1
# ----------------------------------------------------------------------------
# My rewrite of the master pipeline script, this time designed properly.
# This version is for processing new images for the REF cache.
# ============================================================================


# Modules to include

import time
import sys
import os
import glob
import pyfits
import subprocess
from STAP.STAP_comm import psmemchk
from Utils import Constants
from Utils import Pipeline as Pipeline
# from mydjango.jobstats.STAP_API import register_job, register_run

# RS 2012/05/25:  this version not integrated w/Django yet
def register_job(*args):  pass
def register_run(*args):  pass


# ----------------------------------------------------------------------------
#                            The main routine
# ----------------------------------------------------------------------------


def main_runref():
    """Master process of a REF-cache-building sub-pipeline"""

    # Initialize a bunch of important directories that may not exist
    from STAP.STAP_thread import STAP_thread_runref
    etc = Constants.PipelinePath.etc
    scratch = Constants.PipelinePath.scratch
    scratchetc = Constants.PipelinePath.scratchetc
    subprocess.call("mkdir -p {0}".format(scratch), shell=True)
    subprocess.call("rm -rf {0}/sub*".format(scratch), shell=True)
    subprocess.call("rm -rf {0}/Sky*".format(scratch), shell=True)
    subprocess.call("rm -rf {0}".format(scratchetc), shell=True)
    subprocess.call("cp -r {0} {1}".format(etc,scratchetc), shell=True)

    # Set up the queue for pipeline processes
    Q = Pipeline.LowMemQueue(name="Main Pipeline Queue", processors=32)

    start_run_time = time.time()        # start time for run, use as run ID
    os.chdir(scratch)                   # scratch should be our working path

    donefnames = [ ]                    # NEW filenames already subtracted
    fieldccdlock = { }                  # Field+CCD combinations in progress
    loop_delay = 30                     # Minimum duration of loop in seconds

    def locktag(thing):
        return "{0:04d}|{1:02d}".format(thing.field,thing.ccd)

    freshfnames = [f for f in
                   glob.glob(Constants.PipelinePath.new_glob_str)
                   if f not in donefnames] # [:16]

    # The main event loop
    while 1:

        start_loop_time = time.time()
        print "subpipe_runref.py:  running main loop at local time ",
        print time.asctime(time.localtime(start_loop_time))

        # Check for new image data arriving from markle
        # There should be one NEW image in this directory for each exposure
        # of each CCD over the course of a night.
        # This assumes there's another process to download mosaic images
        # from markle, split and overscan-subtract them.
        freshfnames = [f for f in
                       glob.glob(Constants.PipelinePath.new_glob_str)
                       if f not in donefnames]
        """
        freshfnames = [f for f in freshfnames if f not in donefnames]
        """
        # Make sure they've finished copying!
        # RS:  Do them all at once; since each call to ensurecopy() means
        # at least a 5-second delay, sequential calls are much slower!
        freshfnames = ensurecopy(freshfnames)
        memtot, pmemtot = psmemchk('python subpipe_runref.py')
        print "Running {0:d}/{1:d} REFs now; {2:d} NEW images waiting to run;"\
               .format(len(fieldccdlock),Q.N_processors,len(freshfnames)),
        print "memory use {0:.1f} MB ({1:.1f}%)".format(memtot, pmemtot)
        sys.stdout.flush()

        # Start as many new jobs as we have empty slots
        for newfname in freshfnames:
            new = Constants.FilenamesNew(newfname)

            # If there's no room left, don't submit any more
            if len(fieldccdlock) >= Q.N_processors:  break

            # RS 2012/02/23:  Don't process the same field+CCD combination
            # at the same time, it'll screw up candidate cross-referencing.
            fieldccdtag = locktag(new)
            if fieldccdtag in fieldccdlock:  continue

            # RS 2012/04/17:  If we're just processing references,
            # we don't need to find a pre-existing REF.
            jobid = new.basefname.replace(".fits","")
            workflow = STAP_thread_runref(name=jobid,
                    workdir=scratch+"/"+jobid, logfile=jobid+".log")
            Qargs = (workflow, new.absfname)

            # RS 2013/09/05:  HACK -- If we already processed this reference
            # in a previous run, skip it.  Also, find a more elegant way to
            # do this later.
            workflow._setup_workflow(new.absfname)
            if all([os.path.exists(fn) for fn in workflow._copy_back]):
                print "Skipping {0}, we already ran it".format(workflow.name)
                donefnames.append(new.absfname)
                continue

            # Try to submit the job to the queue.
            status = Q.start_job(*Qargs)
            if status == 0:
                print "Submitted job {0} to queue".format(workflow.name)
                fieldccdlock[fieldccdtag] = jobid
                donefnames.append(new.absfname)
            elif status == 1:
                print "Couldn't submit job ", workflow.name,
                print "; will try again later"
            elif status == -1:
                print "Fatal error setting up for job", workflow.name
                donefnames.append(new.absfname)
            sys.stdout.flush()

        # Wait for jobs to finish, if it makes any sense to do so
        mid_loop_time = time.time()
        loop_wait = min(loop_delay/2, loop_delay - (mid_loop_time - start_loop_time))
        if loop_wait > 0:  time.sleep(loop_wait)

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

            # Release the lock on this field+CCD combination.
            for fieldccdtag in fieldccdlock:
                if fieldccdlock[fieldccdtag] == job_out["name"]:
                    del fieldccdlock[fieldccdtag]
                    break

        # If necessary, hang out till the end of this delay cycle
        sys.stdout.flush()
        if len(fieldccdlock) < Q.N_processors and len(freshfnames) > 0:
            continue
        end_loop_time = time.time()
        loop_wait = loop_delay-(end_loop_time-start_loop_time)
        if loop_wait > 0:  time.sleep(loop_wait)


# ----------------------------------------------------------------------------
#                        Various helpful functions
# ----------------------------------------------------------------------------


def ensurecopy(files, delay=5, timeout=60, path=None):
    """Ensure that files copied from remote systems have finished copying
    
    For file in files, check to make sure its size is not changing with time.
    Check every delay seconds, up to timeout, then return a list of files
    which have finished copying.  Excludes files which have disappeared while
    we were waiting (we probably didn't want them anyways).
    """

    t = 0
    (badfiles,fsize0,fsize1) = ([ ],{ },{ })
    if path != None:
        files = ["{0}/{1}".format(os.path.normpath(path),
                                  os.path.basename(f)) for f in files]

    # Get initial sizes
    for f in files:
        try:
            fsize0[f] = 1
            fsize1[f] = os.path.getsize(f)
        except OSError as e:
            badfiles.append(f)

    # If we have *only* bad files, return
    if files == badfiles:
        return []

    # Actually run the loop now
    while fsize0 != fsize1 and t < timeout:
        time.sleep(delay)
        t = t + delay
        fsize0 = fsize1
        for f in files:
            try:
                fsize1[f] = os.path.getsize(f)
            except OSError:
                badfiles.append(f)

    # Return all not-bad files that are done copying
    goodfiles = [f for f in files if fsize0[f] == fsize1[f]
                                  and f not in badfiles]
    return goodfiles


# ----------------------------------------------------------------------------
#                        This is a callable script
# ----------------------------------------------------------------------------


if __name__ == "__main__":
    main_runref()
