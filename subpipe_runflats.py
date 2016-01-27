#!/usr/bin/env python

# ============================================================================
# RS 2012/02/27:  SkyMapper Subtraction Pipeline v1
# ----------------------------------------------------------------------------
# My rewrite of the master pipeline script, this time designed properly.
# This version is for processing new twilight flat images into superflats
# and bad pixel masks (which may be used in the pipeline soon!).
# ============================================================================


# Modules to include

import time
import sys
import os
import os.path
import glob
import pyfits
import subprocess
from Utils import Constants
from Utils import Pipeline as Pipeline
# from mydjango.jobstats.STAP_API import register_job, register_run
from datetime import datetime, timedelta
# RS 2012/05/25:  this version not integrated w/Django yet
def register_job(*args):  pass
def register_run(*args):  pass


# ----------------------------------------------------------------------------
#                            The main routine
# ----------------------------------------------------------------------------


def main_runflats(date):
    """Master process of a flat-field-building sub-pipeline"""
        
    # Initialize a bunch of important directories that may not exist
    from STAP.STAP_thread import STAP_thread_runflats
    etc = Constants.PipelinePath.etc
    scratch = Constants.PipelinePath.scratch
    scratchetc = Constants.PipelinePath.scratchetc
    subprocess.call("mkdir -p {0}".format(scratch), shell=True)
    subprocess.call("rm -rf {0}/flat*".format(scratch), shell=True)
    subprocess.call("rm -rf {0}".format(scratchetc), shell=True)
    subprocess.call("cp -r {0} {1}".format(etc,scratchetc), shell=True)

    # Set up the queue for pipeline processes
    Q = Pipeline.Queue(name="FLATS Queue", processors=12)

    start_run_time = time.time()        # start time for run, use as run ID
    os.chdir(scratch)                   # scratch should be our working path

    runningjobs = { }                   # CAL jobs running now
    loop_delay = 30                     # Minimum duration of loop in seconds

    def jobname(fn):
        return "flat_{0}_{1:02d}".format(fn.filter, fn.ccd)

    # We should already have all the files we need.  Group by filter and CCD.
    imglist = { }
    date_object = datetime.strptime(date, '%Y-%m-%d')+timedelta(1)
    nextdate=date_object.strftime('%Y-%m-%d')
    files=glob.glob(Constants.PipelinePath.cal+"/newflats/*/*/*%sT[12]*.fits"%(date))
    files1=glob.glob(Constants.PipelinePath.cal+"/newflats/*/*/*%sT0*.fits"%(nextdate))
    if files1: 
        for filename in files1:files.append(filename)
    for fn in files:
        header=pyfits.getheader(fn)
        if header['MEDPIX']<5000. or header['MEDPIX']>50000:
            continue
        names = Constants.FilenamesNewflat(fn)
        jobid = jobname(names)
        if jobid not in imglist:  imglist[jobid] = [ ]
        imglist[jobid].append(fn)
        

    badid=[]
    for jobid in imglist:
        print imglist[jobid]
        if len(imglist[jobid])<5:
            print "Less than 5 images found for jobid",jobid
            badid.append(jobid)
    for bid in badid:
        del imglist[bid]


    # The main event loop
    while len(imglist) > 0:

        start_loop_time = time.time()
        print "subpipe_runflats.py:  running main loop at local time ",
        print time.asctime(time.localtime(start_loop_time))
        print "Running {0:d}/{1:d} CALs now; {2:d} CAL groups waiting to run"\
            .format(len(runningjobs),Q.N_processors,len(imglist)-len(runningjobs))
        sys.stdout.flush()

        # Start as many new jobs as we have empty slots
        for jobid in imglist:
            
            # Don't submit this one if it's already running
            if jobid in runningjobs:  continue
            # If there's no room left, don't submit any more
            if len(runningjobs) >= Q.N_processors:  break
            workflow = STAP_thread_runflats(name=jobid,
                       workdir=scratch+"/"+jobid, logfile=jobid+".log")
            Qargs = (workflow, imglist[jobid])
            
            status = Q.start_job(*Qargs)
            if status == 0:
                print "Submitted job {0} to queue".format(workflow.name)
                runningjobs[jobid] = "running"
            elif status == 1:
                print "Couldn't submit job ", jobid
            elif status == -1:
                print "Fatal error setting up for job", jobid

        sys.stdout.flush()
        time.sleep(loop_delay/2)

        # Deal with output of terminated processes
        while not Q.empty():
            wf = Q.finish_job()
            job_out=wf.output
            # print "Retrieved output from job", job_out["name"]
            
            # Figure out what happened to the job
            jobid = job_out["name"]
            print jobid,
            if "exception" in job_out["stages"][-1]:
                print "failed in stage", job_out["stages"][-1]["name"],
                print "from an uncaught exception:"
                print job_out["stages"][-1]["exception"]
            elif "results" in job_out["stages"][-1]:
                print "has successfully completed"
            else:
                print "probably got killed before finishing"

            # Remove this group of images from the hash of "pending" lists,
            # and from the list of jobs currently running.
            del imglist[jobid]
            del runningjobs[jobid]

        # If necessary, hang out till the end of this delay cycle
        sys.stdout.flush()
        if len(runningjobs) < Q.N_processors and len(imglist)>0:  continue
        if len(imglist)>0:
            end_loop_time = time.time()
            loop_wait = loop_delay-(end_loop_time-start_loop_time)
            if loop_wait > 0:  
                print 'waiting for next loop after',loop_wait
                time.sleep(loop_wait)

    return
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
    main_runflats(sys.argv[1])
