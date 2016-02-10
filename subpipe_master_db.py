#!/usr/bin/env python

# ============================================================================
# RS 2012/02/27:  SkyMapper Subtraction Pipeline v1
# ----------------------------------------------------------------------------
# My rewrite of the master pipeline script, this time designed properly.
# ============================================================================


# Modules to include

import time
import sys
import os
import glob
import ephem
import pyfits
import argparse
import datetime
import subprocess
from math import sqrt
from STAP.STAP_comm import psmemchk
from STAP.STAP_thread import STAP_thread, STAP_thread_runref
from Utils import Constants
from Utils import Pipeline as Pipeline
from Utils.Catalog import SkybotObject, SDSSObject
from Utils.ProfilerStopwatch import ProfilerStopwatch
from mydjango.jobstats.STAP_API import register_job, register_run
from mydjango.followup.STAP_API import register_xlist
import mydjango.followup.models as fu
import mydjango.jobstats.models as js
from django.db.models import Count
from Utils.Photometry import calc_fluxlim_img
# RS 2012/05/25:  this version not integrated w/Django yet
# def register_job(*args):  pass
# def register_run(*args):  pass
# def register_xsient(*args):  pass
from django.db import transaction
import numpy as np
import re
# ----------------------------------------------------------------------------
#                            The main routine
# ----------------------------------------------------------------------------


def main():
    """Top level ("master process") of the SkyMapper subtraction pipeline"""

    # Get command-line arguments
    today_utc = datetime.datetime.utcnow().strftime("%Y-%m-%d")
    parser = argparse.ArgumentParser(description=
            "Correct headers of incoming images based on scheduler logs")
    parser.add_argument('--utdate', default=None,
            help="path to images (default: %(default)ss)'")
    parser.add_argument('--searchdb', action='store_true', default=False,
            help="search from jobstats database?")
    parser.add_argument('--nproc', type=int, default=32,
            help="number of processors to use (default: 32)")
    parser.add_argument('--waitfornew', action='store_true', default=False,
            help="wait for new files from given UT date?")
    parser.add_argument('--testing', action='store_true', default=False,
            help="test run, don't register or submit any jobs?")
    parser.add_argument('--refonly', action='store_true', default=False,
            help="builing refs only")
    parser.add_argument('--noref', action='store_true', default=False,
            help="do not submit ref jobs")
    args = parser.parse_args()
    print "Using pyfits version", pyfits.__version__

    # Initialize a bunch of important directories that may not exist
    etc = Constants.PipelinePath.etc
    scratch = Constants.PipelinePath.scratch
    scratchetc = Constants.PipelinePath.scratchetc
    subprocess.call("mkdir -p {0}".format(scratch), shell=True)
    subprocess.call("rm -rf {0}/sub*".format(scratch), shell=True)
    subprocess.call("rm -rf {0}/ref*".format(scratch), shell=True)
    subprocess.call("rm -rf {0}/Sky*".format(scratch), shell=True)
    subprocess.call("rm -rf {0}".format(scratchetc), shell=True)
    subprocess.call("cp -r {0} {1}".format(etc,scratchetc), shell=True)

    # Pause lock file name
    pauselockfname = Constants.PipelinePath.pauselockfname

    # In particular, sort out glob strings for *local* dates for today UT.
    # Maybe someday I'll use ssh to query for header strings on markle,
    # which would be a more right answer than grepping on the filename...
    if not args.searchdb:
        if args.utdate is None:
            today = ephem.now()
        else:
            today = ephem.Date(args.utdate.replace("-","/"))
        tomorrow = ephem.Date(today + 1.0)
    # today, tomorrow = ephem.Date(ephem.now() - 1.0), ephem.now()
        today_str = today.datetime().strftime("%Y-%m-%d")
        ngstr = [today.datetime().strftime("%Y-%m-%dT1"),
                 today.datetime().strftime("%Y-%m-%dT2"),
                 tomorrow.datetime().strftime("%Y-%m-%dT0")]
        print "ngstr =", ngstr

    # Set up the queue for pipeline processes
    Q = Pipeline.LowMemQueue(name="Main Pipeline Queue", processors=args.nproc)
    pwatch = ProfilerStopwatch()

    start_run_time = time.time()        # start time for run, use as run ID
    os.chdir(scratch)                   # scratch should be our working path
    if not args.testing:
        register_run(start_run_time)    # register run in Django

    donefnames = [ ]                    # NEW filenames already subtracted
    pointingcache = [ ]                 # Tags of pointings w/cached info
    fieldccdlock = { }                  # Field+CCD combinations in progress
    loop_delay = 30                     # Minimum duration of loop in seconds
    pause_delay = 600                   # Length of pause cycle in seconds

    # Fetch a list of everything worth following before the run started
    # !!! Need to turn off followup for old objects!!!
    xprelist = fu.Transient.objects.select_related().filter(
        type__type__in=Constants.follow_types) #isactive=1

    # Get unique identifier for this run
    runtag = Constants.get_runtag(start_run_time)

    def locktag(thing):
        return "{0}_{1:02d}".format(thing.field, thing.ccd)

    # The main event loop
    idle = False
    while 1:

        start_loop_time = time.time()
        # Check for new image data arriving from markle
        # There should be one NEW image in this directory for each exposure
        # of each CCD over the course of a night.
        # This assumes there's another process to download mosaic images
        # from markle, split and overscan-subtract them.
        if args.searchdb:
            # FY - 2015/07/19 find freshfnames from database search
            # find all pointing with status < 10 (not bad), and haven't been processed (no pipelinejob)
            freshfnames = []
            newpts=js.SkymapperImage.objects.filter(#pointing__program__in=['3PS'],
                status=0,#pipelinejob__isnull=True,
                pointing__dateobs__gte=datetime.datetime(2014,7,1,0,0,0),
                ).exclude(pointing__field__id=0).order_by('pointing__dateobs').select_related('pointing','pointing__field')[:50]
            #.order_by('-pointing__dateobs')
            lennewpts=len(newpts)
            print "Found new images:",lennewpts
            with transaction.commit_on_success():
                for newpt in newpts:
                    freshfname =Constants.FilenamesNew.absolute_dir(newpt.pointing.field.id,newpt.pointing.filter,newpt.ccd)+'/'+newpt.pointing.filename+'_%d.fits'%newpt.ccd
                    if os.path.exists(freshfname):
                        if freshfname not in donefnames:
                            freshfnames.append(str(freshfname))
                    else:
                        newpt.status=99  #non-exist
                        newpt.save()
            if len(freshfnames)==0:
                print "No images on disk?"
                print freshfnames
        else:
            CPP_base_glob = Constants.PipelinePath.new + "/*/*/*/*"
            freshfnames = [f for f in
                           # glob.glob(Constants.PipelinePath.new_glob_str)
                           glob.glob(CPP_base_glob + ngstr[0] + "*.fits") +
                           glob.glob(CPP_base_glob + ngstr[1] + "*.fits") +
                           glob.glob(CPP_base_glob + ngstr[2] + "*.fits")
                           if f not in donefnames]
        #"""
        freshfnames = [f for f in freshfnames if f not in donefnames]
        #"""        

        # Make sure they've finished copying!
        # RS:  Do them all at once; since each call to ensurecopy() means
        # at least a 5-second delay, sequential calls are much slower!
        #freshfnames = ensurecopy(freshfnames)
        memtot, pmemtot = psmemchk('python subpipe_master_db.py')
        if not idle:
            print "subpipe_master.py:  running main loop at local time ",
            print time.asctime(time.localtime(start_loop_time))
            print "Running {0:d}/{1:d} SUBs now; {2:d} NEW images waiting to run;"\
                    .format(len(fieldccdlock),Q.N_processors,len(freshfnames)),
            print "memory use {0:.1f} MB ({1:.1f}%)".format(memtot, pmemtot)
            sys.stdout.flush()
        idle = (len(fieldccdlock) == 0 and len(freshfnames) == 0)
#        if idle and (not args.waitfornew or
#                     ephem.now().datetime().strftime("%Y-%m-%d") != today_str):
        if idle and (not args.waitfornew):
            print "subpipe_master.py:  ran out of stuff to do"
            print "Exiting normally at", time.asctime(time.localtime(start_loop_time))
            return

        # Start as many new jobs as we have empty slots
        for newfname in freshfnames:
            # Check to see if we've paused.  If we have, break out of this
            # loop but continue to retrieve results for finished jobs.
            if os.path.exists(pauselockfname):
                print "subpipe_master.py:  Someone hit the pause switch"
                print "Not starting any new jobs, but collecting old ones."
                sys.stdout.flush()
                break

            # Pick up a new NEW
            new = Constants.FilenamesNew(newfname, runtag=runtag)
            
            # If there's no room left, don't submit any more
            if len(fieldccdlock) >= Q.N_processors:  break

            # RS 2015/06/24:  Field 0000 is actually code for a failed WCS
            # solution, so don't even bother submitting these.
            if int(new.field) == 0:
                print "Skipping", new.basefname,
                print "because there's a problem with the WCS"
                donefnames.append(new.absfname)
                continue

            # RS 2012/02/23:  Don't process the same field+CCD combination
            # at the same time, it'll screw up candidate cross-referencing.
            fieldccdtag = locktag(new)
            if fieldccdtag in fieldccdlock:  continue

            # RS 2013/10/30:  Basic quality control
            maxelong=1.3
            if False: #new.seeing == None or new.elong > maxelong: # or new.seeing < 10.0
                print "Skipping", new.basefname,
                donefnames.append(new.absfname)
                im=js.SkymapperImage.objects.get(pointing__filename='_'.join(new.basefname.split('_')[:-1]),ccd=new.subfield)
                if im.status==0:
                    im.status=80  #poor quality
                    im.save()
                if new.seeing == None:
                    print "because seeing not known"
                # RS 2015/04/28:  Removing this restriction for now because
                # in bad seeing mode we're liable to get a lot of terrible
                # images.  Keeping elongation cut though.
                # elif new.seeing > 10.0:
                #     print "because seeing is terrible",
                #     print "({0} > 10.0)".format(new.seeing)
                elif new.elong > maxelong:
                    print "because elongation is terrible",
                    print "({0} > {1})".format(new.elong,maxelong)
                continue

            # RS 2012/02/27:  Find a REF image for this NEW.  For now, just
            # pick the first one in this filter and area of the sky.
            # RS 2014/07/01:  Check to make sure that the list of pre-existing
            # transients actually come from the field being considered!
            refname = None
            refonly=args.refonly
            if new.program =='MS': refonly=True
            if not refonly:
                xprelist_loc = xprelist.filter(field=new.field)
                refname = pick_ref(new, xprelist=xprelist_loc, verbose=True)
                sys.stdout.flush()
            if refname is None and not refonly and args.noref:
                print "Skipping", new.basefname
                donefnames.append(new.absfname)
            if refname == None or refonly:
                # RS 2013/09/12:  Added feature which runs an image as a REF
                # if no suitable REF exists for this part of the sky.
                #print "Can't find reference image for " + new.basefname
                print "Submitting it as reference-building job instead."
                sub = Constants.FilenamesRef(newfname, runtag=runtag)
                workflow = STAP_thread_runref(
                        name=sub.id, workdir=sub.workdir,
                        logfile=sub.log_file_name)
                Qargs = (workflow, new.absfname)
            else:
                # FY 2015/08: before subtraction job, update xref
                # Ideally this should be updated when database is updated from web, but web server doesn't have permission for the xref files.
                xreffile = "{0}/{1:04d}/{2:02d}/{1:04d}_{2}_xref.fits".format(Constants.PipelinePath.xref,new.field,new.subfield)
                if os.path.exists(xreffile):
                    trans=list(fu.Transient.objects.filter(field=new.field,ccd=new.subfield).select_related('type__type'))
                    hdulist=pyfits.open(xreffile,mode='update')
                    for (key, value) in hdulist[0].header.items():
                        namematch = re.match("NAME(\d{4})", key)
                        if namematch:
                            id=int(namematch.group(1))
                            transtype=[tran.type.type for tran in trans if tran.name==value]
                            if len(transtype)>0:
                                hdulist[0].header.update("TYPE{0:04d}".format(id),
                                                         str(transtype[0]),
                                                         "Type of source with index {0}".format(id))
                    hdulist.close()

                # RS 2012/02/23:  Subs are now identified by start_run_time,
                # (i.e., by the time the script started), not start_loop_time.
                sub = Constants.FilenamesSub(newfname, runtag)
                workflow = STAP_thread(name=sub.id, workdir=sub.workdir,
                                       logfile=sub.log_file_name)
                ref = Constants.FilenamesRef(refname)
                Qargs = (workflow, new, ref, runtag)

            # RS 2012/08/29:  If this is the first time we've seen this
            # pointing, run a once-per-pointing job to cache relevant info.
            myptag = Constants.datetime2str(new.dateobs)
            if myptag not in pointingcache:
                cache_pointing_info(new)
                pointingcache.append(myptag)

            # RS 2013/09/13:  HACK -- If we already processed this workflow
            # in a previous run, skip it.  Also, find a more elegant way to
            # do this later.
            workflow._setup_workflow(*(Qargs[1:]))
            if refname==None or refonly:
                done=all([os.path.exists(fn) for fn in workflow._copy_back if
                          'wcs.fits.stars' in fn])
                if done:
                    print "Skipping {0}, we already ran it".format(workflow.name)
                    donefnames.append(new.absfname)
                    im=js.SkymapperImage.objects.get(pointing__filename='_'.join(new.basefname.split('_')[:-1]),ccd=new.subfield)
                    if im.status==0:
                        im.status=1
                        im.save()
                    continue
            else:
                done=False #all([os.path.exists(fn) for fn in workflow._copy_back if
                          #'wcs.fits.stars' in fn or 'remap.fits.stars' in fn])
                #done=all([os.path.exists(fn) for fn in workflow._copy_back if
                #          'wcs.fits.stars' in fn]) #not redo subtraction for now
                if done:
                    print "Skipping {0}, we already ran it".format(workflow.name)
                    donefnames.append(new.absfname)
                    im=js.SkymapperImage.objects.get(pointing__filename='_'.join(new.basefname.split('_')[:-1]),ccd=new.subfield)
                    if im.status<2:
                        im.status=2
                        im.save()                                
                    continue

            # If we're just testing, don't submit any jobs
            if args.testing:
                print "Not submitting job {0} to queue, "\
                      "just testing setup".format(workflow.name)
                donefnames.append(new.absfname)
                continue

            # Try to submit the job to the queue.
            status = Q.start_job(*Qargs)
            if status == 0:
                print "Submitted job {0} to queue".format(workflow.name)
                fieldccdlock[fieldccdtag] = sub
                donefnames.append(newfname)
            elif status == 1:
                print "Couldn't submit job ", workflow.name,
                print "; will try again later"
            elif status == -1:
                print "Fatal error setting up for job", workflow.name
                donefnames.append(newfname)
            sys.stdout.flush()

        # Wait for a bit
        time.sleep(loop_delay/2)

        # Deal with output of terminated processes
        while not Q.empty():
            # Check the overhead associated with writing to disk
            pwatch.start("Queue emptying cycle")
            wf = pwatch.profile_call(Q.finish_job)
            job_out = wf.output
            # print "Retrieved output from job", job_out["name"]
            if hasattr(wf, 'refnames'):
                refnames = wf.refnames
            else:
                refnames = None

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

            # Register any transients that might have been found.
            xrefout = [o["results"] for o in job_out["stages"]
                       if o["name"] == "xref" and "results" in o]
            if len(xrefout) == 1:
                register_xlist(xrefout[0])

            # Register the job along with its result, then release the lock
            # on this field+CCD combination.
            for fieldccdtag in fieldccdlock:
                sub = fieldccdlock[fieldccdtag]
                if sub.id == job_out["name"]:
                    pwatch.profile_call(
                            register_job, sub, job_out, refmeta=refnames)
                    del fieldccdlock[fieldccdtag]
                    break
            pwatch.stop("Queue emptying cycle")
            print "Profile for", job_out["name"], "--",
            pwatch.report()
            pwatch.reset_all()
            sys.stdout.flush()

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
#                        Various helpful functions
# ----------------------------------------------------------------------------


def pick_ref(new, xprelist=None, verbose=False):
    """Pick a good REF image for this NEW"""
    mypath = Constants.PipelinePath.ref + "/{0}/*_wcs.fits".format(
            new.relative_dir)
    # Pick up both normal FITS files and compressed FITS files
    files = glob.glob(mypath)
    files += glob.glob(mypath + ".gz")
    if verbose:
        print "Found", len(files), "possible REFs in", mypath
    if len(files) == 0:
        print "ERROR: no references found at ALL for " + new.basefname
        return None
    else:
        # Don't settle for a substandard reference.
        refname = None
        # All processed REF images should now have seeing and zeropoint info
        # in their headers.  We *should* know the NEW seeing now as well.
        # Pick a reasonable range of seeing which should result in a good
        # subtraction.  Within this range, pick the deepest image.
        best_dt, best_seeing, best_maglim = -1.0, 99.9, 0.0
        if new.seeing is None:
            seelo, seehi = 0.0, 99.9
        elif new.seeing < 6.01:
            seelo, seehi = 0.0, new.seeing
        else:
            # seelo, seehi = sqrt(new.seeing**2 - 9.0), new.seeing
            seelo, seehi = sqrt(new.seeing**2 - 36.0), new.seeing
        if verbose:  print "seelo, seehi =", seelo, seehi
        # If there is a thing worthy of follow-up discovered on this field
        # and CCD, we need to make sure the ref is no *younger* than the
        # oldest detection of any such thing.
        if xprelist and len(xprelist) > 0:
            latest_allowed = ephem.Date(ephem.Date(min(
                [min([p.dateobs for p in x.pointing_det.get_query_set()])
                 for x in xprelist])) - 14.0)
        else:
            latest_allowed = ephem.Date(ephem.Date(new.dateobs) - 14.0)
        # Now we go.
        for fn in files:
            if verbose:  print "{0}: ".format(fn),
            try:
                hdr = pyfits.getheader(fn)
            except IOError:
                print "-- rejected (trapped IOError / file unreadable)"
                continue
            # Status debugging statements
            kwhash = { 'DATE-OBS': 'None', 'SEEING': 'None', 'ZPMAG': 'None', 'BKGSIG': 'None' }
            # Don't pick an image as REF with no DATE-OBS, SEEING or ZPMAG
            complete = True
            for kw in kwhash:
                if kw in hdr:
                    kwhash[kw] = hdr[kw]
                else:
                    if verbose:
                        print "-- rejected ({0} not in hdr)".format(kw)
                    complete = False
            if not complete:  continue
            # Don't subtract an image from itself, and don't subtract an image
            # from another image taken at a later time (light curve ordering)
            mydateobs = ephem.Date(Constants.str2datetime(hdr['DATE-OBS']))
            if mydateobs > latest_allowed:
                if verbose:  print "-- rejected (taken after date {0})".format(
                        latest_allowed)
                continue
            dt = latest_allowed - mydateobs
            # If we don't know what the NEW's seeing is, pick as the REF
            # the best seeing from among the available REFs
            seeing, zp, bkgsig = hdr['SEEING'], hdr['ZPMAG'], hdr['BKGSIG']
            maglim = -2.5*np.log10(calc_fluxlim_img(bkgsig, CL=0.50)) + zp
            if new.seeing == None and seeing < best_seeing and maglim > best_maglim:
                refname, best_dt, best_seeing, best_maglim = fn, dt, seeing, maglim
                if verbose:  print "-- best so far (unknown seeing)"
            # If we know what the NEW's seeing is, pick as the REF
            # the deepest image with seeing better than the NEW
            elif seelo < seeing < seehi and maglim > best_maglim:
                refname, best_dt, best_seeing, best_maglim = fn, dt, seeing, maglim
                if verbose:  print "-- best so far (deepest w/suitable seeing)"
            elif seeing > seehi:
                if verbose:  print "-- rejected (REF seeing worse than NEW)"
                continue
            # Otherwise, just move on.
            elif seeing <seelo:
                if verbose: print "-- rejected (REF seeing worse than NEW)"
                continue
            else:
                if verbose:  print "-- rejected (worse than current best)"
        # print "Chose {0} with seeing {1}".format(refname, best_seeing)
        if refname == None:
            if verbose:
                print "No suitable reference found for " + new.basefname
            return None
        else:
            if verbose:
                print "Picked reference with dt, seeing, maglim =",\
                       best_dt, best_seeing, best_maglim
            return refname


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
            fsize0[f] = 1
            fsize1[f] = 1
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
                fsize1[f] = 1
                badfiles.append(f)

    # Return all not-bad files that are done copying
    goodfiles = [f for f in files if fsize0[f] == fsize1[f]
                                  and f not in badfiles]
    return goodfiles


def cache_pointing_info(meta):
    """Cache important info for this pointing

    Run once for each pointing synchronously, before subtractions are run.
    This allows us to cut down on I/O to external sources for selected
    applications.  For example, if we're going to need catalog information
    from some catalog on the Web, we can cache it all in this stage.
    """

    smfield = Constants.SMfield_lookup(id=meta.field)
    dtobs = Constants.datetime2str(meta.dateobs)
    print "Running pre-cache for pointing {0} {1}...".format(smfield, dtobs)

    # Query SkyBot
    if not os.path.exists(meta.skybot_cachefname):
        start_cache_time = time.time()
        subprocess.call("mkdir -p {0}".format(
                os.path.dirname(meta.skybot_cachefname)), shell=True)
        roids = SkybotObject.webquery(ra=smfield.ra, dec=smfield.dec,
                                      datetime=dtobs, rm=120,
                                      cachefname=meta.skybot_cachefname)
        end_cache_time = time.time()
        print "Time for SkyBot to run:  {0:.1f} sec".format(
            end_cache_time - start_cache_time)

    # Query SDSS
    if not os.path.exists(meta.sdss_cachefname):
        # Only do this if file doesn't exist
        subprocess.call("mkdir -p {0}".format(
                os.path.dirname(meta.sdss_cachefname)), shell=True)
    # FY - if ra is less than -30, don't bother to check, but write an empty file
        if smfield.dec<-30: os.system('touch {0}'.format(meta.sdss_cachefname))
        else:
            sdss = SDSSObject.webquery(ra=smfield.ra, dec=smfield.dec,
                                       rm=120,
                                       cachefname=meta.sdss_cachefname)

# ----------------------------------------------------------------------------
#                        This is a callable script
# ----------------------------------------------------------------------------


if __name__ == "__main__":
    main()
