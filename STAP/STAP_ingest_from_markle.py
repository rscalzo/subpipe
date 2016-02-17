#! /usr/bin/env python

# ============================================================================
# RS 2012/02/08:  "Ingest" new images from markle
# ----------------------------------------------------------------------------
# This is a clean-up of the previous script of its name to make it simpler,
# to improve security, and to remove the need to have a separate process
# running on markle (or whatever the remote machine is called).
# ============================================================================

import os
import re
import sys
import time
import math
import glob
import ephem
import pyfits
import argparse
import datetime
import subprocess
import threading
import traceback
from STAP_comm import STAP_callexternal, psmemchk
from Utils import Pipeline as Pipeline
from Utils.Constants import PipelinePath as CPP, get_runtag, Imager, decode_obs_id
from Utils.Constants import SMfield_lookup, Filenames, FilenamesNew, FilenamesNewflat
from Utils.TrackableException import TrackableException as IngestFail
# from mydjango.jobstats.STAP_API import register_job, register_run
def register_job(*args):  pass
def register_run(*args):  pass
from mydjango.jobstats.STAP_API import register_pointing,register_image
import mydjango.jobstats.models as js
#only when miner is down
#def register_pointing(*args):  pass

#define home, raw, new and cal path
# RS:  We need these to function down there.
scratch = os.path.normpath(CPP.scratch + "/../subpipe_ingest")

# set minimum available disk space, size of five images
bytes_to_mb = 1.0/(1024**2)
img_size_in_bytes = 585881280
min_space_avail = 5*img_size_in_bytes*bytes_to_mb

#set loop delay
loop_delay=90

class Logger(object):
    def __init__(self, filename="Default.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "a")      
    def write(self, message):                       
        self.terminal.write(message)                                
        self.log.write(message)                                                 
    def flush(self):
        self.terminal.flush()
        self.log.flush()

# Enable logging to output without redirecting stdout from within the shell.
# This is crucial to allow a fetch to run in a window without typing in my
# password all the friggin' time.
sys.stdout = Logger("ingest_many_nights.log")


# ============================================================================
#                              The main event
# ============================================================================


def main():

    # Get command-line arguments
    today_utc = datetime.datetime.utcnow().strftime("%Y-%m-%d")
    default_remote_base = "skymap@markle:/export/markle2/images/"
    default_remote = default_remote_base + today_utc
    parser = argparse.ArgumentParser(description=
            "Correct headers of incoming images based on scheduler logs")
    parser.add_argument('--fetchfrom', default=default_remote,
            help="path to images (default: %(default)ss)'")
    parser.add_argument('--local', action='store_true', default=False,
            help="shortcut:  set --fetchfrom to CPP.raw")
    parser.add_argument('--nproc', type=int, default=8,
            help="number of processors to use (default: 8)")
    parser.add_argument('--keep_workdir', action='store_true', default=False,
            help="keep the contents of /ramdisk/workdir for debugging?")
    parser.add_argument('--production', action='store_true', default=False,
            help="run in production mode, importing from current UT date?")
    parser.add_argument('--obstype', default='MS_*|TFF_*|3PS_42|BAD_42',
            help="only retrieve images with specified observation type and survey id") 
    parser.add_argument('--imageid', default='53|89|805|053|089|169',
            help="only retrieve images with specified imageid") 
    parser.add_argument('--checkdb', action='store_true', default=False,
                        help="check database for already ingested files")
    args = parser.parse_args()
    print "Using pyfits version", pyfits.__version__

    # Paths where stuff is
    if args.local:
        acct, remotepath, today_utc = None, CPP.raw, ""
    elif args.production:
        acct, remotepath = re.search(
                "(\w+\@[^:]+)*:?(.*)", default_remote_base).groups()
        remotepath += today_utc
    elif '@' in args.fetchfrom:
        acct, remotepath = re.search(
                "(\w+\@[^:]+)*:?(.*)", args.fetchfrom).groups()
        match = re.search("(\d{4}-\d{2}-\d{2})", remotepath)
        print "Fetch from: acct {0}, remotepath {1}".format(acct, remotepath)
        if match:  today_utc = match.group(1)
    else: 
        acct, remotepath, today_utc = None, args.fetchfrom, ""
        args.local = True
    # Use half the processors for downloading and half for splitting
    print "Using", args.nproc, "processors for downloading and splitting"
    nproc = args.nproc/2
    if nproc < 1:  nproc = 1
    if nproc > 8:  nproc = 8

    # Set up the queues for pipeline processes
    Qdown = Pipeline.LowMemQueue(name="Download Queue",
            processors=nproc, keep_workdir=False)
    Qsplit = Pipeline.LowMemQueue(name="SplitAmp Queue",
            processors=nproc, keep_workdir=args.keep_workdir)

    donefnames = [ ]    # Files which have already been processed
    manifest = { }      # Manifest of image properties from scheduler logs

    start_run_time = time.time()        # start time for run, use as run ID
    register_run(start_run_time)
    runtag = get_runtag(start_run_time)

    # Create a manifest of images we expect to find
    print "STAP_ingest.py:  initializing manifest at local time ",
    print time.asctime(time.localtime(time.time()))
    if args.local:
        print "remotepath =", remotepath
        myfiles = glob.glob(remotepath + "/Skymapper*.fits")
        matchlist1 = [re.search("(\d{4}-\d{2}-\d{2})T[12]", fn) for fn in myfiles]
        matchlist2 = [re.search("(\d{4}-\d{2}-\d{2})T0", fn) for fn in myfiles]
        datelist1 = set([m.group(1) for m in matchlist1 if m])
        datelist2 = set([(datetime.datetime.strptime(m.group(1),'%Y-%m-%d')-datetime.timedelta(days=1)).strftime('%Y-%m-%d') for m in matchlist2 if m])
        datelist=set(list(datelist1)+list(datelist2))
        print "myfiles = ", myfiles
        print "Fetching logs from", datelist
        schedule_log_list = [ ]
        for date in datelist:
            status, log_sublist = fetch_schedule_logs(date)
            schedule_log_list += log_sublist
    elif not args.production:
        status, schedule_log_list = fetch_schedule_logs(today_utc)
    else:
        status, schedule_log_list = 0, [ ]
    print "schedule_log_list =", schedule_log_list
    # if status !=0: status, schedule_log_list= 0,[]
    if status != 0:
        print "STAP_ingest.py:  Couldn't fetch schedule logs!"
        sys.stdout.flush()
        #return -1
    else:
        for logfile in schedule_log_list:
            if os.path.exists(logfile) and re.search(today_utc, logfile):
                print "Updating manifest from", logfile
            # manifest.update(manifest_from_scheduler_log(logfile))
                manifest = manifest_from_scheduler_log(logfile, manifest)
                sys.stdout.flush()

    remotefnames = [ None ]
    idle = False

    while 1:    # "for-EV-errrrr"

        start_loop_time = time.time()

        # RS 2013/09/19:  If running in "production" mode, check whether
        # UTC has changed, and add to manifest if so.  There may not yet be
        # any schedule logs 
        if args.production:
            acct, remotepath = re.search(
                    "(\w+\@[^:]+)*:?(.*)", default_remote_base).groups()
            today_utc = datetime.datetime.utcnow().strftime("%Y-%m-%d")
            remotepath += today_utc
            status, schedule_log_list = fetch_schedule_logs(today_utc)
            if len(schedule_log_list) < 1 or status != 0:
                print "STAP_ingest.py:  Couldn't download schedule logs"
                print "Waiting 30 minutes and hoping problem is resolved."
                sys.stdout.flush()
                time.sleep(1800)
                continue
            for logfile in schedule_log_list:
                if os.path.exists(logfile) and re.search(today_utc, logfile):
                    print "Updating manifest from", logfile
                    manifest = manifest_from_scheduler_log(logfile, manifest)
                    sys.stdout.flush()

        # RS:  See which new files have appeared on the "remote" system.
        # If this fails for some reason, print an error message and wait
        # for a while for the user to notice and/or fix something.
        status, remotefnames = ls_remote(acct, remotepath, obstype=args.obstype, imageid=args.imageid)
        if status != 0:
            print "STAP_ingest.py:  remote ls failed"
            print "Waiting 30 minutes and hoping problem is resolved."
            sys.stdout.flush()
            time.sleep(1800)
            continue

        #just do those not yet in database
        if args.checkdb:
            fnames=[fname.split('/')[-1].split('.')[0] for fname in remotefnames]
            spnames=list(js.SciencePointing.objects.filter(filename__in=fnames).values_list('filename',flat=True))
            for fname in remotefnames:
                if fname.split('/')[-1].split('.')[0] in spnames:
                    donefnames.append(fname)

        print "Found {0} new files on remote host".format(len(remotefnames))
        remotefnames = [fn for fn in remotefnames if fn not in donefnames]
        print "working on %s remote files"%len(remotefnames)
        print remotefnames

        if not idle: # or len(remotefnames) > 0:
            print "STAP_ingest.py:  running main loop at local time ",
            print time.asctime(time.localtime(start_loop_time))
            print "Found {0} new files on remote host".format(len(remotefnames))
        sys.stdout.flush()

        # For each new file:
        for fn in remotefnames:

            # Set up a download process for the file
            if Qdown.full():  break
            id = os.path.splitext(os.path.basename(fn))[0]
            workflow = STAP_ingest_download(name=id)
            Qargs = (workflow, acct, remotepath, fn)
            status = Qdown.start_job(*Qargs)
            if status == 0:
                print "Started download for", id
                donefnames.append(fn)
            elif status == 1:
                print "Couldn't start download for", id, " -- will retry"
            elif status == -1:
                print "Fatal error setting up for job", id, "-- won't retry"
                donefnames.append(fn)
            sys.stdout.flush()

        # For each file that's done downloading:
        while not Qdown.empty():

            # Don't submit any more amp-splits if Qsplit is full.
            if Qsplit.full():  break

            # Ok, so try grabbing a finished download job
            wf = Qdown.finish_job()
            job_out = wf.output
            if "stages" not in job_out or len(job_out["stages"]) < 1:
                print "Whoops, something REALLY weird with this job."
                print "It bombed so badly its output is hosed"
                continue
            elif "exception" in job_out["stages"][-1]:
                print "Download for", job_out["name"], "threw an exception:"
                print job_out["stages"][-1]["exception"]
                sys.stdout.flush()
                continue
            else:
                print "Download for", job_out["name"], "succeeded"

            # Assume everything's ok and try starting the split process
            id = job_out["name"]
            lfn = job_out["stages"][-1]["results"]
            # Construct a log file name based on the job name we're going to
            # want to give this thing below.
            logfile = id.replace("Skymapper", "ing" + runtag) + ".log"
            workflow = STAP_ingest_splitamp(name=id, workdir=scratch+"/"+id,
                                            logfile=logfile)
            Qargs = (workflow, lfn, manifest, today_utc)
            status = Qsplit.start_job(*Qargs)
            if status == 0:
                print "Started amp-split process for", id
            elif status == 1:
                print "Couldn't start amp-split for", id,
            elif status == -1:
                print "Fatal error setting up for job", id, "-- won't retry"
            sys.stdout.flush()

        # Finally, for each file that's done splitting:
        while not Qsplit.empty():

            wf = Qsplit.finish_job()
            job_out = wf.output
            # FY - add status code
            if "exception" in job_out["stages"][-1]:
                exit_status = js.PipelineExitStatus.objects.get_or_create(
                    description=job_out["stages"][-1]["exception"],
                    stage_fail=job_out["stages"][-1]["name"])[0]
            else:
                exit_status = js.PipelineExitStatus.objects.get_or_create(
                    description="Success")[0]
            
            # Register the job under some job name so we can track what
            # happened to each field that's ingested.
            bfn = job_out["name"]
            id = bfn.replace("Skymapper", "ing" + runtag)
            if "stages" not in job_out or len(job_out["stages"]) < 1:
                print "Whoops, something REALLY weird with this job."
                print "Job {0} bombed so badly its output is hosed".format(id)
                continue
            elif "exception" in job_out["stages"][-1]:
                # RS 2012/07/06:  Even if the WCS failed, we'll want some
                # idea of where the telescope was pointing.  Go ahead and
                # initialize on nearest SkyMapper field to given RA, DEC.
                # The fact that the WCS failed will show this information
                # to be potentially unreliable.
                lfn = job_out["args"][0]
                try:
                    print "WCS failed for", bfn,
                    print "-- trying to guess field from header (RA, DEC)."
                    head = pyfits.getheader(lfn)
                    coo = ephem.Equatorial(head['RA'], head['DEC'])
                    RR = float(coo.ra*180/ephem.pi)
                    DD = float(coo.dec*180/ephem.pi)
                    field = SMfield_lookup(RR, DD).id
                except:
                    print "Nope, couldn't get FIELD_ID from header of", bfn
                    field = 0
            else:
                lfn, field = job_out["stages"][-1]["results"]

            # Whether or not the job succeeded, we should register it in the
            # database of jobs so we can look at it again later.
            try:
                meta = Filenames(lfn, runtag=runtag)
                meta.field = field  #intended or actual
                job_out["name"] = id
                print "Registering image {0} under job {1}".format(bfn, id)
                newpt=register_pointing(meta,proc_status=exit_status)
                #FY - also register images that actually get copied back
                # only if job successful
                if not "exception" in job_out["stages"][-1] and newpt:
                    for ccd in range(1,33):
                        newname=FilenamesNew.absolute_dir(newpt.field.id,newpt.filter,ccd)+'/'+newpt.filename+'_%d.fits'%ccd
                        if os.path.exists(newname):
                            names=FilenamesNew(newname)
                            if ccd==17:
                                image=register_image(names,pointing=newpt,updatepointing=True)
                            else:
                                image=register_image(names,pointing=newpt,updatepointing=False)
                #not actually registering jobs at the moment
                #register_job(meta, job_out)
            except Exception as e:
                print "Warning:  exception thrown while trying to register job"
                exc_type, exc_value, exc_traceback = sys.exc_info()
                traceback.print_exception(exc_type, exc_value, exc_traceback,
                                          limit=100, file=sys.stdout)

            # Clean out raw data files
            if "exception" in job_out["stages"][-1]:
                print "Amp-split for", job_out["name"], \
                    "threw an exception in stage", job_out["stages"][-1]["name"], ":"
                print job_out["stages"][-1]["exception"]
                print "NOT unlinking", lfn, "so you can check why it failed."
            else:
                print "Amp-split for", job_out["name"], "succeeded"
                #if not args.local:
                print "Unlinking", lfn
                os.unlink(lfn)
            sys.stdout.flush()

        # RS 2013/09/16:  Man, do we really need to do all these cross checks?
        idle = (len(remotefnames) == 0 and Qdown.empty() and Qsplit.empty()
                and Qdown.N_jobs == 0 and Qsplit.N_jobs == 0)
        if idle and not args.production:
            break
        # This is the end of the loop; sleep till the next one
        end_loop_time = time.time()
        sleeptime = loop_delay - (end_loop_time-start_loop_time)
        if sleeptime > 0:  time.sleep(sleeptime)

    print "Ran out of stuff to do; exiting.\n"


class STAP_ingest_download(Pipeline.Workflow):
    """Runs a download process.  Really pretty simple."""

    def _setup_workflow(self, acct, remotepath, fn):
        self.banner = False
        rfn = remotepath + "/" + fn
        lfn = CPP.raw + "/" + fn
        self._stagelist = [["ingest_scp", scp_from_remote, [acct, rfn, lfn]],
#                           ["cleanup_header", cleanup_image_header, [lfn]],
                           ]


class STAP_ingest_splitamp(Pipeline.Workflow):
    """Runs an amp-splitting process.  A little more involved."""

    def _setup_workflow(self, lfn, manifest, date):
        self.banner = False
        self.lfn = lfn
        self._copy_over = [CPP.etc + "/" + fn for fn in CPP.SEx_files]
        self._copy_over += [CPP.etc + "/" + CPP.fsatconf, CPP.etc + "/" + CPP.fsatpar]
        self._copy_back = [CPP.raw + "/ingest_logs/" + date + "/" + self.logfile]
        try:
            self.imgtype = Filenames(lfn).imgtype
        except:
            raise IngestFail(msg="Couldn't initialize image metadata")
        # Don't apply WCS to flat-field images!
        if self.imgtype == 'flat':
            self._stagelist = \
            [
                ["ingest_fixheader", correct_image_header, [lfn, manifest]],
                ["ingest_splitamp", run_fits64to32_cal, [lfn]],
                ["ingest_addfieldccd", add_fieldccd_to_header, [lfn]],
            ]
        elif self.imgtype == 'object':
            self._stagelist = \
            [
                ["ingest_fixheader", correct_image_header, [lfn, manifest]],
                ["ingest_splitamp", run_fits64to32_cal, [lfn]],
                ["ingest_applywcs", apply_rough_WCS, [lfn]],
                ["ingest_addfieldccd", add_fieldccd_to_header, [lfn]],
            ]
        else:
            pass
        
    def _cleanup_workflow(self):
        for ccdfits in glob.glob(self.workdir + "/??/*.fits"):
        # for ccdfits in glob.glob(self.workdir + "/17/*.fits"):
            if '_solved.fits' in ccdfits or '_mask.fits' in ccdfits:
                continue
            if self.imgtype == "object":
                names = FilenamesNew(ccdfits)
            elif self.imgtype == "flat":
                names = FilenamesNewflat(ccdfits)
            else:
                # RS 2013/09/13:  If it's not a flat or a science frame,
                # we're not set up to deal with it yet.
                continue
            self._copy_back.append(names.absfname)
            dotdot = os.path.dirname(ccdfits) + "/.."
            # RS 2014/06/26:  Modified to use STAP_callexternal to return
            # status codes and print error messages.
            if doage("mv {0} {1}".format(ccdfits, dotdot))[0] != 0:
                print "PROBLEM moving", ccdfits, "after generation"
            if self.imgtype == "object":
                sexfits = ccdfits + ".stars"
                if os.path.exists(sexfits):
                    self._copy_back.append(names.absfname + '.stars')
                    if doage("mv {0} {1}".format(sexfits, dotdot))[0] != 0:
                        print "PROBLEM moving", sexfits, "after generation"
                maskfits=ccdfits.replace('.fits','_mask.fits')
                if os.path.exists(maskfits):
                    if doage("gzip {0}".format(maskfits))[0] != 0:
                        print "PROBLEM gzipping", sexfits
                        self._copy_back.append(names.absfname.replace('.fits','_mask.fits'))
                    else:
                        maskfits = maskfits + '.gz'
                        self._copy_back.append(names.absfname.replace('.fits','_mask.fits.gz'))
                        if doage("mv {0} {1}".format(maskfits, dotdot))[0] != 0:
                            print "PROBLEM moving", maskfits, "after generation"


# ============================================================================
#                       RS helper functions put here
# ============================================================================


def doage(cmd):
    """Quickly run a command and retrieve its stdout and stderr"""
    print cmd
    sys.stdout.flush()
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                                             stderr=subprocess.PIPE)
    (stdoutstr, stderrstr) = proc.communicate()
    if proc.returncode != 0:  print stderrstr
    sys.stdout.flush()
    return (proc.returncode, stdoutstr, stderrstr)

def fetch_schedule_logs(date):
    """Fetch scheduler logs from remote system and return list of filenames"""
    # make local directory
    local_scheduler_log_path = CPP.raw + "/scheduler_logs"
    #hack
    #return 0, glob.glob("{0}/SM_LOGFILE.{1}*".format(
    #        local_scheduler_log_path, date))
    doage("mkdir -p " + local_scheduler_log_path)
    # copy files on remote system over
    remote_scheduler_log_path = \
        "skymap@markle:/home/skymap/Scheduler/Scheduler_4/TEST_REALLIFE_1"
    status, stdoutstr, stderrstr = doage("scp {0}/SM_LOGFILE.{1}\* {2}".format(
            remote_scheduler_log_path, date, local_scheduler_log_path))
    # see what we copied
    return status, glob.glob("{0}/SM_LOGFILE.{1}*".format(
                             local_scheduler_log_path, date))

def ls_remote(acct, remotepath, obstype='MS_*|TFF_*|3PS_42|BAD_42',imageid='53|89|805|053|089|169'):
    """Determine which files need to be ingested.

        acct:  e.g., "skymap@markle", or None if image is local
        remotepath:  e.g., /export/markle2/images/2012-05-31
        obstype: obstype_surveyId as decoded by decode-Ob_Id.pl
    """
    ismarkle=False
    if acct: 
        # RS 2013/09/21:  Modified to make sure TAROS is done writing the file
        # to disk before attempting to download!  Incomplete files will not
        # have most keywords filled yet; we pick NCCDS here.  Using gethead on
        # files without the requested keyword will print nothing, at least on
        # markle.  I wouldn't swear to this code being portable!
        if re.search("markle", acct):
            ismarkle=True
            cmd = "ssh {0} \'ls -tr {1}/*.fits | xargs gethead -au NCCDS"\
                .format(acct, remotepath)
            # RS 2014/06/26:  In fact it isn't portable.  If we're trying to get
            # data off of raijin, gethead won't be available, but it won't matter
            # because raijin will only have data that are completely written out.
            # i.e., processing archived data is not real-time, by definition.
        else:
            cmd = "ssh {0} \'cd {1}; ls *.fits".format(acct, remotepath)

        cmd += " | awk '\\''{if ($2!=\"___\") print $1}'\\''"
        if ismarkle and len(obstype)>0:
            #add awk command to call decode-Ob_Id.pl                                                                                                                         
            matchstrs=[]
            for obstypeid in obstype.split('|'):
                parts=obstypeid.split('_')
                if parts[1]=='*':
                    matchstrs.append("/Observation type:[ ]+"+parts[0]+"/")
                else:
                    matchstrs.append("/Observation type:[ ]+"+parts[0]+".*[ ]+Survey id:[ ]+"+parts[1]+"/")
            awkstr=" || ".join(matchstrs)
            cmd +=" | awk -F'\\''_'\\'' '\\''{print $0} {system(\"/home/skymap/Scheduler/Scheduler_4/utils/decode-Ob_Id.pl \"$2)}'\\''"
            cmd +=" | xargs -L 6"
            cmd +=" | awk '\\''"+awkstr+" {print $1}'\\''"
            cmd +="\'"
        else:
            cmd +="\'"
    else:
        # RS 2014/06/26:  This is for data already downloaded locally.
        cmd = "ls -tr {0}/*.fits | xargs gethead -au NCCDS".format(remotepath)
        cmd += " | awk '{if ($2!=\"___\") print $1}'"

                    
    status, stdoutstr, stderrstr = doage(cmd)
    if len(obstype)>0 and ismarkle:
        filelist = [line for line in stdoutstr.strip().split("\n")
                    if re.search('Skymapper_.*T(\d{2})', line)]
    elif len(obstype)>0:
        filelist=[]
        obstypes=obstype.split('|')
        obstypes1=[obst.split('_')[0] for obst in obstypes if obst.split('_')[1]=='*']
        obstypes2=[obst for obst in obstypes if obst.split('_')[1]!='*']
        for line in stdoutstr.strip().split("\n"):
            if not re.search('Skymapper_.*T(\d{2})', line): continue
            obscode=decode_obs_id(line.split('/')[-1].split('_')[1])
            if (obscode[0]+'_%d'%obscode[1]) in obstypes2 or obscode[0] in obstypes1:
                filelist.append(line)
    else:
        # Return the individual filenames we found.
        # The regex below should match:
        filelist = [line for line in stdoutstr.strip().split("\n")
                    # FY - since June 2014, flat fiels changed from 11 to 53
                    # original SN science test data = 805, flat fields = 53,
                    # main survey data = 87, SN survey data = 89,
                    # test WCS set for Chris Wolf (2013/03/18) = 133
                    if re.search('Skymapper_('+imageid+').*T(\d{2})', line)]

    return status, filelist[::]


# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#                           Workflow stage functions
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .


def scp_from_remote(acct, rfn, lfn):
    """Copy a file from the remote server"""
    # Does the file already exist on local disk?
    if os.path.exists(lfn):  return lfn
    # Do we have enough disk space left?
    diskcheck = os.statvfs(CPP.raw)
    diskfree = (diskcheck.f_frsize*diskcheck.f_bavail)*bytes_to_mb
    if diskfree < min_space_avail:
        raise IngestFail(msg="not enough space left on local disk")
    # We don't already have it and we have enough space, so try to copy
    if acct:  cmd = "scp {0}:{1} {2}".format(acct, rfn, lfn)
    else:  cmd = "cp {0} {1}".format(rfn, lfn)
    status = doage(cmd)[0]
    if status != 0:
        raise IngestFail("scp failed")
    temp=os.system('chmod 644 %s'%lfn)
    return lfn


def cleanup_image_header(lfn):
    """Clean up empty or duplicated keywords."""
    #For now, just the DUMMY placeholder for SkyDice
    with pyfits.open(lfn,mode='update') as hdulist:
        header=hdulist[0].header
        if 'DUMMY' in header: del header['DUMMY']
        hdulist.flush()

    return lfn


def correct_image_header(lfn, manifest):
    """Corrects image header quantities to agree with the scheduler logs."""
    default_val = { "FILTNAME": '',
                    "ROTSKYPA": 0.0,
                    "RA"      : '00:00:00.00',
                    "DEC"     : '+00:00:00.0', }
    # RS:  Only the obsseq is guaranteed to be a unique ID; the survey ID
    # depends on things the user puts in that may not be verified.
    # Since not all manifest images will be on disk, the usual Filenames
    # metadata wrapper (which reads the FITS header) may not work.
    # So we just grep the obsseq from the filename.
    basename = os.path.basename(lfn)
    match = re.search("(\d+)_(\d+)_(\d+-\d+-\d+T\d+:\d+)", lfn)
    if match == None:
        match = re.search("(\d+)_(\d+-\d+-\d+T\d+:\d+)", lfn)
        if match == None:
            print "correct_headers:  can't extract obsseq from", basename
            return
        else:
            id = "{0}_{1}".format(match.group(1), match.group(2))
    else:
        id = "{0}_{1}".format(match.group(1), match.group(3))
    print "correct_headers:  Looking for id", id, "in manifest"
    if id not in manifest:
        print "WARNING:  {0} NOT found in manifest".format(id)
        print_manifest(manifest)
        return
    # This used to be in a "try" clause, but now the pipeline will catch
    # any exceptions, so just go ahead with it.
    with pyfits.open(lfn, mode='update') as hdulist:
        header = hdulist[0].header
        print "correct_headers:  Checking", basename, "against manifest"
        for kw in ['FILTNAME','RA','DEC','ROTSKYPA']:
            if kw not in manifest[id]:  continue
            if kw not in header:  header.update(kw, default_val[kw])
            if header[kw] == manifest[id][kw]:  continue
            # If the keyword value disagrees with the manifest, correct it
            if (kw == 'RA'  and header[kw][0:6] == manifest[id][kw][0:6] or
                kw == 'DEC' and header[kw][0:7] == manifest[id][kw][0:7] or
                kw == 'ROTSKYPA'
                    and abs(float(header[kw])-float(manifest[id][kw]))<0.1):
                continue
            print "correct_headers:  corrected {0} from {1} to {2} in {3}"\
                  .format(kw, header[kw], manifest[id][kw], basename)
            header.add_history(" correct_headers: "
                               "corrected {0} from {1} to {2}".format
                               (kw, header[kw], manifest[id][kw]))
            header.update(kw, manifest[id][kw])
        sys.stdout.flush()

def run_fits64to32_cal(lfn, workdir="."):
    """Run fits64to32_cal on an image."""
    names = Filenames(lfn)        
    calbin = CPP.bin + "/fits64to32_cal"
    # RS 2013/09/14:  Looks like FY took out the "-align" here.
    # cmd = "{0} {1} {2} -align".format(calbin, lfn, workdir)
    cmd = "{0} {1} {2} ".format(calbin, lfn, workdir)
    # If this image is a flat field, we don't need to flatten it
    if names.imgtype == 'flat':
        print "This is a flat field image, no need to flatten it!"
    # If it's a science image, try to flatten it
    elif names.imgtype == 'object':
        #full-well saturation mask
        fsatname = names.basefname + ".fsat"
        fsatcmd = "sex {0} -c {1} -PARAMETERS_NAME {2} -CHECKIMAGE_NAME {3}".format(lfn, CPP.fsatconf,CPP.fsatpar,fsatname)
        status, stdoutstr, stderrstr = doage(fsatcmd)
        if status == 0:
            cmd = "{0} -fsatmap {1}".format(cmd,fsatname)

        cmd = "{0} -domask".format(cmd)
        filtname = names.filter
        # Find the flat field for this filter, if it exists.
        # If there are several, choose the most recently built flat.
        flatpath = CPP.flats + '/' + filtname
        files = glob.glob(flatpath + "/17/flat_*_17.fits")
        doflat=False
        if len(files) >=1:
            #find closest flat
            flatoff=[]
            for file in files:
                try:
                    mm=re.search("/\d+/flat.*_(\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2})_.*17.fits",file)
                    timeoffset=datetime.datetime.strptime(mm.group(1),'%Y-%m-%dT%H:%M:%S')-names.dateobs
                    flatoff.append(abs(timeoffset.days+timeoffset.seconds/24./3600.))
                except:
                    flatoff.append(100)
            if min(flatoff)<100:
                minflat=flatoff.index(min(flatoff))
                file=files[minflat]
                            
                mm = re.search("/\d+/(flat.*)_17.fits", file)
                flatbase = mm.group(1)
                # If flats exist for all the CCDs, call fits64to32 with flat fields.
                flatfnames = ["{0}/{1:02d}/{2}_{3:d}.fits".format(flatpath,ccd,flatbase,ccd) for ccd in range(1,33)]
                # print "Checking flat", flatbase
                flatsfound = [os.path.exists(fn) for fn in flatfnames]
                if all(flatsfound):
                    doflat=True
                else:
                    print "Missing/incomplete flat field in", flatbase
                    # If one or more CCDs is/are missing flats for this field and filter,
                    # don't flatten the image, just overscan-subtract and align amps.
        if doflat:
            cmd = "{0} -flatpath {1} -flatbase {2}".format(
                cmd, flatpath, flatbase)
        else:
            print "No proper flat fields available for filter", filtname
            print "Splitting without applying flatfield"
            cmd = "{0} -align".format(cmd)

        # Find the bad pixel mask, if it exists.
        # If there are several, choose the most recently built.
        inmaskpath = CPP.masks + '/' + filtname
        files = glob.glob(inmaskpath + "/17/badpix_*_17.fits")
        if len(files) < 1:
            print "No bad pixel mask available for filter", filtname
        else:
            files.sort()
            files.reverse()
            for file in files:
                try:
                    mm=re.search("/\d+/badpix_(\d{4}-\d{2}-\d{2})_17.fits",file)
                    fobsseq=mm.group(1)
                    if fobsseq < names.obsseq: break
                except:
                    pass

            mm = re.search("/\d+/(badpix.*)_17.fits", file)
            inmaskbase = mm.group(1)
            # If bad pixel masks exist for all the CCDs, call fits64to32 with
            # input mask
            inmasknames = ["{0}/{1:02d}/{2}_{3:d}.fits".format(inmaskpath,ccd,inmaskbase,ccd) for ccd in
                           range(1,33)]
            masksfound = [os.path.exists(fn) for fn in inmasknames]
            if all(masksfound):
                cmd = "{0} -maskpath {1} -maskbase {2}".format(
                    cmd, inmaskpath, inmaskbase)
            else:
                print "Missing/incomplete bad pixel masks in", inmaskbase

    # If it's any other kind of image, or IMAGETYP is corrupt, pitch it.
    else:
        raise IngestFail("unknown image type")

    # Run the external command and handle problems if they arise.
    status, stdoutstr, stderrstr = doage(cmd)
    if status == 10:
        print "WARNING:  image might be empty (shutter problem?)"
        print "Some amps might be ok, so I'll try to copy those back."
    elif status != 0:
        print "----------"
        print stdoutstr.strip()
        print stderrstr.strip()
        print "----------"
        raise IngestFail("run_fits64to32_cal failed to process image")

def apply_rough_WCS(lfn, workdir="."):
    """Apply tangent-plane WCS solution to all 32 CCDs of a mosaic"""
    basename = os.path.splitext(os.path.basename(lfn))[0]
    cmd = CPP.bin + '/SM-ASTROMETRY.py {0} -impath {1}'.format(
            basename, workdir)
    status, stdoutstr = STAP_callexternal(
            cmd, combinestderr=True, getstdout=True)
    # RS 2014/07/10:  Print SM-ASTROMETRY.py output regardless of whether
    # it fails or not.  This will help Main Survey debugging efforts.
    print "----------"
    print stdoutstr.strip()
    print "----------"
    sys.stdout.flush()
    # RS 2012/06/19:  This way of doing it is kind of hackneyed, but quickest
    # since SM-ASTROMETRY.py is being treated here as an external program.
    # If SM-ASTROMETRY.py failed, parse the exit status code to find out why.
    if status != 0:
        if status == -1:
            msg = "Invalid/corrupt image input to SM-ASTROMETRY.py"
        elif status == -2:
            msg = "SExtractor failed on image"
        elif status == -3:
            msg = "Not enough objects found on image to solve for WCS"
        elif status == -4:
            msg = "astrometry.net terminated abnormally on image"
        elif status == -5:
            msg = "astrometry.net failed to solve image"
        else:
            msg = "Unknown exception thrown in SM-ASTROMETRY.py (status = {0})"\
                    .format(status)
        raise IngestFail(msg)

def add_fieldccd_to_header(lfn):
    """Add some crucial FITS keywords to each CCD header"""
    field = 0
    bfn = os.path.splitext(os.path.basename(lfn))[0]
    for ext in range(1,33):
        # RS:  Is this extension empty?  (the others might not be)
        # We can tell by checking whether fits64to32_cal finished its loop.
        src = "{ccd:02d}/{basefn}_{ccd:d}.fits".format(ccd=ext, basefn=bfn)
        # Open the FITS header
        with pyfits.open(src, mode='update') as hdulist:
            head = hdulist[0].header
            imagetype = head['IMAGETYP']
            # Sanity checks
            # RS 2013/09/14:  Since FY turned the "-align" option off,
            # even good object frames won't have this keyword anymore.
            # if 'AMPSCL_L' not in head:                
            #     print "Chucking {0}, it seems to be empty".format(src)
            #     os.unlink(src)
            #     continue
            # RS 2013/09/14:  If this isn't a science frame, it won't have an
            # associated field and there's no point in going further.
            # elif imagetype != 'object':
            if imagetype != 'object':
                continue
            elif not all([k in head for k in ('ROTSKYPA','CRVAl1','CRVAL2')]):
                print "Chucking {0}, WCS failed for some reason".format(src)
                os.unlink(src)
                continue
            # RS 2012/06/18:  We should only be here for 'object' frames.
            # Look up the field in the SkyMapper field config information
            crval1, crval2 = float(head['CRVAL1']), float(head['CRVAL2'])
            field = SMfield_lookup(crval1, crval2).id
            # Work out which SUBFIELD we're in.  SUBFIELD is meant to be a
            # rotationally-invariant CCD number on the sky for each pointing.
            # It will equal the CCD number if ROTSKYPA = 0, and will be
            # something else (see mapping in Constants) for ROTSKYPA = 180.
            flip = int(math.cos(head['ROTSKYPA']*math.pi/180)*1.001)
            # ROTSKYPA = 0
            if flip == 1:  subfield = ext
            # ROTSKYPA = 180
            elif flip == -1:  subfield = Imager.ccd_map[ext]
            # some kind of error
            else:
                print "Chucking {0}, WCS failed for some reason".format(src)
                os.unlink(src)
                continue
            # Add FIELD_ID and SUBFIELD to header
            head.update('FIELD_ID', field,
                        "Unique SkyMapper field ID used by scheduler")
            head.update('SUBFIELD', subfield,
                        "Subfield (SkyMapper CCD for ROTSKYPA=0)")
            # FY: add skymapper observation type
            try:
                head.update('PROGRAM',decode_obs_id(bfn.split('_')[1])[0],
                            "SkyMapper observation program type")
            except:
                pass
        # RS 2012/07/12:  As a final step, for 'object' files rename so that
        # the extension is the subfield, not the physical CCD number.
        if imagetype == 'object' and ext != subfield:
            newsrc = src.replace("{ccd:d}.fits".format(ccd=ext),
                                 "{ccd:d}.fits".format(ccd=subfield))
            subprocess.call(["mv", src, newsrc])
            #rename the mask as well
            masksrc=src.replace('.fits','_mask.fits')
            if os.path.exists(masksrc):
                newmsrc = masksrc.replace("{ccd:d}_mask.fits".format(ccd=ext),
                                         "{ccd:d}_mask.fits".format(ccd=subfield))
                subprocess.call(["mv", masksrc, newmsrc])
            
    return lfn, field


# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#                          Manifest-related functions
# . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .


def manifest_from_scheduler_log(logfname, manifest={ }):
    """Creates manifest of images and image properties from scheduler log"""
    filter, ra, dec = 'indef', 'indef', 'indef'
    with open(logfname) as log:
        for line in log:
            # here's the line with what the scheduler told TAROS to do
            mm = re.search("(\S+ \S+) .* id (\d+); .* filter (.); rot (\d+);"
                           ".*RA/dec. (\S+) / (\S+) \((\S+) / (\S+)\)", line)
            if mm != None:
                obsseq = mm.group(1)[:19].replace(' ','T')
                obsid = int(mm.group(2))
                imgid = "{0}_{1}".format(obsid, obsseq[:-3])
                if imgid in manifest:  continue
                filter = mm.group(3).strip()
                rotskypa = float(mm.group(4))
                ra_str, dec_str = mm.group(5,6)
                ra_rad, dec_rad = float(mm.group(7)), float(mm.group(8))
                # normalize format of RA and DEC
                ra_groups = [float(f) for f in ra_str.split(":")]
                ra  = "{0:02.0f}:{1:02.0f}:{2:04.1f}".format(*ra_groups)
                dec_groups = [float(f) for f in dec_str.split(":")]
                if dec_rad > 0:
                    dec_groups = ['+'] + dec_groups
                else:
                    dec_groups = ['-'] + dec_groups
                    dec_groups[1] = abs(dec_groups[1])
                dec = "{0}{1:02.0f}:{2:02.0f}:{3:02.0f}".format(*dec_groups)
                # add to manifest
                print "Adding FILTNAME etc. of {0} to manifest".format(imgid)
                manifest[imgid] = { 'FILTNAME' : filter,
                                    'ROTSKYPA' : rotskypa,
                                    'RA'       : ra,
                                    'DEC'      : dec,
                                    'ra_rad'   : ra_rad,
                                    'dec_rad'  : dec_rad }
            # here's the line with the name of the corresponding image
            mm = re.search("(\S+ \S+).*TAROS response: DONE (\d+).*"
                           "(Skymapper.*)\.fits", line)
    return manifest

def print_manifest(manifest):
    files = manifest.keys()
    files.sort()
    print "Manifest contents:"
    for file in files:
        header = manifest[file]
        print file, header["RA"], header["DEC"], header["FILTNAME"]


# ============================================================================
#                    This file is of course executable
# ============================================================================


if __name__ == "__main__":
    main()
