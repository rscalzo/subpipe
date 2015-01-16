#!/usr/bin/env python

# ============================================================================
# RS 2012/02/20:  Implementation of Pipeline class hierarchy
# ----------------------------------------------------------------------------
# This code defines an interface for a modular pipeline.  It covers all
# the stuff we'd know and love pipelines to do, including:
# -- "Each stage [of the pipeline] should know its inputs and outputs."
# -- "Each stage should know where to find its inputs and outputs."
#     (actually, the naming conventions etc. belong in the Constants module)
# -- "Each stage should be able to report whether it was successful,
#     and which, if any, of its outputs it produced."
# -- "Each stage, on failure, should catch all exceptions and pass back
#     terse error messages to the main process." (possibly Django-ized.)
# -- "Each stage should be able to time itself, and profile itself if it's
#     our Python code, so we can track performance."
# -- "Each stage should be able to run stand-alone with appropriate command
#     line arguments, (just separate __main__() from algorithm)."
# The main pipeline process/class should be able to connect these components
# in a modular/automatic way to form a pipeline with as many stages as needed
# for the user to do any desired (re-)processing workflow.
# ============================================================================


import os
import sys
import time
import datetime
import traceback
import subprocess
import processing
from .TrackableException import TrackableException


class Workflow(object):
    """Workflow class for a sequence of steps which can be queued in parallel

    This represents a list of executable functions.  Each stage is run in
    sequence, using the output of each to determine the input of the next.
    Although we can know the functions we're going to call in advance,
    it's easier to determine the (possibly complicated) mappings between
    inputs and outputs to different stages at run time.  We therefore put
    the work of defining the linkages to a method called _setup_workflow,
    which is called at the beginning of the "run" method and defines these
    attributes (empty lists by default):

        self._stagelist:  List of pipeline stages to run, in sequence.
            Items:  [name, func, args, kwargs]
            When executed, calls func(*args,**kwargs).
        self._copy_over:  List of files to copy to working directory
            before starting workflow (optional).

    We also allow the definition of a method called _cleanup_workflow,
    which is called at the end of the "run" method and defines this attribute
    (empty list by default):

        self._copy_back:  List of files to copy from working directory
            after completing workflow (optional).

    Note that we define _copy_back only after the workflow has run, since the
    names of many of the files to copy back may not be known in advance and
    may not even be well-described as a glob.
    """

    # TODO:  I don't really like the division of labor between __init__ and
    # _setup_workflow; it doesn't seem quite natural.  Why not have it all
    # in one or the other?  Will have to think a bit.  I'm not sure it makes
    # that much difference, it just bugs me.
    def __init__ (self, name = "Generic Workflow",
                logfile = None, workdir = ".", banner=True):
        self.name = name
        self.workdir = workdir
        self.logfile = logfile
        self.banner = banner
        self._stagelist = [ ]
        self._copy_over = [ ]
        self._copy_back = [ ]

    def _setup_workflow(self, *args, **kwargs):
        pass

    def _cleanup_workflow(self, stage_failed=None):
        pass

    def run(self, *args, **kwargs):
        """Runs the workflow

        Pipeline.Workflow output is a dict, structured as follows:

            name:       Name of the workflow.
            workdir:    Working directory of the workflow.
            copy_back:  Mirror of self._copy_back so that the Queue
                        or user knows the names of important output files.
            stages:     A list of dicts, one for each stage, with keywords:
            {
                name:       Name of the stage.
                t_start:    Timestamp for start of the stage.
                t_finish:   Timestamp for stage finishing.
                usr_time:   User time elapsed from start to finish.
                cpu_time:   CPU time elapsed from start to finish.
                walltime:   Wall time elapsed from start to finish.
                results:    Return value of func, if any.
                exception:  Exception message raised during execution, if any.
                exc_type:   Type of exception raised during execution, if any.
                traceback:  Traceback for terminating exception, if any.
            }
        """

        stdout_tmp, stderr_tmp = sys.stdout, sys.stderr
        self.output = { "name": self.name,
                        "workdir": self.workdir,
                        "logfile": self.logfile,
                        "args": args,
                        "kwargs": kwargs,
                        "copy_back": [ ],
                        "stages": [ ], }
        try:
            self._setup_workflow(*args, **kwargs)
        except Exception as e:
            print "Problem running setup phase for workflow", self.name
            exc_type, exc_value, exc_traceback = sys.exc_info()
            msg = "{0}: {1}".format(exc_type.__name__, exc_value)
            self.output["stages"].append({ "name" : "setup",
                                           "exception": msg })
            return self.output

        # If a working directory was specified, change directory there
        os.chdir(self.workdir)
        # Redirect stdout to the logfile if one was specified
        if self.logfile:  sys.stdout = open(self.logfile,"wt")

        stage_failed = None
        for stage in self._stagelist:
            # Time the stage.  Also consider "import resource" to do this:
            # http://docs.python.org/release/2.5.2/lib/module-resource.html
            t_start = datetime.datetime.now()
            start_times = os.times()
            self.output["stages"].append({ "name" : stage[0] })
            # Catch exceptions so we can exit gracefully if we're running
            # this thing in parallel pipeline mode.
            try:
                if self.banner:  self.print_banner(stage[0])
                if len(stage) > 3:  stage_kwargs = stage[3]
                else:  stage_kwargs = { }
                results = stage[1](*(stage[2]),**(stage_kwargs))
                self.output["stages"][-1]["results"] = results
            except:
                exc_type, exc_value, exc_traceback = sys.exc_info()
                traceback.print_exception(exc_type, exc_value, exc_traceback,
                                          limit=100, file=sys.stdout)
                print "{0}:  failed at stage {1}".format(self.name, stage[0])
                exc_name = exc_type.__name__
                if issubclass(exc_type, TrackableException):  msg = exc_value
                else:  msg = "unexpected {0}, see log file".format(exc_name)
                msg = "{0}: {1}".format(exc_name, msg)
                self.output["stages"][-1]["exception"] = msg
                stage_failed = stage[0]
            end_times = os.times()
            t_finish = datetime.datetime.now()
            # RS 2012/06/18:  Get *children's* user and system time -- looks
            # like this relates to child *processes*, which is what we want.
            dt_usr  = end_times[2] - start_times[2]
            dt_sys  = end_times[3] - start_times[3]
            dt_wall = end_times[4] - start_times[4]
            self.output["stages"][-1].update({ "t_start"  : t_start,
                                               "t_finish" : t_finish,
                                               "usr_time" : dt_usr,
                                               "cpu_time" : dt_sys,
                                               "walltime" : dt_wall, })
            if self.banner:
                print
                print "Elapsed user time in {0}:  {1:.3f} sec".format(
                        stage[0], dt_usr)
                print "Elapsed CPU time in {0}:  {1:.3f} sec".format(
                        stage[0], dt_sys)
                print "Elapsed wall time in {0}:  {1:.3f} sec".format(
                        stage[0], dt_wall)
                print; print
            sys.stdout.flush()
            if stage_failed:  break

        # Stuff might go wrong in the clean-up stage too, for example the
        # image compression stage.  So trap and report those as well.
        try:
            if self.banner:  self.print_banner("cleanup")
            self._cleanup_workflow()
        except Exception as e:
            print "Problem running cleanup phase for workflow", self.name
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_traceback,
                                      limit=100, file=sys.stdout)
        self.output["copy_back"] = self._copy_back
        # Close the logfile, if one was opened
        if self.logfile:  sys.stdout.close()
        sys.stdout, sys.stderr = stdout_tmp, stderr_tmp
        # Return entire Workflow object with filled "output" attribute
        return self

    def print_banner(self, stagename):
        """Prints a banner in the log output with the stage name"""
        now = time.asctime(time.localtime(time.time()))
        mybanner = "Executing stage {0} at {1}".format(stagename,now)
        print "+{0:-^76s}+".format("")
        print "|{0: ^76s}|".format(mybanner)
        print "+{0:-^76s}+".format("")
        print
        sys.stdout.flush()


class Queue(object):
    """Parallel queue class for workflows including file management

    This represents the main pipeline process.  Its job is as follows:

        Copy the input files from wherever they are into working_dir.
            In the context of routine pipeline operations this will be a
            subpath of /ramdisk, but one could run it in '.' instead for
            debugging purposes.
        Run a specified Pipeline.Workflow in a processing.Queue.
        Collect the outputs from working_dir and copy them to their
            final destinations on disk.

    All input and output files must be specified by their absolute initial
    and final pathnames.  It's not Pipeline.Queue's job to make sure that
    the files are where the user says they should be, since that would
    require detailed global knowledge of the problem being solved.

    Inputs:
        name            Name for the queue, if desired.
        processors      Number of processes to run in parallel.
        keep_workdir    Boolean:  on fail, don't clean up workdir
    """

    def __init__(self, name = "Generic Pipeline Queue",
                 processors = min(16,processing.cpuCount()),
                 keep_workdir=False):
        """Initializes the queue and starts worker threads"""
        processing.freezeSupport()
        self.name = name
        self.N_processors = processors
        self.N_jobs = 0
        # Sorting out names:  These Queues are processing.Queues...
        self.Qin = processing.Queue(self.N_processors)
        self.Qout = processing.Queue(self.N_processors)
        self.workers = [ ]
        self.keep_workdir = keep_workdir
        # ... while "worker" belongs to this class (Pipeline.Queue).
        for i in range(self.N_processors):
            self.workers.append(processing.Process
                    (target=Queue._worker, args=(self.Qin,self.Qout)))
            self.workers[-1].start()

    def __del__(self):
        """Kills the worker threads when the object is destructed"""
        print "{0}:  Terminating worker processes...".format(self.name)
        # for w in self.workers:
        #     self.Qin.put('STOP')
        for w in self.workers:
            w.terminate()
            w.join()
        print "{0}:  Joining queue threads...".format(self.name)
        self.Qin.close()
        self.Qin.jointhread()
        self.Qout.close()
        self.Qout.jointhread()

    @staticmethod
    def _worker(input, output):
        """A wrapper allowing Pipeline.Workflow instances to be queued"""
        # This is actually trickier than I thought at first.  To accomodate
        # different argument structures we need to parse the input to _worker.
        for job, args, kwargs in iter(input.get, 'STOP'):
            output.put(job.run(*args,**kwargs))

    def _copystuff(self, src, destdir, cleanup=False):
        """Copy files to/from "workdir" (e.g. local scratch)

        Copies everything it can without complaint.  Prints a warning only
        if some file system error occurs:  that is, if the file exists but
        the copy fails for some reason, or if the cleanup fails.
        If a workflow fails because one of its inputs wasn't copied, it'll
        raise an exception.  If a workflow fails for some other reason and
        doesn't produce one or more output files, we'll have as much of its
        output copied back as we can find.
        """
        src, destdir = os.path.normpath(src), os.path.normpath(destdir)
        # Check 1:  Does source file exist?
        if not os.path.exists(src):
            print "{0}:  warning, {1} doesn't actually exist".format(
                self.name, src)
            return 1
        # Check 2:  Does destination directory exist?  If not, mkdir -p it.
        if not os.path.exists(destdir):
            if subprocess.call(["mkdir", "-p", destdir]) != 0:
                print "{0}:  warning, couldn't mkdir -p {1}".format(
                   self.name, destdir)
                return 2
        # Check 3:  Is destdir the src dir?  (cp will fail if so)
        if os.path.samefile(os.path.dirname(src), destdir):
            print "{0}:  warning, {1} and {2} are the same dir".format(
                    self.name, os.path.dirname(src), destdir)
            return 3
        # Copy the file and check 4: that it succeeded
        if subprocess.call(["cp", src, destdir]) != 0:
            print "{0}:  warning, couldn't copy {1} to {2}".format(
                   self.name, src, destdir)
            # Check whether this is because we're out of disk space.
            # Define "out" as "less than 0.5 GB = 512 MB left".
            min_space_avail = 512
            bytes_to_mb = 1.0/(1024**2)
            diskcheck = os.statvfs(destdir)
            diskfree = (diskcheck.f_frsize*diskcheck.f_bavail)*bytes_to_mb
            if diskfree < min_space_avail:
                print "Looks like we've run out of space on the target disk."
            return 4
        # Optional:  Clean up behind us
        if cleanup:
            if subprocess.call(["rm", src]) != 0:
                print "{0}:  warning, couldn't clean up {1}".format(
                       self.name, src)
                return 5
        return 0

    def full(self):
        """Checks to see if the input queue is full"""
        # return self.Qin.full()
        return self.N_jobs >= self.N_processors

    def empty(self):
        """Checks to see if the output queue is empty"""
        return self.Qout.empty() # and self.N_jobs == 0

    def start_job(self, workflow, *inputs, **kwinputs):
        """Starts a job running in the queue
        
        First tries to mkdir -p workflow.workdir.  If this fails, there will
        be no working space, so just die.  Otherwise, copy input files to
        workflow.workdir and start the job in the queue.
        """
        if self.full():  return 1
        # TODO:  For now we have to call workflow._setup_workflow() here,
        # even though it's also called in workflow.run(), so that our
        # Pipeline.Queue knows what to copy *over*.  This is clumsy and what
        # we really want is for the Workflow to do the copying, since it knows
        # what it needs; but there still has to be a throttle on the rate of
        # data access.  Figure out a good solution to this.
        try:
            workflow._setup_workflow(*inputs, **kwinputs)
        except:
            # TODO:  Need to think of some better way to handle this...
            print "start_job:  problem setting up workflow..."
            exc_type, exc_value, exc_traceback = sys.exc_info()
            traceback.print_exception(exc_type, exc_value, exc_traceback,
                                      limit=100, file=sys.stdout)
            return -1
        try:
            status=0
            for file in workflow._copy_over:
                status = self._copystuff(file, workflow.workdir)
                if status > 1:  break
            if status == 2:
                print "start_job:  problem creating workdir; permissions?"
                return 1
            elif status == 4:
                print "start_job:  problem copying files over; scrapping"
                subprocess.call(["rm", "-rf", workflow.workdir])
                return 1
        except:
            print "start_job:  problem copying stuff over, disk may be full"
            subprocess.call(["rm", "-rf", workflow.workdir])
            return 1
        self.Qin.put((workflow, inputs, kwinputs), block=False)
        self.N_jobs += 1
        return 0

    def finish_job(self):
        """Remove a job's output from the queue
        
        After taking the job's output, try to copy things back.  The items
        in workflow._copy_back will all have absolute paths, so look for
        files with the same basenames in workflow.workdir.
        """
        if self.empty():  return { }
        workflow = self.Qout.get(block=False)
        # Since the original workflow object is out of this scope, the output
        # files to copy back have been mirrored in the process's output dict.
        for file in workflow._copy_back:
            dest, src = os.path.split(file)
            src = "{0}/{1}".format(workflow.workdir, src)
            status = self._copystuff(src, dest)
            if status > 3:
                print "copystuff failed to copy {0} to {1} with exit code {2}"\
                        .format(src, dest, status)
        # Clean up working directory on ramdisk (if not default value),
        # if no failure occurred or if self.keep_workdir was not specified.
        if "stages" not in workflow.output or len(workflow.output["stages"]) < 1:
            print "Wow, something is REALLY wrong with this job"
        elif "exception" in workflow.output["stages"][-1] and self.keep_workdir:
            print "Not cleaning up workdir so you can debug as necessary."
            pass
        elif workflow.workdir != ".":
            if subprocess.call(["rm", "-rf", workflow.workdir]) != 0:
                print "{0}:  warning, couldn't clean up {1}".format(
                       self.name, workflow.workdir)
        # Return results
        self.N_jobs -= 1
        return workflow


class LowMemQueue(Queue):
    """Parallel queue class for workflows including file management

    Similar class to Queue above, except it doesn't have standing worker
    threads, which may be subject to memory leaks.  Instead it simply creates
    one process per Workflow.
    """

    def __init__(self, name = "Generic Pipeline LowMemQueue",
                 processors = min(16,processing.cpuCount()),
                 keep_workdir=False):
        """Initializes the queue, but doesn't start any worker threads"""
        processing.freezeSupport()
        self.name = name
        self.N_processors = processors
        self.N_jobs = 0
        # Sorting out names:  These Queues are processing.Queues...
        self.Qin = processing.Queue(self.N_processors)
        self.Qout = processing.Queue(self.N_processors)
        self.workers = [ ]
        self.keep_workdir = keep_workdir

    @staticmethod
    def _worker(input, output):
        """A wrapper allowing Pipeline.Workflow instances to be queued"""
        # Just take one item from input and add one item to output.
        job, args, kwargs = input.get()
        output.put(job.run(*args,**kwargs))

    def start_job(self, workflow, *inputs, **kwinputs):
        s = super(LowMemQueue, self).start_job(workflow, *inputs, **kwinputs)
        if s == 0:
        # Now that that's done, just start the process
            self.workers.append(processing.Process
                    (target=LowMemQueue._worker, args=(self.Qin,self.Qout)))
            self.workers[-1].start()
        return s

    def finish_job(self):
        s = super(LowMemQueue, self).finish_job()
        if bool(s):
        # Now that that's done, terminate the process
            if not self.workers[0].isAlive():
                self.workers[0].join()
                del self.workers[0]
        return s
