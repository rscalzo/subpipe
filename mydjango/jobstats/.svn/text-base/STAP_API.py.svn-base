#!/usr/bin/env python

# ============================================================================
# RS 2012/05/18:  Django API for pipeline
# ----------------------------------------------------------------------------
# This allows the pipeline to interact with jobstats.models.
# ============================================================================

import datetime
import numpy as np
from Utils import Constants
import mydjango.jobstats.models as js


def SMfield_lookup(ra, dec):
    """Find SkyMapper field corresponding to (ra, dec) in decimal degrees."""
    
    # There's no nice way to do this until and unless the scheduler starts
    # including field information.  We therefore need to search the entire
    # database for fields.
    # if SMfield_list == None:
    SMfield_list = list(js.SkymapperField.objects.all())
    # Calculate the arc distance to all field centers; return the closest.
    # TODO:  Make this into a proper arc distance
    arcdist = np.sqrt([(ra-f.ra)**2 + (dec-f.dec)**2 for f in SMfield_list])
    return SMfield_list[arcdist.argmin()]

def register_pointing(names):
    """Registers a pointing in Django.  Called by STAP_ingest.py.
    
    names:  A Constants.Filenames instance."""

    # Since the Filenames structure from Constants.py already accessed these
    # (immutable) attributes via the header, just read them from there.
    night = js.ObsNight.objects.get_or_create(ut_date=names.dateobs.date())[0]

    if names.obstype == 'object': # hasattr(names, 'field'):
        field = js.SkymapperField.objects.get(id=names.field)
        pointing = js.SciencePointing.objects.get_or_create(
                filename=Constants.get_basename(names.basefname),
                night=night, field=field, filter=names.filter,
                obstype=names.obstype, dateobs=names.dateobs,
                exptime=names.exptime)[0]
        pointing.airmass=names.airmass
        pointing.seeing=names.seeing
        pointing.seewid=names.seewid
        pointing.elong=names.elong
        pointing.elongwid=names.elongwid
        pointing.zp=names.zpmag
        pointing.save()
        return pointing
    else:
        return js.SkymapperPointing.objects.get_or_create(
                filename=Constants.get_basename(names.basefname),
                night=night, filter=names.filter, obstype=names.obstype,
                dateobs=names.dateobs, exptime=names.exptime)[0]

def register_run(start_run_time):
    """Registers a pipeline run in Django.  Called by subpipe_master.py.
    
    start_run_time:  A timestamp as returned by time.time().
    """

    runtag = Constants.get_runtag(start_run_time)
    start_run_datetime = datetime.datetime.fromtimestamp(start_run_time)
    run = js.PipelineRun.objects.get_or_create(start_time=start_run_datetime,
            svn_version=Constants.svn_version(), runtag=runtag)

def register_job(meta, output, refmeta=None):
    """Registers a subtraction job in Django.  Called by subpipe_master.py.

    RS 2012/07/06:  Because we have a proliferating number of "mini-pipelines"
    designed to carry out different tasks, PipelineJob now tracks any job from
    one of these mini-pipelines.  They are distinguished by different prefixes
    to the job name:  ingest ("ing"), superflat building ("cal"), reference
    cache maintenance ("ref"), and subtractions ("sub").

    meta:  A Constants.Filenames instance with relevant job metadata.
    output:  The outputs structure from the corresponding Workflow.
    """

    # Total up the time from all the pipeline job stages
    # We could track time from each stage separately, but that starts making
    # the database more complicated than I want it to be.  If some jobs are
    # slower than others, we can figure out why by looking at the log files.
    stages = output['stages']
    totaltimes = { }
    for label in ('usr_time', 'cpu_time', 'walltime'):
        totaltimes[label] = sum([stages[i][label] for i in range(len(stages))])

    # Find the corresponding pointing, run, and exit status if they exist.
    pointing = register_pointing(meta)
    run = js.PipelineRun.objects.get(runtag=meta.runtag)
    if "exception" in stages[-1]:
        exit_status = js.PipelineExitStatus.objects.get_or_create(
                description=stages[-1]["exception"],
                stage_fail=stages[-1]["name"])[0]
    else:
        exit_status = js.PipelineExitStatus.objects.get_or_create(
            description="Success")[0]

    # Register the job in the database
    if refmeta is None:
        # Not a subtraction job
        job = js.PipelineJob.objects.get_or_create(
                  jobname=output['name'], ccd=meta.ccd,
                  start_time=stages[0]['t_start'],
                  end_time=stages[-1]['t_finish'],
                  usr_time=totaltimes['usr_time'],
                  cpu_time=totaltimes['cpu_time'],
                  walltime=totaltimes['walltime'],
                  pointing=pointing, run=run, exit_status=exit_status)
    else:
        # Subtraction job; add REF pointing and other fun things
        refpointing = register_pointing(refmeta)
        minra, mindec = 0.0, 0.0
        maxra, maxdec = 0.0, 0.0
        for s in stages:
            if s is None:  continue
            try:
                if 'results' in s and 'minra' in s['results']:
                    minra, mindec = s['results']['minra'], s['results']['mindec']
                    maxra, maxdec = s['results']['maxra'], s['results']['maxdec']
            except:
                continue
        job = js.SubtractionJob.objects.get_or_create(
                  jobname=output['name'], ccd=meta.ccd,
                  start_time=stages[0]['t_start'],
                  end_time=stages[-1]['t_finish'],
                  usr_time=totaltimes['usr_time'],
                  cpu_time=totaltimes['cpu_time'],
                  walltime=totaltimes['walltime'],
                  pointing=pointing,
                  refpointing=refpointing,
                  minra=minra, mindec=mindec,
                  maxra=maxra, maxdec=maxdec,
                  run=run, exit_status=exit_status)
    return job
