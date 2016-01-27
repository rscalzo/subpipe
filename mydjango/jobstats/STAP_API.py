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


def register_pointing(names,update=False,addimage=False,proc_status=None):
    """Registers a pointing in Django.  Called by STAP_ingest.py.
    
    names:  A Constants.Filenames instance."""

    # Since the Filenames structure from Constants.py already accessed these
    # (immutable) attributes via the header, just read them from there.
    #night = js.ObsNight.objects.get_or_create(ut_date=names.dateobs.date())[0]

    if names.obstype == 'object': # hasattr(names, 'field'):
        night = js.ObsNight.objects.get_or_create(ut_date=names.dateobs.date())[0]
        field = js.SkymapperField.objects.get(id=names.field)
        pointing,isnew = js.SciencePointing.objects.get_or_create(
                filename=Constants.get_basename(names.basefname),
                night=night, field=field, filter=names.filter,
                obstype=names.obstype, dateobs=names.dateobs,
                exptime=names.exptime)
        if proc_status:
                pointing.proc_status=proc_status
                pointing.save()
        #find out what kind of observation it is from filename
        #only update if new
        if isnew or update:
            #get program, surveyid, sequence
            try:
                obs_id=Constants.decode_obs_id(pointing.filename.split('_')[1])
                pointing.program=obs_id[0]
                pointing.surveyid=obs_id[1]
                pointing.sequence=obs_id[2]
            except:
                if hasattr(names, 'program'):
                    pointing.program=names.program
            pointing.airmass=names.airmass
            pointing.seeing=names.seeing
            pointing.seewid=names.seewid
            pointing.elong=names.elong
            pointing.elongwid=names.elongwid
            pointing.zp=names.zpmag
            pointing.bkgsig=names.bkgsig
            pointing.bkgmed=names.bkgmed
            pointing.maglim50=names.maglim50
            pointing.maglim95=names.maglim95
            pointing.save()
        if isnew:
            #update field lastobs
            if pointing.program=='3PS':
                field.lastobs=field.sciencepointing_set.filter(program='3PS').order_by('-dateobs')[0].dateobs
                field.save()
            #also register 32 CCD images
            if addimage:
                ccds=range(1,33)
                for ccd in ccds:
                    image,isnew=js.SkymapperImage.objects.get_or_create(pointing=pointing,ccd=ccd)
                    if isnew:
                        image.seeing=names.seeing
                        image.seewid=names.seewid
                        image.elong=names.elong
                        image.elongwid=names.elongwid
                        image.zp=names.zpmag
                        image.bkgmed=names.bkgmed
                        image.bkgsig=names.bkgsig
                        image.maglim50=names.maglim50
                        image.maglim95=names.maglim95
                        image.save()

        return pointing
    else:
        #FY - 2015-08-20, there was an error, so at this point, don't bother with nonscience frames
        return None
        return js.SkymapperPointing.objects.get_or_create(
                filename=Constants.get_basename(names.basefname),
                night=night, filter=names.filter, obstype=names.obstype,
                dateobs=names.dateobs, exptime=names.exptime)[0]
    
def register_image(meta,pointing=None,update=True,isref=False,updaterefim=False,updatepointing=False,refondisk=1,refcheck=True):
    if meta.subfield==0:
        return None
    if pointing is None:
        pointing=register_pointing(meta,addimage=False)
    if isref:
        refim,isnew=js.SkymapperImage.objects.get_or_create(pointing=pointing,ccd=meta.subfield)
        #image,newref=js.SubtractionRef.objects.get_or_create(pointing=pointing,ccd=meta.subfield,skymapperimage__ptr=refim)
        image,newref=js.SubtractionRef.objects.get_or_create(image=refim)
        if refcheck:
            if refim.seeing > Constants.MAX_REFSEEING or refim.elong > Constants.MAX_REFELONG or refim.maglim95 < Constants.MIN_REFDEPTH[refim.pointing.filter]:
                refim.status= Constants.FLAG_BADREF
                refim.save()
        if updaterefim:
            if refim.status==0:
                refim.status=1
                refim.save()
        if newref:
            image.ondisk=refondisk
            image.save()
    else:
        image,isnew=js.SkymapperImage.objects.get_or_create(pointing=pointing,ccd=meta.subfield)
    if isnew or update:
        image.seeing=meta.seeing
        image.seewid=meta.seewid
        image.elong=meta.elong
        image.elongwid=meta.elongwid
        image.zp=meta.zpmag
        image.bkgmed=meta.bkgmed
        image.bkgsig=meta.bkgsig
        image.maglim50=meta.maglim50
        image.maglim95=meta.maglim95
        image.minra=meta.minra
        image.maxra=meta.maxra
        image.mindec=meta.mindec
        image.maxdec=meta.maxdec       
        image.save()
    if updatepointing:
        pointing.seeing=image.seeing
        pointing.seewid=image.seewid
        pointing.elong=image.elong
        pointing.elongwid=image.elongwid
        pointing.zp=image.zp
        pointing.bkgsig=image.bkgsig
        pointing.bkgmed=image.bkgmed
        pointing.maglim50=image.maglim50
        pointing.maglim95=image.maglim95
        pointing.save()
    return image


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
    pointing = register_pointing(meta,addimage=False)
    image = register_image(meta,pointing=pointing)
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
            jobname=output['name'], ccd=meta.subfield,
            start_time=stages[0]['t_start'],
            end_time=stages[-1]['t_finish'],
            usr_time=totaltimes['usr_time'],
            cpu_time=totaltimes['cpu_time'],
            walltime=totaltimes['walltime'],
            pointing=pointing, image=image,
            run=run, exit_status=exit_status)
        image.status=1
        if "exception" in stages[-1]:
            image.status+=10
        else:
            refimage=register_image(meta,pointing=pointing,update=True,isref=True)
        image.save()
    else:
        # Subtraction job; add REF pointing and other fun things
        refpointing = register_pointing(refmeta,addimage=False)
        refimage=register_image(refmeta,pointing=refpointing,update=False,isref=True)
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
            jobname=output['name'], ccd=meta.subfield,
            start_time=stages[0]['t_start'],
            end_time=stages[-1]['t_finish'],
            usr_time=totaltimes['usr_time'],
            cpu_time=totaltimes['cpu_time'],
            walltime=totaltimes['walltime'],
            pointing=pointing,
            image=image,
            refpointing=refpointing,
            ref=refimage,
            minra=minra, mindec=mindec,
            maxra=maxra, maxdec=maxdec,
            run=run, exit_status=exit_status)
        image.status=2
        if "exception" in stages[-1]:
            image.status+=10
        image.save()
    return job
