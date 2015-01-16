from django.db import models

# Create your models here.


# ----------------------------------------------------------------------------
#                SQL models for pipeline job status logging
# ----------------------------------------------------------------------------


class ObsNight(models.Model):
    """A night's worth of observations on a given UT date."""
    ut_date = models.DateField(primary_key=True)

    def __unicode__(self):
        return '{0:04d}-{1:02d}-{2:02d}'.format(
                self.ut_date.year, self.ut_date.month, self.ut_date.day)
    class Meta:
        ordering=['-ut_date']


class SkymapperField(models.Model):
    """One of the 3900 SkyMapper fields.
    
    (ra,dec) = the field center in decimal degrees.
    """
    id  = models.IntegerField(primary_key=True)
    ra  = models.FloatField()
    dec = models.FloatField()

    def __unicode__(self):
        return "{0} [{1:9.5f},{2:9.5f}]".format(self.id, self.ra, self.dec)
    class Meta:
        ordering=['id']


class SkymapperPointing(models.Model):
    """An individual pointing at a particular field."""
    id      = models.AutoField(primary_key=True)
    night   = models.ForeignKey(ObsNight)
    filter  = models.CharField(max_length=1)
    dateobs = models.DateTimeField()
    filename = models.CharField(max_length=60)
    exptime = models.FloatField(default=0.0)
    obstype = models.CharField(max_length=8)
    status  = models.IntegerField(default=0)

    def __unicode__(self):
        return "{0:.1f}-sec exposure in filter {1} at {2}".format(
                self.exptime, self.filter, self.dateobs)
    class Meta:
        ordering=['dateobs']


class SciencePointing(SkymapperPointing):
    """Specifically for science pointings, not biases or flats"""
    field   = models.ForeignKey(SkymapperField)
    airmass  = models.FloatField(default=1.0)
    seeing   = models.FloatField(default=99.9)
    seewid   = models.FloatField(default=1.0)
    elong    = models.FloatField(default=99.9)
    elongwid = models.FloatField(default=1.0)
    zp       = models.FloatField(default=0.0)

    def __unicode__(self):
        return "{0:.1f}-sec exposure of field {1} in filter {2} at {3}".format(
            self.exptime, self.field_id, self.filter, self.dateobs)


class PipelineRun(models.Model):
    """A particular run of the pipeline.
    
    This will correspond either to a production run or a test run,
    with a particular pipeline version as encoded in Utils.Constants.
    """
    id          = models.AutoField(primary_key=True)
    runtag      = models.CharField(max_length=30)
    svn_version = models.CharField(max_length=16)
    start_time  = models.DateTimeField()


class PipelineExitStatus(models.Model):
    """Exit status from a particular pipeline stage.
    
    These will be created as needed, as new failure modes are discovered.
    """
    id          = models.AutoField(primary_key=True)
    stage_fail  = models.CharField(max_length=30)
    description = models.CharField(max_length=255)

    def __unicode__(self):
        return self.description
    class Meta:
        ordering=['id']


class PipelineJob(models.Model):
    """A single job run by the pipeline, e.g., a subtraction.

    RS 2012/07/06:  Because we have a proliferating number of "mini-pipelines"
    designed to carry out different tasks, PipelineJob now tracks any job from
    one of these mini-pipelines.  They are distinguished by different prefixes
    to the job name:  ingest ("ing"), superflat building ("cal"), reference
    cache maintenance ("ref"), and subtractions ("sub").

    Example name:  sub20120218_025714_3279_g_2011-08-25T23:01:01_17
    
    The first part of the name will always be the prefix followed by a runtag
    describing when the run was started.  There will usually also be a time
    stamp of the form "YYYY-MM-DDThh:mm:ss" corresponding to the "obsseq",
    or the date and time at which the Ob: request for the primary image was
    sent to TAROS by the scheduler.
    """
    jobname     = models.CharField(primary_key=True,max_length=60)
    start_time  = models.DateTimeField()
    end_time    = models.DateTimeField()
    usr_time    = models.FloatField()
    cpu_time    = models.FloatField()
    walltime    = models.FloatField()
    ccd         = models.IntegerField()
    pointing    = models.ForeignKey(SkymapperPointing)
    run         = models.ForeignKey(PipelineRun)
    exit_status = models.ForeignKey(PipelineExitStatus)

    def __unicode__(self):
        return self.jobname
    class Meta:
        ordering=['jobname']


class SubtractionJob(PipelineJob):
    """A PipelineJob carrying extra indices + QA info for subtractions.

    RS 2013/02/21:  I pulled this class off to the side in order to be able
    to quickly sort and search on QA information relevant to subtractions.
    """
    refpointing = models.ForeignKey(SciencePointing)
    minra = models.FloatField(default=0.0)
    maxra = models.FloatField(default=0.0)
    mindec = models.FloatField(default=0.0)
    maxdec = models.FloatField(default=0.0)
