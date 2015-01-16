from django.db import models
from mydjango.jobstats.models import *
from django.contrib.auth.models import User
import datetime

# Create your models here.

# ----------------------------------------------------------------------------
#              SQL models for the follow-up website interface
# ----------------------------------------------------------------------------


class TransientType(models.Model):
    """Transient types.  Always a slippery thing."""
    id        = models.AutoField(primary_key=True)
    type      = models.CharField(max_length=15)
    supertype = models.ForeignKey('self', blank=True, null=True)

    def __unicode__(self):
        return self.type
    class Meta:
        ordering=['id']


class FollowUpStatus(models.Model):
    id     = models.AutoField(primary_key=True)
    status = models.CharField(max_length=15)

    def __unicode__(self):
        return self.status
    class Meta:
        ordering=['id']


class Transient(models.Model):
    """Known transients found by SkyMapper.
    
    Currently defined as anything we want to keep track of
    (~ things with > N coincident detections).
    """
    id                = models.AutoField(primary_key=True)
    name              = models.CharField(max_length=30)
    ra                = models.FloatField()
    dec               = models.FloatField()
    field             = models.ForeignKey(SkymapperField)
    ccd               = models.IntegerField()
    n_det             = models.IntegerField(default=0)
    n_obs             = models.IntegerField(default=0)
    pointing_det      = models.ManyToManyField(SciencePointing)
    pointing_obs      = models.ManyToManyField(SciencePointing,related_name="+")
    type              = models.ForeignKey(TransientType)
    follow_up_status  = models.ForeignKey(FollowUpStatus)

    def __unicode__(self):
        return self.name
    class Meta:
        ordering=['id']
    
    def save(self, *args, **kwargs):
        if self.id:
            self.n_det = self.pointing_det.count()
            self.n_obs = self.pointing_obs.count()
        super(Transient, self).save(*args, **kwargs)

    def savelog(self, *args, **kwargs):
        """Log change to change table"""
        if self.id:
            self.n_det = self.pointing_det.count()
            self.n_obs = self.pointing_obs.count()
            
        super(Transient, self).save() # Call the "real" save() method.

        if 'user' in kwargs:
            user=kwargs['user']
        else:
            username='unknown'
            (user,isnew)=User.objects.get_or_create(username=username)
        if 'comment' in kwargs:
            comment=kwargs['comment']
        else:
            comment=''
        change = Transient_Change.objects.create(transient=self,time=datetime.datetime.now(),
                                           user=user,comment=comment,
                                           type=self.type,follow_up_status=self.follow_up_status)


class FollowUpLog(models.Model):
    id        = models.AutoField(primary_key=True)
    transient = models.ForeignKey(Transient)
    time      = models.DateTimeField()
    jd        = models.FloatField()
    instrument= models.CharField(max_length=20)
    isspec    = models.BooleanField(default=0)
    isphoto   = models.BooleanField(default=0)
    filter    = models.CharField(max_length=10)

    def __unicode__(self):
        return self.jd
    class Meta:
        ordering=['id']

    
class Transient_Change(models.Model):
    id        = models.AutoField(primary_key=True)
    transient = models.ForeignKey(Transient)
    time      = models.DateTimeField()
    user      = models.ForeignKey(User)
    type      = models.ForeignKey(TransientType)
    follow_up_status  = models.ForeignKey(FollowUpStatus)
    comment   = models.CharField(max_length=255)

    def __unicode__(self):
        return "{0} changed at {1}".format(self.transient.name,self.time)
    class Meta:
        ordering=['id']
    

class CommentLog(models.Model):
    id        = models.AutoField(primary_key=True)
    transient = models.ForeignKey(Transient)
    jd        = models.FloatField()
    user      = models.ForeignKey(User)
    type      = models.ForeignKey(TransientType)
    comment   = models.CharField(max_length=255)

    def __unicode__(self):
        return self.jd
    class Meta:
        ordering=['id']
