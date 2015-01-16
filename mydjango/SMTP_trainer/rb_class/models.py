from django.db import models

# Create your models here.

class Users(models.Model):
    name=models.CharField(max_length=20)
    userid=models.IntegerField()
    def __unicode__(self):
        return self.name
    class Meta:
        ordering=['userid']

class Candidate(models.Model):
    name=models.CharField(max_length=15)
    date=models.DateField(blank=True, null=True)
    thumb_new=models.ImageField(upload_to='/')
    thumb_ref=models.ImageField(upload_to='/')
    thumb_diff=models.ImageField(upload_to='/')
    # basic
    xsub=models.FloatField(blank=True, null=True)
    ysub=models.FloatField(blank=True, null=True)
    asub=models.FloatField(blank=True, null=True)
    bsub=models.FloatField(blank=True, null=True)
    esub=models.FloatField(blank=True, null=True)
    thsub=models.FloatField(blank=True, null=True)
    fwhmsub=models.FloatField(blank=True, null=True)
    f4sub=models.FloatField(blank=True, null=True)
    df4sub=models.FloatField(blank=True, null=True)
    f8sub=models.FloatField(blank=True, null=True)
    df8sub=models.FloatField(blank=True, null=True)
    flagsub=models.IntegerField(blank=True, null=True)
    starsub=models.FloatField(blank=True, null=True)
    # Quantities you need pixel data to calculate
    # Mostly counts of some property in a square box around the candidate
    n2sig3=models.IntegerField(blank=True, null=True)
    n3sig3=models.IntegerField(blank=True, null=True)
    n2sig5=models.IntegerField(blank=True, null=True)
    n3sig5=models.IntegerField(blank=True, null=True)
    #bad pixel mask
    nmask=models.IntegerField(blank=True, null=True)
    #global
    Rfwhm=models.FloatField(blank=True, null=True)
    goodcn=models.FloatField(blank=True, null=True)
    subconv=models.BooleanField(blank=True)
    #REF
    xref=models.FloatField(blank=True, null=True)
    yref=models.FloatField(blank=True, null=True)
    aref=models.FloatField(blank=True, null=True)
    bref=models.FloatField(blank=True, null=True)
    eref=models.FloatField(blank=True, null=True)
    thref=models.FloatField(blank=True, null=True)
    fwhmref=models.FloatField(blank=True, null=True)
    f4ref=models.FloatField(blank=True, null=True)
    df4ref=models.FloatField(blank=True, null=True)
    flagref=models.IntegerField(blank=True, null=True)
    starref=models.FloatField(blank=True, null=True)
    refsrc=models.BooleanField(blank=True)
    nndref=models.FloatField(blank=True, null=True)
    #NEW
    xnew=models.FloatField(blank=True, null=True)
    ynew=models.FloatField(blank=True, null=True)
    anew=models.FloatField(blank=True, null=True)
    bnew=models.FloatField(blank=True, null=True)
    enew=models.FloatField(blank=True, null=True)
    thnew=models.FloatField(blank=True, null=True)
    fwhmnew=models.FloatField(blank=True, null=True)
    f4new=models.FloatField(blank=True, null=True)
    df4new=models.FloatField(blank=True, null=True)
    flagnew=models.IntegerField(blank=True, null=True)
    starnew=models.FloatField(blank=True, null=True)
    newsrc=models.BooleanField(blank=True)
    nndnew=models.FloatField(blank=True, null=True)   
    # Other
    apsig4=models.FloatField(blank=True, null=True)
    apsig8=models.FloatField(blank=True, null=True)
    normrms=models.FloatField(blank=True, null=True)
    normfwhm=models.FloatField(blank=True, null=True)
    Ranew=models.FloatField(blank=True, null=True)
    Renew=models.FloatField(blank=True, null=True)
    Dthnew=models.FloatField(blank=True, null=True)
    Raref=models.FloatField(blank=True, null=True)
    Reref=models.FloatField(blank=True, null=True)
    Dthref=models.FloatField(blank=True, null=True)
    Rfref=models.FloatField(blank=True, null=True)
    Rfnew=models.FloatField(blank=True, null=True)
    
    # scores
    score1=models.FloatField(default=-1)
    score2=models.FloatField(default=-1)
    score3=models.FloatField(default=-1)
    score4=models.FloatField(default=-1)
    score5=models.FloatField(default=-1)
    score6=models.FloatField(default=-1)
    score7=models.FloatField(default=-1)
    score8=models.FloatField(default=-1)
    score9=models.FloatField(default=-1)
    score10=models.FloatField(default=-1)

    def __unicode__(self):
        return self.name
    class Meta:
        ordering=['name']

class ImageDate(models.Model):
    date=models.DateField()
    ncandidate=models.IntegerField(default=0)
    nscore1=models.IntegerField(default=0)
    nscore2=models.IntegerField(default=0)
    nscore3=models.IntegerField(default=0)
    nscore4=models.IntegerField(default=0)
    nscore5=models.IntegerField(default=0)
    nscore6=models.IntegerField(default=0)
    nscore7=models.IntegerField(default=0)
    nscore8=models.IntegerField(default=0)
    nscore9=models.IntegerField(default=0)
    nscore10=models.IntegerField(default=0)
    def __unicode__(self):
        return '%04d-%02d-%02d'%(self.date.year,self.date.month,self.date.day)
    class Meta:
        ordering=['-date']
                


