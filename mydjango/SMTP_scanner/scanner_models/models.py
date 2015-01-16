from django.db import models

# Create your models here.
class Scanner_Date(models.Model):
    id_date = models.AutoField(primary_key=True)
    date = models.DateField()
    begin_mjd = models.FloatField()
    end_mjd = models.FloatField()
    photometric = models.BooleanField()
    n_exposure = models.IntegerField(default=0)
    n_subtraction = models.IntegerField(default=0)
    n_candidate = models.IntegerField(default=0)
    n_new_candidate = models.IntegerField(default=0)
    n_known_sn = models.IntegerField(default=0)
    n_sn_nova = models.IntegerField(default=0)
    n_varstar = models.IntegerField(default=0)
    n_AGN = models.IntegerField(default=0)
    n_astroid = models.IntegerField(default=0)
    n_other = models.IntegerField(default=0)
    n_scanned = models.IntegerField(default=0)
    
    def __unicode__(self):
        return '%04d-%02d-%02d'%(self.date.year,self.date.month,self.date.day)
    class Meta:
        ordering=['-date']

class Scanner_Candidate(models.Model):
    id_candidate = models.AutoField(primary_key=True)
    name = models.CharField(max_length=20)
    score = models.IntegerField(blank=True, null=True)
    ra = models.FloatField(blank=True, null=True)
    dec = models.FloatField(blank=True, null=True)
    first_mjd = models.FloatField(blank=True, null=True)
    last_mjd = models.FloatField(blank=True, null=True)
    n_obs_night = models.IntegerField(blank=True, null=True)
    type_best = models.CharField(max_length=10)
    type_best_comment = models.CharField(max_length=20)
    type_auto = models.CharField(max_length=10)
    type_auto_comment = models.CharField(max_length=20)
    fup_status = models.BooleanField(default=0)
    fup_type = models.CharField(max_length=10)
    def __unicode__(self):
        return self.name
    class Meta:
        ordering=['id_candidate']

class Scanner_ObjectType(models.Model):
    id_type = models.AutoField(primary_key=True)
    type = models.CharField(max_length=15)
    def __unicode__(self):
        return self.type
    class Meta:
        ordering=['id_type']

    
