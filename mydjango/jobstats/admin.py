from mydjango.jobstats.models import ObsNight, \
        SkymapperField, SkymapperPointing, SciencePointing, \
        PipelineRun, PipelineExitStatus, PipelineJob, SubtractionJob, SkymapperImage
from django.contrib import admin

class ObsNightAdmin(admin.ModelAdmin):
    list_display = ['ut_date']
    search_fields = ['ut_date']

class SkymapperFieldAdmin(admin.ModelAdmin):
    list_display = ['id', 'ra', 'dec','label','rank']
    search_fields = ['id','label','rank']

class SkymapperPointingAdmin(admin.ModelAdmin):
    list_display = ['filter', 'dateobs', 'filename', 'exptime']
    search_fields = ['filter', 'filename']

class SciencePointingAdmin(admin.ModelAdmin):
    list_display = ['field', 'filter', 'dateobs', 'filename', 'seeing', 'elong', 'zp']
    search_fields = ['filter', 'filename']

class SkymapperImageAdmin(admin.ModelAdmin):
    list_display = ['pointing', 'ccd', 'seeing', 'seewid','elong','elongwid','zp','bkgmed','bkgsig','status']
    search_fields = ['pointing', 'status']

class PipelineRunAdmin(admin.ModelAdmin):
    list_display = ['runtag', 'svn_version', 'start_time']
    search_fields = ['runtag']

class PipelineJobAdmin(admin.ModelAdmin):
    list_display = ['jobname', 'start_time', 'end_time', 'walltime',
                    'pointing', 'ccd', 'exit_status']
    search_fields = ['jobname']

class SubtractionJobAdmin(admin.ModelAdmin):
    list_display = ['jobname', 'walltime', 'minra', 'maxra', 'mindec', 'maxdec']
    search_fields = ['jobname']

admin.site.register(ObsNight, ObsNightAdmin)
admin.site.register(SkymapperField, SkymapperFieldAdmin)
admin.site.register(SkymapperPointing, SkymapperPointingAdmin)
admin.site.register(SciencePointing, SciencePointingAdmin)
admin.site.register(SkymapperImage,SkymapperImageAdmin)
admin.site.register(PipelineRun, PipelineRunAdmin)
admin.site.register(PipelineJob, PipelineJobAdmin)
admin.site.register(SubtractionJob, SubtractionJobAdmin)
admin.site.register(PipelineExitStatus)
