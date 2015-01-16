from mydjango.followup.models import Transient, Transient_Change, \
        TransientType
from django.contrib import admin

class TransientAdmin(admin.ModelAdmin):
    list_display = ['name', 'ra', 'dec', 'field', 'ccd', 'type']
    search_fields = ['name']
    list_filter = ['type']

class TransientChangeAdmin(admin.ModelAdmin):
    list_display = ['transient', 'time', 'user', 'type', 'comment']
    search_fields = ['transient']
    list_filter = ['user', 'type']

class TransientTypeAdmin(admin.ModelAdmin):
    list_display = ['type', 'supertype']

admin.site.register(Transient, TransientAdmin)
admin.site.register(TransientType, TransientTypeAdmin)
admin.site.register(Transient_Change, TransientChangeAdmin)
