from django.conf.urls.defaults import patterns, include, url
from views import SMTP_rb_trainer, SMTP_rb_trainer_date, SMTP_trainer_register
from django.contrib.auth.views import login, logout
import settings

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
                       (r'^trainer/$',SMTP_rb_trainer),
                       (r'^trainer/(\d{4})/(\d+)/(\d+)/$',SMTP_rb_trainer_date),
                       (r'^accounts/register/$', SMTP_trainer_register),
                       (r'^accounts/login/$',  login, {'template_name':'trainer_login.html'}),
                       (r'^accounts/logout/$', logout, {'template_name':'trainer_logout.html'}),
                       (r'^thumbs/(?P<path>.*)$','django.views.static.serve',{'document_root': settings.MEDIA_ROOT}),
                       
    # Examples:
    # url(r'^$', 'SMTP_trainer.views.home', name='home'),
    # url(r'^SMTP_trainer/', include('SMTP_trainer.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
                           url(r'^admin/', include(admin.site.urls)),
)
