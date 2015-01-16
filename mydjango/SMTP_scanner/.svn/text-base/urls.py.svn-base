from django.conf.urls.defaults import patterns, include, url
from views import SMTP_scanner_home,SMTP_scanner_date,SMTP_scanner_candidate,SMTP_scanner_register
from django.contrib.auth.views import login, logout
import settings

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
                       (r'^scanner/$',SMTP_scanner_home),
                       (r'^scanner/(\d{4})/(\d+)/(\d+)/$',SMTP_scanner_date),
                       (r'^scanner/cand_(\S+)/$',SMTP_scanner_candidate),
                       (r'^accounts/register/$', SMTP_scanner_register),
                       (r'^accounts/login/$',  login, {'template_name':'scanner_login.html'}),
                       (r'^accounts/logout/$', logout, {'template_name':'scanner_logout.html'}),
                       (r'^scanner/media/(?P<path>.*)$','django.views.static.serve',{'document_root': settings.MEDIA_ROOT}),

    # Examples:
    # url(r'^$', 'SMTP_scanner.views.home', name='home'),
    # url(r'^SMTP_scanner/', include('SMTP_scanner.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
                       url(r'^admin/', include(admin.site.urls)),
)
