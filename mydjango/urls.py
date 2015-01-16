from django.conf.urls.defaults import patterns, include, url
from django.contrib.auth.views import login, logout
import settings

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('mydjango.views',
                       # Examples:
                       # url(r'^$', 'mydjango.views.home', name='home'),
                       # url(r'^mydjango/', include('mydjango.foo.urls')),
                       #
                       #home, in .views
                       url(r'^$','smt_home',name='home'),
                       url(r'^/$','smt_home', name='home2'),
                       url(r'^external/pessto/','page_pessto',name='page_pessto'),
                       #url(r'^accounts/register/$', 'smt_register',name='smt_register'),
                       )

#jobstats, in .jobstats.views
urlpatterns += patterns('mydjango.jobstats.views',                       
                        url(r'^jobstats/$', 'index', name='jobstats_home'),
                        url(r'^jobstats/list_jobs/run=(\d{8}_\d{6})/$', 'listjobsbyrun',name='joblist'),
                        url(r'^jobstats/sublogfile/(\w{3}\d{8}_\d{6}.*)/$', 'displaylog',name='joblogs'),
                        # url(r'^jobstats/coveragemap/$'),
                        # url(r'^jobstats/cadencemap/$'),
                        # url(r'^jobstats/controlmap/$'),
                        )

#accounts stuff
urlpatterns += patterns('',
                        url(r'^accounts/login/$',login,{'template_name':'smt_login.html'}),
                        url(r'^accounts/logout/$', logout, {'template_name':'smt_logout.html'}),
                        url(r'^accounts/password/reset/$', 
                            'django.contrib.auth.views.password_reset', 
                            {'template_name':'smt_pw_reset_form.html',
                             'email_template_name':'smt_pw_reset_email.html',
                             'post_reset_redirect' : 'http://www.mso.anu.edu.au/skymapper/smt/accounts/password/reset/done/'},
                            name="password_reset"),
                        url(r'^accounts/password/reset/done/$',
                            'django.contrib.auth.views.password_reset_done',
                            {'template_name':'smt_pw_reset_done.html'},
                            name='password_reset_done'),
                        url(r'^accounts/password/reset/(?P<uidb36>[0-9A-Za-z]+)-(?P<token>.+)/$', 
                            'django.contrib.auth.views.password_reset_confirm', 
                            {'template_name':'smt_pw_reset_confirm.html',
                             'post_reset_redirect' : 'http://www.mso.anu.edu.au/skymapper/smt/accounts/password/done/'}),
                        url(r'^accounts/password/done/$', 
                            'django.contrib.auth.views.password_reset_complete',{'template_name':'smt_pw_done.html'}),
                        )                       
#media stuff
urlpatterns += patterns('',
                        url(r'^media/(?P<path>.*)$','django.views.static.serve',{'document_root': settings.MEDIA_ROOT}),
                        )                        
#candidate pages
urlpatterns += patterns('mydjango.followup.views',
                        url(r'^transients/all/$','transients_all',name='trans_all'),
                        url(r'^transients/latest/$','transients_date_latest',name='trans_latest'),
                        url(r'^transients/(\d{4})/(\d+)/(\d+)/$','transients_date',name='trans_date'),
                        url(r'^transients/(SKY_J\d+[+-]\d+)/$','transient_info',name='trans_info'),           
                        url(r'^transients/(SMTJ\d+[+-]\d+)/$','transient_info',name='trans_info'),
                        )

urlpatterns += patterns('',
                        # Uncomment the admin/doc line below to enable admin documentation:
                        # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),
                        
                        # Uncomment the next line to enable the admin:
                        url(r'^admin/', include(admin.site.urls)),
                        )

