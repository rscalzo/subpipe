from django.http import HttpResponse,HttpResponseRedirect, Http404
from django.shortcuts import render_to_response
from django.template import RequestContext
#from django.contrib import auth
from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django.contrib.auth.forms import UserCreationForm

import datetime
import mydjango.jobstats.models as js
import mydjango.followup.models as fu
from django.core.urlresolvers import reverse

#views for main smt pages

def smt_home(request):
    template="smt_home.html"
    
    #grab the run tags
    smt_runs=js.PipelineRun.objects.values('runtag')

    #get all the ObsNight
    smt_nights=js.ObsNight.objects.values('ut_date')
    #get all the transient names
    smt_trans=fu.Transient.objects.values('name')
    
    #fill in some statistics
    n_nights=smt_nights.count()
    n_pointings=js.SciencePointing.objects.count()
    n_fields=js.SciencePointing.objects.values('field').distinct().count()
    n_jobs=js.PipelineJob.objects.count()
    n_trans=smt_trans.count()
    
    #area of sky covered
    #number of transients per type
    
    smt_stats={'n_nights':n_nights,'n_pointings':n_pointings,'n_fields':n_fields,'n_jobs':n_jobs,'n_trans':n_trans}
    
    #figure out last night of observation
    if n_nights>0:
        last_obsnight=smt_nights.order_by('-ut_date')[0]['ut_date']
    else:
        last_obsnight=None

    #jump to certain night/transient
    if request.method=="POST":
        requestkeys=request.POST.keys()
        for requestkey in requestkeys:
            #if requestkey in 'smt_logs_night':
            #    night=request.POST.get(requestkey)
            #    return HttpResponseRedirect(reverse('jobstats_date', args=tuple(night.split('-'))))
            if requestkey in 'smt_trans_night':
                night=request.POST.get(requestkey)
                return HttpResponseRedirect(reverse('trans_date',args=tuple(night.split('-'))))
            if requestkey in 'smt_transient':
                transname=request.POST.get(requestkey)
                return HttpResponseRedirect(reverse('trans_info',args=(transname,)))
            if requestkey in 'smt_run':
                runtag=request.POST.get(requestkey)
                return HttpResponseRedirect(reverse('joblist',args=(runtag,)))

    return render_to_response(template,{'smt_stats':smt_stats,'last_obsnight':last_obsnight,'smt_runs':smt_runs,
                                        'smt_nights':smt_nights,'smt_trans':smt_trans},context_instance=RequestContext(request))


def smt_register(request):
    if request.method=="POST":
        #add user account
        form = UserCreationForm(request.POST)
        if form.is_valid():
            form.save()
            
            username=form.cleaned_data['username']
            password=form.cleaned_data['password1']
            user=authenticate(username=username,password=password)
            
            if user is not None:
                if user.is_active:
                    login(request, user)
                    message = "Congratulations! You have been registered. <a href=\"%s\">Go back to SMT home ?</a>"%(reverse('home'))
                else:
                    message = "urh... you account might be disabled."
            else:
                message = "There was an error automatically logging you in. Try <a href=\"%saccounts/login/\"> logging in </a> manually."%(reverse('home'))
                
            return render_to_response('smt_register_success.html',{
                    'username':username,
                    'message':message,},context_instance=RequestContext(request))
    else:        
        form=UserCreationForm()            
        
    template='smt_register.html'
    return render_to_response(template,{'form':form},context_instance=RequestContext(request))
    
 
@login_required
def page_pessto(request):
    return render_to_response('smt_pessto.html',context_instance=RequestContext(request))
