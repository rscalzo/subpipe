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
from django.utils import simplejson
import ephem
from mydjango.settings import MEDIA_URL, MEDIA_ROOT
from Utils.Constants import decode_obs_id
#views for main smt pages

def smt_home(request):
    template="smt_home.html"
    
    #grab the run tags
    smt_runs=js.PipelineRun.objects.values('runtag').order_by('-start_time')

    #get all the ObsNight
    smt_nights=js.ObsNight.objects.values('ut_date').order_by('-ut_date')
    #get all the transient names
    smt_trans=fu.Transient.objects.filter(type__type__in=['Cand','SN','Ia','II','Ibc','IIn','SLSN','FUP']).values('name').order_by('-id')
    #smt_trans=fu.Transient.objects.values('name').order_by('-id')

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
        last_obsnight=smt_nights[0]['ut_date']
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


#transient list for public page
from bin.transients_for_public import PublicTarget
        
#public page
#no login required
def smt_public(request):
    template='smt_public.html'
    publicfile='public.csv'
    publicurl=MEDIA_URL+'/'+publicfile
    publiclocal=MEDIA_ROOT+'/'+publicfile
    #parse the public table
    transients=PublicTarget.read_ascii_file(publiclocal)
    
    headers=[]
    for f in PublicTarget._ascii_fields:
        headers.append(f[-1])

    rows=[]
    for trans in transients:
        rowhtml=''
        rowdata=[]
        for f in PublicTarget._ascii_fields:
            if 'URL' in f[0]:
                val='<img src=\'%s\' alt=\'%s\' width=150px/>'%(eval('trans.'+f[0]),f[-1])
            else:
                val=eval('trans.'+f[0])
            rowdata.append(val)
        rows.append({'html':rowhtml,'data':rowdata})

    return render_to_response(template,{'publicurl':publicurl,'headers':headers,'rows':rows},context_instance=RequestContext(request))


@login_required
def smt_map(request):
    """plot skymap of observations"""
    
    template='smt_map.html'
    
    #default values
    progtypes=['BAD','3PS']
    
    filt='g'
    
    today=datetime.date.today()
    fromnight=today+datetime.timedelta(-1)
    tonight=today+datetime.timedelta(-1)
    
    #or submitted value
    if request.method=="POST":
        requestkeys=request.POST.keys()
        if 'plot_bad' in requestkeys:
            if not 'BAD' in progtypes: 
                progtypes.append('BAD')
        else:
            if 'BAD' in progtypes: 
                progtypes.remove('BAD')
        if 'plot_3ps' in requestkeys:
            if not '3PS' in progtypes: 
                progtypes.append('3PS')
        else:
            if '3PS' in progtypes: 
                progtypes.remove('3PS')
        if 'plot_ms' in requestkeys:
            if not 'MS' in progtypes: 
                progtypes.append('MS')
        else:
            if 'MS' in progtypes: 
                progtypes.remove('MS')
        if 'plot_filt' in requestkeys:
            filt=request.POST.get('plot_filt')
        if 'from' in requestkeys:
            fromnightstr=request.POST.get('from')
            try:
                fromnight = datetime.datetime.strptime(fromnightstr, '%Y-%m-%d').date()
            except:
                pass
        if 'to' in requestkeys:
            tonightstr=request.POST.get('to')
            try:
                tonight = datetime.datetime.strptime(tonightstr, '%Y-%m-%d').date()
            except:
                pass


     #select observations
    sp=js.SciencePointing.objects.filter(filter=filt,night__ut_date__gte=fromnight).filter(night__ut_date__lte=tonight).exclude(field__id=0)
    temp=len(sp)

    sp_dict={}
    sp_dict['fields']=[]
    for sp_obj in sp:
         #MS, 3PS or BAD
        progtype=decode_obs_id(sp_obj.filename.split('_')[1])[0]
        if not progtype in progtypes:
            continue
        sp_dict['fields'].append({'ra_deg':sp_obj.field.ra,'dec_deg':sp_obj.field.dec,'status':progtype.lower()+'_ingested','selected':0})

    sp_json = simplejson.dumps(sp_dict,skipkeys=True,sort_keys=True)
    return render_to_response(template,
                              {'progtypes':progtypes,
                               'fromnight':fromnight,'tonight':tonight,'filter':filt,
                               'fields':sp_json},
                              context_instance=RequestContext(request))


@login_required
def smt_cadence(request):
    """Plot cadence of last 100 days.
    color sequence for each field
    cadence histrogram
    """
    
    template='smt_cadence.html'
    
    #default values
    today=datetime.date.today()
    fromnight=datetime.date(2015,7,1)
    #fromnight=today+datetime.timedelta(-100)
    tonight=today

    if request.method=="POST":
        requestkeys=request.POST.keys()
        if 'from' in requestkeys:
            fromnightstr=request.POST.get('from')
            try:
                fromnight = datetime.datetime.strptime(fromnightstr, '%Y-%m-%d').date()
            except:
                pass
        if 'to' in requestkeys:
            tonightstr=request.POST.get('to')
            try:
                tonight = datetime.datetime.strptime(tonightstr, '%Y-%m-%d').date()
            except:
                pass

    #select observations                                                                    
    #exclude MS
    sp=js.SciencePointing.objects.filter(night__ut_date__gte=fromnight,night__ut_date__lte=tonight,program='3PS').exclude(field__id=0).select_related('field__id')
    temp=len(sp)
    
    fids=list(set(sp.values_list('field__id',flat=True)))

    #will translate fids to y-axis
    fids.sort()
    ticks=[]

    #fill cadence and field dates dict
    filters=['v','g','r','i']
    cadence={}
    fdates={}
    for filt in filters:
        cadence[filt]=[]
        fdates[filt]={'label':filt,'data':[]}
    
    filtoffset=[0.3,0.1,-0.1,-0.3]
    for fidnum,fid in enumerate(fids):
        ticks.append([fidnum,'%04d'%fid])
        for filtnum,filt in enumerate(filters):
            obs=sp.filter(field__id=fid,filter=filt)
            dates=list(obs.values_list('dateobs',flat=True))
            if len(dates)==0:
                continue
            
            dates.sort()
            mjds=[ephem.Date(date) + 15019.5 for date in dates]
            for mjd in mjds: fdates[filt]['data'].append([mjd,fidnum+filtoffset[filtnum]])
            
            if len(mjds)>1:
                cad=[mjds[did+1]-mjd for did,mjd in enumerate(mjds[:-1])]              
                for icad in cad: cadence[filt].append(icad)
                

    #make cadence histrogram
    hist={}
    nbin=10.
    for filt in filters:
        hist[filt]={'label':filt,'data':[]}
        if len(cadence[filt])==0:continue
        maxcad=round(max(cadence[filt])+0.5)
        mincad=round(min(cadence[filt])-0.5)
#        binsize=round((maxcad-mincad)/nbin+0.5)
        binsize=1.
        bins=range(int(mincad),int(maxcad+binsize),int(binsize))
        for bid,binl in enumerate(bins[:-1]):
            inbin=[icad for icad in cadence[filt] if icad>=binl and icad<bins[bid+1]]
            #hist[filt]['data'].append([(binl+bins[bid+1])/2.,len(inbin)])
            hist[filt]['data'].append([binl,len(inbin)])

    cad_json = simplejson.dumps(hist,skipkeys=True,sort_keys=True)
    fdates_json = simplejson.dumps(fdates,skipkeys=True,sort_keys=True)
    return render_to_response(template,
                              {'fromnight':fromnight,'tonight':tonight,
                               'cad_hist':cad_json,
                               'fdates':fdates_json,
                               'ticks':ticks},
                              context_instance=RequestContext(request))

