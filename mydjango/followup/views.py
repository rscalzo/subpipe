from django.http import HttpResponse,HttpResponseRedirect, Http404
from django.shortcuts import render_to_response
from django.template import RequestContext
#from django.contrib import auth
from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger
from django.core.exceptions import ObjectDoesNotExist
from django.utils import simplejson
from mydjango.pysrc.SortHeaders import SortHeaders
from django.core.urlresolvers import reverse

from Utils.Constants import PipelinePath as CPP
from Utils.Constants import FilenamesSub, FilenamesXsient, ThumbNames
from Utils.Constants import follow_types
#from Utils.TargetVisuals import makefinder
import datetime
import pyfits
import numpy as np
import os
import re
import mydjango.jobstats.models as js
import mydjango.followup.models as fu

from django.db.models import Min, Max
import ephem

#views for transient pages

#def requires_login(view):
#    def new_view(request, *args, **kwargs):
#        if not request.user.is_authenticated():
#            return HttpResponseRedirect('/accounts/login/?next=%s'%request.path)
#        return view(request, *args, **kwargs)
#    return new_view

class trans_table():
    """create headers and rows to list transients from Transient table query"""
    def __init__(self,request,transients):
        list_headers=(
            ('Name','name'),
            ('Type','type'),
            ('RA','ra'),
            ('DEC','dec'),
            ('First Detect','first_detected'),
            ('Last Detect','last_detected'),
            ('Last Observed','last_observed'),
            ('# of Detections','n_det'),
            ('Field','field'),
            ('CCD','ccd'),
            ('Followup Status','follow_up_status')
            )
        sort_headers=SortHeaders(request,list_headers,default_order_field=7,default_order_type='desc')
        
        transients=transients.order_by(sort_headers.get_order_by())
        n_trans=transients.count()
        
        rows=[{'html':'','data':['<a href=\'%s\'>%s</a>'%(reverse('trans_info',args=(trans.name,)),trans.name),
                                 trans.type.type,trans.ra,trans.dec,
                                 #trans.first_detected,
                                 #trans.last_detected,
                                 #trans.last_observed,
                                 trans.pointing_det.order_by('dateobs')[0].night,
                                 trans.pointing_det.order_by('-dateobs')[0].night,
                                 trans.pointing_obs.order_by('-dateobs')[0].night,
                                 trans.n_det,
                                 trans.field.id,trans.ccd,
                                 trans.follow_up_status]}
              for trans in transients]
        
        for index,trans in enumerate(transients):
            if 'SN' in trans.type.type or 'Trans' in trans.type.type:
                rows[index]['html']='class=type_sn'
            elif 'AGN' in trans.type.type:
                rows[index]['html']='class=type_agn'
            else:
                rows[index]['html']='class=type_other'
                
        self.headers=list(sort_headers.headers())
        self.rows=rows
        
@login_required
def transients_all(request):
    template='smt_trans_all.html'
    #display all smt transients
    transients=fu.Transient.objects.all().select_related('pointing_det','pointing_obs','field','type')#.annotate(first_detected=Min('pointing_det__dateobs')).annotate(last_detected=Max('pointing_det__dateobs')).annotate(last_observed=Max('pointing_obs__dateobs'))
    n_trans=len(transients)
    typed=[trans for trans in transients if trans.type.type!=u'?']
    n_typed=len(typed)
    
    stats={'n_trans':n_trans,'n_typed':n_typed}
    
    trans_tb=trans_table(request,transients)
    return render_to_response(template,{'stats':stats,'headers':trans_tb.headers,'rows':trans_tb.rows},context_instance=RequestContext(request))

@login_required
def transients_date(request,year,month,day):
    template='smt_trans_date.html'
    #display transients that were detected on selected date
 
    #check if date has observation
    try: 
        obsnight=js.ObsNight.objects.get(ut_date__year=year,ut_date__month=month,ut_date__day=day)
    except ObjectDoesNotExist:
        raise Http404
    

    #get some stats for that obsnight
    pointings=js.SciencePointing.objects.filter(night=obsnight,program='3PS')
    n_pointings=pointings.count()
    n_fields=pointings.values('field').distinct().count()
    
    n_jobs=0
    n_goodjobs=0
    n_refjobs=0
    n_goodrefjobs=0
    for pointing in pointings:
        jobs=js.PipelineJob.objects.filter(jobname__contains='sub',pointing=pointing).select_related('exit_status')
        n_jobs=n_jobs+jobs.count()
        n_goodjobs=n_goodjobs+jobs.filter(exit_status__description='Success').count()
        refjobs=js.PipelineJob.objects.filter(jobname__contains='ref',pointing=pointing).select_related('exit_status')
        n_refjobs=n_refjobs+refjobs.count()
        n_goodrefjobs=n_goodrefjobs+refjobs.filter(exit_status__description='Success').count()
        
    #OK, now all the transients detected on that day
    transients=fu.Transient.objects.filter(pointing_det__night=obsnight).distinct().select_related('pointing_det','pointing_obs','field','type')#.annotate(first_detected=Min('pointing_det__dateobs')).annotate(last_detected=Max('pointing_det__dateobs')).annotate(last_observed=Max('pointing_obs__dateobs'))
    n_trans=len(transients)
    typed=[trans for trans in transients if trans.type.type!=u'?']
    n_typed=len(typed)
    
    stats={'n_pointings':n_pointings,'n_fields':n_fields,
           'n_jobs':n_jobs,'n_goodjobs':n_goodjobs,
           'n_refjobs':n_refjobs,'n_goodrefjobs':n_goodrefjobs,
           'n_trans':n_trans,'n_typed':n_typed}
    
    trans_tb=trans_table(request,transients)
    return render_to_response(template,{'obsnight':obsnight,'stats':stats,'headers':trans_tb.headers,'rows':trans_tb.rows},context_instance=RequestContext(request))

@login_required
def transients_date_latest(request):
    #display transients that were detected on the last Obsnight
    last_night=js.ObsNight.objects.order_by('-ut_date')[0]
    return transients_date(request,last_night.ut_date.year,last_night.ut_date.month,last_night.ut_date.day)
    

@login_required
def transient_info(request,trans_name):
    template='smt_trans_info.html'

    #check if transient name exists
    try: 
        trans=fu.Transient.objects.get(name=trans_name)
    except ObjectDoesNotExist:
        raise Http404
    
    old_type=trans.type.type

    try:
        ptrans=fu.Transient.objects.filter(id__lt=trans.id).reverse()[0]
        prevtrans=ptrans.name
    except: # ObjectDoesNotExist:
        prevtrans=None

    try:
        ntrans=fu.Transient.objects.filter(id__gt=trans.id)[0]
        nexttrans=ntrans.name
    except: # ObjectDoesNotExist:
        nexttrans=None
       
    xsientnames = FilenamesXsient(trans.name,trans.field.id,trans.ccd)

    #get lcdata
    lcfile=xsientnames.abssublcfname
    lctable=pyfits.getdata(lcfile)
    detind=np.where(lctable['FLUX_AP4_ERR'] >=0)[0]
    dettable=lctable[detind]
    limind=np.where(lctable['FLUX_AP4_ERR'] < 0)[0]
    limtable=lctable[limind]
        
    #get all trans types
    trans_types=fu.TransientType.objects.all()
    
    #get history file
    #histfile=xsientnames.abshistfname
    #trans_chist=[]
    #if os.path.exists(histfile):
    #    with open(histfile) as hand:
    #        histlines=hand.readlines()
    #    #parse
    #    for lineind,line in enumerate(histlines):
    #        line=line.strip()
    #        if '#' in line:
    #            userline=histlines[lineind+1]
    #            commentline=histlines[lineind+2]
    #            (ctime,cuser,ctype)=re.match('(.*):\s(.*)\sthinks\s.*\sis\sa\s(.*).\n',userline).groups()
    #            cuser=cuser.replace('user ','')
    #            (ccomm,)=re.match('Comment:\s(.*)\n',commentline).groups()
    #            trans_chist.append({'ctime':ctime.strip(),'cuser':cuser.strip(),'ctype':ctype.strip(),'ccomm':ccomm.strip()})
        

    #user input
    if request.method=="POST":
        #user input type and comment
        user_type=request.POST.get("user_type")
        user_comment=request.POST.get("user_comment")
        user_name=request.user.username
        if user_type in 'other':
            user_type=request.POST.get("user_type_custom")

        #update history file
        #add file lock when writing
        #with open(histfile,'a') as hand:
        #    current_time=datetime.datetime.now().strftime("%d. %B %Y %I:%M %p")
        #    hand.write('#\n')
        #    hand.write('%s: user %s thinks it is a %s.\n'%(current_time,user_name,user_type))
        #    hand.write('Comment: %s \n'%(user_comment))
        #    hand.write('\n')
        
        #update Transient table
        #is this a known type?
        trans_type=fu.TransientType.objects.get_or_create(type=user_type)[0]
        trans.type=trans_type
        trans.savelog(user=request.user,comment=user_comment)
        # FY 2015/08: update the XREF file
        #xreffile=xsientnames.absxrefname
        #hdulist=pyfits.open(xreffile,mode='update')
        #for (key, value) in hdulist[0].header.items():
        #    if value==trans.name:
        #        namematch = re.match("NAME(\d{4})", key)
        #        id=int(namematch.group(1))
        #        hdulist[0].header.update("TYPE{0:04d}".format(id),
        #                                 str(trans.type.type),
        #                                 "Type of source with index {0}".format(id))
        #hdulist.close()
        
        # RS 2014/01/31:  In another awful hack, add to the end of the
        # "followup" file in etc if it's a new Cand or similar.
        # FY 2015/08: use database to track followup
        if user_type in follow_types:
            sf=js.SkymapperField.objects.get(id=trans.field.id)
            if sf.rank==1:
                if sf.label<4:
                    sf.cadence=2.
                    sf.cwidth=sf.cadence*0.1
                    sf.pmult=2.
                    sf.minalt=28.
                    sf.exp='v:120,g:75,r:75,i:60'
                    sf.rank=2
                    sf.save()
                else:
                    sf.cadence=3.
                    sf.cwidth=sf.cadence*0.1
                    sf.pmult=2.
                    sf.minalt=28.
                    sf.exp='v:100,g:100,r:100,i:100'
                    sf.rank=2
                    sf.save()
            # also generate a finding chart for the transient
            # makefinder(trans.name)
        #FY 2015/08: de-followup
        if old_type in follow_types and (not user_type in follow_types):
            sf=js.SkymapperField.objects.get(id=trans.field.id)
            othertypes=list(set(sf.transient_set.all().values_list('type__type',flat=True)))
            goodf=False
            for otype in othertypes:
                if otype in follow_types:
                    goodf=True
            if not goodf and sf.rank>1:
                if sf.label<4:
                    sf.cadence=4.
                    sf.cwidth=1.
                    sf.pmult=1.
                    sf.minalt=40.
                    sf.exp='g:45,r:45,i:30'
                    sf.rank=1
                    sf.save()
                else:
                    sf.cadence=4.
                    sf.cwidth=1.
                    sf.pmult=1.
                    sf.minalt=40.
                    sf.exp='g:100,r:100'
                    sf.rank=1
                    sf.save()

    #changelogs                                                                                                                
    changes=fu.Transient_Change.objects.filter(transient=trans)

    dict_keys=['RUNTAG','FILTER','SM_FIELD','SUBFIELD','SM_CCD','JD_OBS','X_IMAGE','Y_IMAGE','ALPHA_J2000','DELTA_J2000','EXPTIME',
               'FLUX_AP4','FLUX_AP4_ERR','FLAGS','RBSCORE','APCORR04','APCERR04','ZPMAG','ZPMAGERR','MAG_CAL','MAG_CAL_ERR']
    dict_keys.insert(5,'DATE_OBS')
    dict_keys.append('THUMB_NEW')
    dict_keys.append('THUMB_REF')
    dict_keys.append('THUMB_DIFF')
    
    lc_dict = map_lctable_to_dict(dettable,dict_keys,islimit=False,objname=trans.name)
    lim_dict = map_lctable_to_dict(limtable,dict_keys,islimit=True,objname=trans.name)
    lc_dict.update(lim_dict)

    lc_json = simplejson.dumps(lc_dict,skipkeys=True,sort_keys=True)
    return render_to_response(template,{'trans_types':trans_types,'trans':trans,'changes':changes,
                                        #'trans_chist':trans_chist,
                                        'lc_keys':dict_keys,'lcjson':lc_json,
                                        'prevtrans':prevtrans,'nexttrans':nexttrans},context_instance=RequestContext(request))
            
def map_lctable_to_dict(lctable,dict_keys,islimit=False,objname=''):
    #map to a dict object

    #build json to be used for interactive plotting
    lc_json_dict={}
    uniqfilters=list(set(lctable['FILTER']))
    n_uf=len(uniqfilters)
    for uf in uniqfilters:
        if islimit: uf = uf+'_lim'
        lc_json_key=uf
        lc_json_val={'label':uf,'data':[]}
        for dict_key in dict_keys:
              lc_json_val[dict_key.lower()]=[]
        lc_json_dict[lc_json_key]=lc_json_val
    
    hasobsid='OBSID' in lctable.dtype.names
    for lcrow in lctable:
        filter=lcrow['FILTER']
        if islimit: filter=filter+'_lim'
        mjd=float(lcrow['JD_OBS']-2400000.5)
        mag=float(round(lcrow['MAG_CAL']*1e2)/1e2)
        mag_err=float(round(lcrow['MAG_CAL_ERR']*1e2)/1e2)
        if lcrow['MAG_CAL_ERR'] < 0. or lcrow['MAG_CAL_ERR'] > 10. or islimit:
            mag_err=0.
        lc_json_dict[filter]['data'].append([mjd,mag,mag_err])
        if hasobsid:
            subid=FilenamesSub.subid_from_fields(
                runtag=lcrow['RUNTAG'],field=int(lcrow['SM_FIELD']),filter=lcrow['FILTER'],
                obsid=lcrow['OBSID'],obsseq=lcrow['OBSSEQ'], subfield=int(lcrow['SUBFIELD']))
        else:
            subid=FilenamesSub.subid_from_fields(
                runtag=lcrow['RUNTAG'],field=int(lcrow['SM_FIELD']),filter=lcrow['FILTER'],
                obsseq=lcrow['OBSSEQ'], subfield=int(lcrow['SUBFIELD']))
        thumbs=ThumbNames(objname,subid)
        for dict_key in dict_keys:
            if dict_key in 'DATE_OBS':
                #convert jd to datestr
                datetp=ephem.Date(lcrow['JD_OBS']-2415020).tuple()
                dict_val='%04d-%02d-%02dT%02d:%02d:%02d'%tuple(map(round,datetp))
            elif dict_key in 'THUMB_NEW':
                dict_val=thumbs.new
            elif dict_key in 'THUMB_REF':
                dict_val=thumbs.ref
            elif dict_key in 'THUMB_DIFF':
                dict_val=thumbs.diff
                print dict_val
            else:
                dict_val=lcrow[dict_key]
                #make sure values are JSON serializable
                if type(dict_val) is np.float64: 
                    dict_val=float(round(dict_val*1e5)/1e5)
                elif (type(dict_val) is float) or (type(dict_val) is np.float32):
                    dict_val=round(dict_val*1e2)/1e2
                elif type(dict_val) is np.int16:
                    dict_val=int(dict_val)
            lc_json_dict[filter][dict_key.lower()].append(dict_val)
                         
    return lc_json_dict
