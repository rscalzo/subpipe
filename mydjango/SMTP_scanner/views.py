from django.http import HttpResponse,HttpResponseRedirect, Http404
from django.shortcuts import render_to_response
from django.template import RequestContext
#from django.contrib import auth
from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger

import datetime, os
#from scanner_models.models import Scanner_Date,Scanner_Candidate
from scanner_models.models import Scanner_ObjectType
import pyfits
import numpy as np
import glob

from STAP.STAP_tools.datetime2mjd import mjd2datetime, datetime2mjd
from STAP.STAP_display import radec2str

DAILYDIR=os.environ['CANDDAILY']
CANDPATH=os.environ['CANDPATH']

def requires_login(view):
    def new_view(request, *args, **kwargs):
        if not request.user.is_authenticated():
            return HttpResponseRedirect('/scanner/accounts/login/?next=%s'%request.path)
        return view(request, *args, **kwargs)
    return new_view

@login_required
def SMTP_scanner_home(request):
    template='scanner_home.html'
    #found n_night, n_field and n_candidate from database
    smtp_stats={'n_night':0,'n_field':0,'n_candidate':0}

    dailyfiles=glob.glob('%s/cand_*.fits'%(DAILYDIR))
    smtp_nights=[]
    if len(dailyfiles)>0:
        for dailyfile in dailyfiles:
            datestr=os.path.basename(dailyfile).split('.')[0].split('_')[1]
            smtp_nights.append('-'.join((datestr[:-4],datestr[-4:-2],datestr[-2:])))

    candfiles=glob.glob('%s/cand_*.fits'%(CANDPATH))
    smtp_candidates=[]
    if len(candfiles)>0:
        for candfile in candfiles:
            candname=os.path.basename(candfile).split('.')[0].replace('cand_','')
            smtp_candidates.append(candname)

    if request.method=="POST":
        requestkeys=request.POST.keys()
        for requestkey in requestkeys:
            if requestkey in 'smtp_night':
                night=request.POST.get(requestkey)
                return HttpResponseRedirect('/scanner/%s/%s/%s/'%(tuple(night.split('-'))))
            if requestkey in 'smtp_candidate':
                candname=request.POST.get(requestkey)
                return HttpResponseRedirect('/scanner/cand_%s/'%candname)
    
    return render_to_response(template,{'smtp_stats':smtp_stats,'smtp_nights':smtp_nights,
                                        'smtp_candidates':smtp_candidates},
                              context_instance=RequestContext(request))

@login_required
def SMTP_scanner_date(request,year,month,day):
    template='scanner_date.html'
    valid=1
    try:
        intyear=int(year)
        intmonth=int(month)
        intday=int(day)
        if intyear < 2011 or intyear > 2015:valid=0
        if intmonth <1 or intmonth > 12:valid=0
        if intday <1 or intday >31: valid=0
    except:
        valid=0
        
    if valid==0:
        raise Http404


    #access database for some stats
    n_exp=0
    n_sub=0
    #==============================
    #date=datetime.date(intyear,intmonth,intday)
    #dateobject=Scanner_Date.objects.filter(date=date)
    #==============================

    datestr='%04d%02d%02d'%(intyear,intmonth,intday)
    datedisplay='%04d-%02d-%02d'%(intyear,intmonth,intday)

    #get daily candidate table
    dailyfile='%s/cand_%s.fits'%(DAILYDIR,datestr)
    if not os.path.exists(dailyfile):
        raise Http404        
    dailycand=pyfits.getdata(dailyfile)

    #make daily stats
    ncand=len(dailycand)
    nscan=0
    for candind in range(ncand):
        if dailycand[candind]['type_best_user'].strip() not in '': nscan=nscan+1
            
    dailystats={'date_display':datedisplay,'n_exposure':n_exp,
                'n_subtraction':n_sub,'n_candidate':ncand,
                'n_scanned':nscan}
    
    #sort
    list_headers = (
        ('rbscore','rbscore'),
        ('name','name'),
        ('ra','ra'),
        ('dec','dec'),
        ('filter','filter'),
        ('mag','mag'),
        ('mag_err','mag_err'),
        ('fwhm','fwhm'),
        ('field','field'),
        ('ccd','ccd'),
        ('type_auto','type_auto'),
        ('type_best','type_best'),
        ('type_best_user','type_best_user'),
        )
    sort_headers=SortHeaders(request,list_headers,default_order_field=0,default_order_type='desc')
    sortbystr=sort_headers.get_order_by()
    if sortbystr[0]=='-':
        sortedindex=dailycand[sortbystr[1:]].argsort()[::-1]
    else:
        sortedindex=dailycand[sortbystr].argsort()

    #render
    return render_to_response(template,{'dailycand':dailycand[sortedindex],'dailystats':dailystats,
                                        'headers':list(sort_headers.headers())},
                              context_instance=RequestContext(request))

@login_required
def SMTP_scanner_candidate(request,candname):
    template='scanner_candidate.html'

    #find candfile 
    candfile='%s/cand_%s.fits'%(CANDPATH,candname)
    if not os.path.exists(candfile):
        raise Http404

    #open variable length table, can't get data at once
    hdulist=pyfits.open(candfile)
    candtb=hdulist[1].data

    mjdset=list(set(np.floor(candtb[0]['mjd'])))
    mjdset.sort()
    mjdset.reverse()
    nights=[]
    for mjd in mjdset:
        nights.append(mjd2datetime(mjd))

    type_best=''.join(list(candtb[0]['type_best']))
    if type_best in '': type_best=''.join(list(candtb[0]['type_auto']))
    else: type_best=type_best.split('|')[-2]
    candstats={'name':candtb[0]['name'],'n_night':candtb[0]['n_night'],
               'ra_best':candtb[0]['ra_best'],'dec_best':candtb[0]['dec_best'],
               'lc':'lc_figs/cand_%s_lc.png'%candname,
               'type_best':type_best,
               'nights':nights}

    pagedate=request.GET.get('page')
    try:
        pagedate=pagedate.replace('/','')
        pageparts=map(int,pagedate.split('-'))
        pagemjd=datetime2mjd(datetime.datetime(pageparts[0],pageparts[1],pageparts[2],0,0,0))
    except:
        pagemjd=mjdset[0]

    candind=np.where(np.floor(candtb[0]['mjd'])==pagemjd)[0]
    if len(candind) ==0:
        pagemjd=mjdset[0]
        candind=np.where(np.floor(candtb[0]['mjd'])==pagemjd)[0]

    candstats['shownight']=mjd2datetime(pagemjd)

    canddata=[]
    allsubname=''.join(list(candtb[0]['subname']))
    for candindi in candind:
        radecstr=radec2str(candtb[0]['ra'][candindi],candtb[0]['dec'][candindi])
        subname='sub'+allsubname.split('sub')[candindi+1]
        thumbname='thumbs/%s/%s_%s'%(subname,subname,radecstr)
        canddata.append({'mjd':candtb[0]['mjd'][candindi],'filter':candtb[0]['filter'][candindi],
                         'mag':candtb[0]['mag'][candindi],'mag_err':candtb[0]['mag_err'][candindi],'mlim':candtb[0]['mlim'][candindi],
                         'fwhm':candtb[0]['fwhm'][candindi],
                         'thumb_new':'%s_new.png'%thumbname,
                         'thumb_ref':'%s_ref.png'%thumbname,
                         'thumb_diff':'%s_diff.png'%thumbname})
    
    hdulist.close()

    if request.method=="POST":
        type=str(request.POST.get("usertype"))+'|'
        comment=str(request.POST.get("usercomment"))+'|'
        #update candfile
        hdulist=pyfits.open(candfile,mode='update')
        hdulist[1].data[0]['type_best']=''.join(hdulist[1].data[0]['type_best'])+type
        hdulist[1].data[0]['type_best_user']=''.join(hdulist[1].data[0]['type_best_user'])+str(request.user.username)+'|'
        hdulist[1].data[0]['type_best_user_comment']=''.join(hdulist[1].data[0]['type_best_user_comment'])+comment
        hdulist[1].update()
        #refresh type_best
        type_best=''.join(list(hdulist[1].data[0]['type_best']))
        if type_best in '': type_best=''.join(list(hdulist[1].data[0]['type_auto']))
        else: type_best=type_best.split('|')[-2]
        candstats['type_best']=type_best
        #update dailyfile?
        hdulist.flush()
        hdulist.close()
        submitstr=str(request.POST.get("submit"))
        if 'go back' in submitstr:
            date=submitstr.split()[-1].split('-')
            if len(date)==3:
                return HttpResponseRedirect('/scanner/%s/%s/%s'%(date[0],date[1],date[2]))
    
    objtypes=Scanner_ObjectType.objects.all()
    if len(objtypes)==0:
        objtypes=[{'type':'SN'},{'type':'AGN'},{'type':'VarStar'},{'type':'Asteroid'},{'type':'Junk'}]
        
    return render_to_response(template,{'candstats':candstats,'canddata':canddata,'objtypes':objtypes},
                              context_instance=RequestContext(request))

def SMTP_scanner_register(request):
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
                    message = "Congratulations! You have been registered. You can <a href=\"/scanner/\"> start scanning the candidates </a> now."
                else:
                    message = "urh... you account might be disabled."
            else:
                message = "There was an error automatically logging you in. Try <a href=\"/scanner/login/\"> logging in </a> manually."
                
            return render_to_response('scanner_register_success.html',{
                'username':username,
                'message':message,},context_instance=RequestContext(request))
    else:        
        form=UserCreationForm()            
        
    template='scanner_register.html'
    return render_to_response(template,{'form':form},context_instance=RequestContext(request))
    

ORDER_VAR = 'o'
ORDER_TYPE_VAR = 'ot'
class SortHeaders:
    """
    Handles generation of an argument for the Django ORM's
    ``order_by`` method and generation of table headers which reflect
    the currently selected sort, based on defined table headers with
    matching sort criteria.

    Based in part on the Django Admin application's ``ChangeList``
    functionality.
    """
    def __init__(self, request, headers, default_order_field=None,
            default_order_type='asc', additional_params=None):
        """
        request
            The request currently being processed - the current sort
            order field and type are determined based on GET
            parameters.

        headers
            A list of two-tuples of header text and matching ordering
            criteria for use with the Django ORM's ``order_by``
            method. A criterion of ``None`` indicates that a header
            is not sortable.

        default_order_field
            The index of the header definition to be used for default
            ordering and when an invalid or non-sortable header is
            specified in GET parameters. If not specified, the index
            of the first sortable header will be used.

        default_order_type
            The default type of ordering used - must be one of
            ``'asc`` or ``'desc'``.

        additional_params:
            Query parameters which should always appear in sort links,
            specified as a dictionary mapping parameter names to
            values. For example, this might contain the current page
            number if you're sorting a paginated list of items.
        """
        if default_order_field is None:
            for i, (header, query_lookup) in enumerate(headers):
                if query_lookup is not None:
                    default_order_field = i
                    break
        if default_order_field is None:
            raise AttributeError('No default_order_field was specified and none of the header definitions given were sortable.')
        if default_order_type not in ('asc', 'desc'):
            raise AttributeError('If given, default_order_type must be one of \'asc\' or \'desc\'.')
        if additional_params is None: additional_params = {}

        self.header_defs = headers
        self.additional_params = additional_params
        self.order_field, self.order_type = default_order_field, default_order_type

        # Determine order field and order type for the current request
        params = dict(request.GET.items())
        if ORDER_VAR in params:
            try:
                new_order_field = int(params[ORDER_VAR])
                if headers[new_order_field][1] is not None:
                    self.order_field = new_order_field
            except (IndexError, ValueError):
                pass # Use the default
        if ORDER_TYPE_VAR in params and params[ORDER_TYPE_VAR] in ('asc', 'desc'):
            self.order_type = params[ORDER_TYPE_VAR]

    def headers(self):
        """
        Generates dicts containing header and sort link details for
        all defined headers.
        """
        for i, (header, order_criterion) in enumerate(self.header_defs):
            th_classes = []
            new_order_type = 'asc'
            if i == self.order_field:
                th_classes.append('sorted %sending' % self.order_type)
                new_order_type = {'asc': 'desc', 'desc': 'asc'}[self.order_type]
            yield {
                'text': header,
                'sortable': order_criterion is not None,
                'url': self.get_query_string({ORDER_VAR: i, ORDER_TYPE_VAR: new_order_type}),
                'class_attr': (th_classes and ' class="%s"' % ' '.join(th_classes) or ''),
            }

    def get_query_string(self, params):
        """
        Creates a query string from the given dictionary of
        parameters, including any additonal parameters which should
        always be present.
        """
        params.update(self.additional_params)
        return '?%s' % '&amp;'.join(['%s=%s' % (param, value) \
                                     for param, value in params.items()])

    def get_order_by(self):
        """
        Creates an ordering criterion based on the current order
        field and order type, for use with the Django ORM's
        ``order_by`` method.
        """
        return '%s%s' % (
            self.order_type == 'desc' and '-' or '',
            self.header_defs[self.order_field][1],
        )
    
