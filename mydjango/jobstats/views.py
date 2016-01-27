# Create your views here.

from django.shortcuts import render_to_response
from django.template import RequestContext
from django.contrib.auth.decorators import login_required

import re
from Utils.Constants import PipelinePath as CPP, FilenamesSub, FilenamesRef
import mydjango.jobstats.models as js
import mydjango.followup.models as fu
from mydjango.jobstats.SortHeaders import SortHeaders
from django.core.urlresolvers import reverse


def index(request):
    return render_to_response("base.html", {})

def listjobsbyrun(request, runtag):
    """Renders a list of jobs by run.
    
    RS:  At the moment I've designed things so that you must specify a run,
    because that's likely to be a common use case and because the query is
    run using the reverse relation for PipelineRun into PipelineJob.
    I can always lift this restriction if we think it'd be useful."""

    # Put the column headers into a dict form the template can understand.
    colheads = [ ( "Job Name", "jobname" ),
                 ( "Start Time", "start_time" ),
                 ( "Wall Time", "walltime" ),
                 ( "Exit Stage", "exit_status__stage_fail" ),
                 ( "Status", "exit_status__description" ), ]

    # Sort out any other filters we might have specified.
    allowed_filters = [attr for (name, attr) in colheads]
    allowed_filters += ['walltime__gt', 'walltime__lt']
    jobfilters = dict([(key, val) for (key, val) in request.GET.items()
                       if key in allowed_filters])

    # Get all the jobs associated with a particular run, and use the
    # "order_by" field to order the QuerySet accordingly.
    sort_headers = SortHeaders(request, colheads,
                               default_order_field=0,
                               default_order_type='asc',
                               additional_params=jobfilters)
    run = js.PipelineRun.objects.get(runtag=runtag)
    job_list = run.pipelinejob_set.select_related('exit_status')
    if len(jobfilters) < 1:  job_list = job_list.all()
    else:  job_list = job_list.filter(**jobfilters)
    job_list = job_list.order_by(sort_headers.get_order_by())

    # Put the data themselves into a dict form the template can understand.
    rows = [ ({ 'html': "class=nowrap", 'data': job.jobname, },
              { 'html': "class=nowrap", 'data': job.start_time, },
              { 'html': "class=nowrap", 'data': "{0:.1f} sec".format(job.walltime), },
              { 'html': "class=nowrap", 'data': job.exit_status.stage_fail, },
              { 'data': job.exit_status.description }) for job in job_list ]

    # Add some useful links to the data in the columns.
    for r in rows:
        # The user can click on a job name to view the log file.
        r[0]['url'] = reverse('joblogs',args=(r[0]['data'],))
        # The user can filter on failure conditions just by clicking on them.
        for i in (3, 4):
            tmpdict = dict(request.GET.items())
            tmpdict.update({ allowed_filters[i]: r[i]['data'] })
            r[i]['url'] = '?{0}'.format('&amp;'.join(
                ['{0}={1}'.format(key, val) for key, val in tmpdict.items()]))

    # Render the page
    return render_to_response("listjobsby_run.html",
            { "colheads": colheads, "rows": rows, "runtag": runtag,
              "headers": list(sort_headers.headers()),
              "jobfilters" : jobfilters },
              context_instance=RequestContext(request))

def displaylog(request, jobname):
    """Displays the log file for the given subtraction."""

    # Figure out where the log file is on disk.  Then just read the whole
    # darn file in and display it.
    if 'ref' in jobname:
        sub=FilenamesRef('',subid=jobname)
    else:
        sub = FilenamesSub(subid=jobname)
    logcontent = [ ]
    with open(sub.absolute_dir + '/' + sub.log_file_name) as logfile:
        for line in logfile:  logcontent.append(line.strip())
    return render_to_response("displaylog.html",
            { "jobname": jobname, "logcontent": logcontent },context_instance=RequestContext(request))


@login_required
def smt_field(request,fid):
    pass
    template='smt_field.html'
    #find field by id
    #update schedule requirement
    #exp,cadence,cwidth,minalt,pmult
    try:
        sf=js.SkymapperField.objects.get(id=fid)
    except ObjectDoesNotExist:
        raise Http404
  
    return render_to_response(template,{'field':sf},context_instance=RequestContext(request))

def get_ra(request_ra):
    ra=None
    try:
        trans_ra=float(request_ra)
    except:
        try:
            if ':' in request_ra:
                parts=map(float,request_ra.split(':'))
            else:
                parts=map(float,request_ra.split())
            trans_ra=(parts[0]+parts[1]/60.+parts[2]/3600.)*15.
        except:
            pass
    return trans_ra

def get_dec(request_dec):
     try:
         trans_dec=float(request_dec)
     except:
         sign=1.
         if '-' in request_dec[0]:
             sign=-1.
             request_dec=request_dec[1:]
         try:
             if ':' in request_dec:
                 parts=map(float,request_dec.split(':'))
             else:
                 parts=map(float,request_dec.split())
             trans_dec=(parts[0]+parts[1]/60.+parts[2]/3600.)*sign
         except:
             pass 
     return trans_dec


CCDW=2048.*0.5/3600.
CCDH=4096.*0.5/3600.
GAPX=101.*0.5/3600.
GAPY0=64.*0.5/3600.
GAPY1=383.*0.5/3600.

def xyoffset(ccd_offset):
    if ccd_offset[0]<0:
        xoffset=(ccd_offset[0]+0.5)*(CCDW+GAPX)
    else:
        xoffset=(ccd_offset[0]-0.5)*(CCDW+GAPX)
    if ccd_offset[1]==-2:
        yoffset=(ccd_offset[1]+0.5)*CCDH-GAPY0/2.-GAPY1
    elif ccd_offset[1]==-1:
        yoffset=(ccd_offset[1]+0.5)*CCDH-GAPY0/2.
    elif ccd_offset[1]==1:
        yoffset=(ccd_offset[1]-0.5)*CCDH+GAPY0/2.
    else:
        yoffset=(ccd_offset[1]-0.5)*CCDH+GAPY0/2.+GAPY1
       
    return xoffset,yoffset

CCD_XYOFF=[(-4,-1),(-4,-2),(+1,-1),(+1,-2),(-3,-1),(-3,-2),(+2,-1),(+2,-2),
           (-2,-1),(-2,-2),(+3,-1),(+3,-2),(-1,-1),(-1,-2),(+4,-1),(+4,-2),
           (+1,+1),(+1,+2),(-4,+1),(-4,+2),(+2,+1),(+2,+2),(-3,+1),(-3,+2),
           (+3,+1),(+3,+2),(-2,+1),(-2,+2),(+4,+1),(+4,+2),(-1,+1),(-1,+2)]
CCD_OFFSETS=[xyoffset(offset) for offset in CCD_XYOFF]

@login_required
def smt_field_lookup(request):
    import ephem
    import numpy as np

    template='smt_field_lookup.html'
    #ra/dec search

    #default
    trans_ra=None
    trans_dec=None
    trans_rad=1./3600.
    
    f_ra=None
    f_dec=None
    f_rad=1.25

    translist=[]
    fields=[]

    if request.method=="POST":
        requestkeys=request.POST.keys()
        if "trans_ra" in requestkeys:
            request_ra=request.POST.get("trans_ra")
            trans_ra=get_ra(request_ra)
        if "trans_dec" in requestkeys: 
            request_dec=request.POST.get("trans_dec")
            trans_dec=get_dec(request_dec)
        if "trans_rad" in requestkeys: trans_rad=float(request.POST.get("trans_rad"))
        if "f_ra" in requestkeys: f_ra=get_ra(request.POST.get("f_ra"))
        if "f_dec" in requestkeys: f_dec=get_dec(request.POST.get("f_dec"))
        if "f_rad" in requestkeys: f_rad=float(request.POST.get("f_rad"))
        
    if any([trans_ra is None, trans_dec is None, trans_rad is None]):
        translist=[]
    else:
        SMtrans=list(fu.Transient.objects.all())
        translist=[]
        ephra=ephem.degrees(str(trans_ra))
        ephdec=ephem.degrees(str(trans_dec))
        for f in SMtrans:
            arcdist=ephem.separation((ephra,ephdec),(ephem.degrees(str(f.ra)), ephem.degrees(str(f.dec)))) * 180/ephem.pi
            if arcdist<=trans_rad:
                f.dist=arcdist
                translist.append(f)
    
    if any([f_ra is None, f_dec is None, f_rad is None]):
        fields=[]
    else:
        SMfield_list=list(js.SkymapperField.objects.all())
        fields=[]
        ephra=ephem.degrees(str(f_ra))
        ephdec=ephem.degrees(str(f_dec))
        for f in SMfield_list:
            arcdist=ephem.separation((ephra,ephdec),(ephem.degrees(str(f.ra)), ephem.degrees(str(f.dec)))) * 180/ephem.pi
            if arcdist<=f_rad:
                f.dist=arcdist
                #look up which ccd it's on
                xoff=(f_ra-f.ra)*np.cos(f_dec*ephem.pi/180.)
                yoff=(f_dec-f.dec)
                ccd_ind=[ccdid+1 for ccdid, xyoff in enumerate(CCD_OFFSETS)
                         if abs(xoff-xyoff[0])< CCDW/2. and abs(yoff-xyoff[1])<CCDH/2. ]
                if len(ccd_ind)>0:
                    f.ccd=ccd_ind[0]
                else:
                    f.ccd=0
                fields.append(f)
                
    return render_to_response(template,{'trans_ra':trans_ra,'trans_dec':trans_dec,'trans_rad':trans_rad,
                                        'f_ra':f_ra,'f_dec':f_dec,'f_rad':f_rad,
                                        'translist':translist,'fieldlist':fields},context_instance=RequestContext(request))

