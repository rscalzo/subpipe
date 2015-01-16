# Create your views here.

from django.shortcuts import render_to_response
from django.template import RequestContext

import re
from Utils.Constants import PipelinePath as CPP, FilenamesSub, FilenamesRef
import mydjango.jobstats.models as js
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
