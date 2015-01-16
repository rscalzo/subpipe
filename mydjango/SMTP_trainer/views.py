from django.http import HttpResponse,HttpResponseRedirect, Http404
from django.shortcuts import render_to_response
from django.template import RequestContext
from rb_class.models import Users, Candidate, ImageDate
from django.contrib import auth
from django.contrib.auth.decorators import login_required
from django import forms
from django.contrib.auth.forms import UserCreationForm
import datetime
from django.core.paginator import Paginator, EmptyPage, PageNotAnInteger

def requires_login(view):
    def new_view(request, *args, **kwargs):
        if not request.user.is_authenticated():
            return HttpResponseRedirect('/accounts/login/?next=%s'%request.path)
        return view(request, *args, **kwargs)
    return new_view

@login_required
def SMTP_rb_trainer(request):
    template='trainer_rb.html'
    try:
        currentuser=Users.objects.get(name=request.user.username)
    except:
        #user is not in the database and thus can not submit data
        return HttpResponseRedirect('/accounts/login/?next=%s'%request.path)

    imagedates=ImageDate.objects.all()
    for imagedate in imagedates:
        userid=str(currentuser.userid)
        exec("nscore=Candidate.objects.filter(date=imagedate.date,score"+userid+"__gt=-1).count()")
        exec("imagedate.nscore"+userid+"=nscore")
        imagedate.save()

    return render_to_response(template,{'trainer_dates':imagedates,'userid':userid},
                              context_instance=RequestContext(request))

@login_required
def SMTP_rb_trainer_date(request,year,month,day):
    template='trainer_rb_date.html'
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

    try:
        currentuser=Users.objects.get(name=request.user.username)
    except:
        #user is not in the database and thus can not submit data
        return HttpResponseRedirect('/accounts/login/?next=%s'%request.path)

    date=datetime.date(intyear,intmonth,intday)
    candidates=Candidate.objects.filter(date=date)
        
    if candidates.count() ==0:
        raise Http404

    paginator = Paginator(candidates, 50)

    page = request.GET.get('page')
    try:
        candspage = paginator.page(page)
    except PageNotAnInteger:
        # If page is not an integer, deliver first page.
        candspage = paginator.page(1)
    except EmptyPage:
        # If page is out of range (e.g. 9999), deliver last page of results.
        candspage = paginator.page(paginator.num_pages)
        
    userid=str(currentuser.userid)
    if request.method=="POST":
        scorekey=request.POST.keys()
        scorekey=[item for item in scorekey if "score" in item][0]
        scoreid=scorekey.split('_')[-1]
        candname=scorekey.replace('_'+scoreid,'')
        score=int(request.POST.get(scorekey,"-1"))
        for candidate in candidates: 
#        for candidate in candspage:
            if candidate.name in candname: 
                exec("candidate."+scoreid+" = %d "%score)
                candidate.save()

    return render_to_response(template,{'trainer_candidates':candspage,'userid':userid,'date':date},
                              context_instance=RequestContext(request))

def SMTP_trainer_register(request):
    if request.method=="POST":
        #add user account
        form = UserCreationForm(request.POST)
        if form.is_valid():
            new_user=form.save()
            if not Users.objects.filter(name=new_user.username):
                if Users.objects.count()>0:
                    lastuser=Users.objects.order_by('-userid')[0]
                    lastid=lastuser.userid
                else: lastid=0
                if lastid <10:
                    row=Users.objects.create(name=new_user.username,userid=lastid+1)
                    return HttpResponseRedirect('/trainer/')
                else:
                    print "enough classifiers"
                    return HttpResponseRedirect('/accounts/login/?next=/trainer/')
            else:
                return HttpResponseRedirect('/trainer/')
        else:
            print "form invalid"
    else:        
        form=UserCreationForm()            
    
    template='trainer_register.html'
    return render_to_response(template,{'form':form},context_instance=RequestContext(request))
    
