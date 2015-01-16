#!/usr/bin/env python

import os
os.environ['DJANGO_SETTINGS_MODULE']='mydjango.SMTP_trainer.settings'

import pyfits
import numpy as np
from mydjango.SMTP_trainer.rb_class.models import Candidate, ImageDate
import glob
import datetime
from STAP.STAP_display import radec2str

subpath = os.environ["SUBSUBPATH"]
newpath = os.environ["SUBNEWPATH"]
refpath = os.environ["SUBREFPATH"]
thumbpath = os.environ["SUBTHUMBS"]

candfiles = glob.glob(subpath+'/*/*/*/*.cands')
print "# of candidate files",len(candfiles)

flags=['xsub','ysub','asub','bsub','esub','thsub','fwhmsub','f4sub','df4sub','f8sub','df8sub','flagsub','starsub',\
       'n2sig3','n3sig3','n2sig5','n3sig5','nmask','Rfwhm','goodcn','subconv',\
       'xref','yref','aref','bref','eref','thref','fwhmref','f4ref','df4ref','flagref','starref','refsrc','nndref',\
       'xnew','ynew','anew','bnew','enew','thnew','fwhmnew','f4new','df4new','flagnew','starnew','newsrc','nndnew',\
       'apsig4','apsig8','normrms','normfwhm','Ranew','Renew','Dthnew','Raref','Reref','Dthref','Rfref','Rfnew']
zp=20.

#cleanup the table
Candidate.objects.all().delete()
ImageDate.objects.all().delete()


lastdate=0
i, imax = 0, len(candfiles)
for candfile in candfiles:
    i += 1
    print "Processing candfile ", candfile, " ({0} of {1})".format(i,imax)
    subname=os.path.basename(candfile).replace('.cands','')
    imname=subname.replace('.fits','')

    dates=map(int,imname.split('_')[-2].split('T')[0].split('-'))
    date=datetime.date(dates[0],dates[1],dates[2])
    if date!=lastdate:
        if ImageDate.objects.filter(date=date).count()==0:
            row=ImageDate.objects.create(date=date)
    lastdate=date

    # RS 2012/05/11:  There are going to be a *lot* of junk detections.
    # In order to make this a less tedious experience for the users, we can
    # instead take a random subset of all Bogus detections.
    keepfrac = 0.1

    try:
        candtb=pyfits.getdata(candfile)
    except:
        print "problem?",candfile
        continue
    for cand in candtb:
        # RS 2012/02/21:  Temp fix for now, gotta get those globs outta here.
        # RS 2012/05/11:  Non-descriptive comment.  I think I enabled this
        # line because there was just *so* much junk and it was senseless to
        # make everyone scan as much junk as I did.  See "keepfrac" above.
        if cand['rbscore'] < 5e-3 and np.random.uniform() > keepfrac:
            continue
        # RS 2012/02/21:  Re-enabled the "_x_y" convention for trainer
        # thumbnail filenames.  It takes too long to calculate (x,y) from
        # (ra,dec) for all the Bogus candidates, when we already even know
        # what (x,y) is for each of them!
        # candname=imname+'_%s'%(radec2str(cand['rasub'],cand['decsub']))
        candname=imname+'_%04.0f_%04.0f'%(cand['xsub'],cand['ysub'])
        dict={'name':candname, 'date':date,
              'thumb_new':imname+'/'+candname+'_new.png',
              'thumb_ref':imname+'/'+candname+'_ref.png',
              'thumb_diff':imname+'/'+candname+'_diff.png'}

        for flag in flags:
            dict[flag]=cand[flag]

        for scoreid in range(10):
            # RS 2012/02/23:  I don't necessarily want new users to see how
            # I've ranked things...
            # exec("dict[\'score%d\'] = %f"%(scoreid+1,cand['rbscore']))
            exec("dict[\'score%d\'] = %f"%(scoreid+1,0))
        exec("dict[\'score%d\'] = %f"%(1,cand['rbscore']))
                    
        row=Candidate.objects.create(**dict)


imagedates=ImageDate.objects.all()
for imagedate in imagedates:
    ncand=Candidate.objects.filter(date=imagedate.date).count()
    imagedate.ncandidate=ncand
    imagedate.save()

    for ind in range(10):
        scoreid="score%d"%(ind+1)
        exec("nscore=Candidate.objects.filter(date=imagedate.date,"+scoreid+"__gt=-1).count()")
        exec("imagedate.n"+scoreid+"=nscore")
        imagedate.save()

    
