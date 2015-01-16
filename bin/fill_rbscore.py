import glob,pyfits,os,subprocess

os.environ['DJANGO_SETTINGS_MODULE']='mydjango.SMTP_trainer.settings'

from mydjango.SMTP_trainer.rb_class.models import Candidate
import datetime
import numpy as np
from STAP.STAP_display import radec2str

date=datetime.date(2011,8,25)
db=Candidate.objects.filter(date=date)

subfiles=glob.glob('%s/sub/*/*/*/sub*.fits.cands'%os.environ['SUBPIPEDATA'])

print "# of subfiles:",len(subfiles)
ct=0
ncand=0
for subfile in subfiles:
    #print subfile
    imname=os.path.basename(subfile).replace('.fits.cands','')
    
    hdulist=pyfits.open(subfile,mode='update')
    tb=hdulist[1].data
    ncand=ncand+len(tb)

    for tbrow in tb:
        find=0
        try:
            #candname=imname+'_%04.0f_%04.0f'%(tbrow['xsub'],tbrow['ysub'])
            candname=imname+'_%s'%(radec2str(tbrow['rasub'],tbrow['decsub']))
            cand=db.get(name=candname)
            find=1
        except:
            pass
                    
        if find==1:
            #two scanner at the moment
            score=(cand.score1+cand.score2)/2.
            if score > 0: ct=ct+1
            tbrow['rbscore']=score
        else:
            print "not found:",imname+'_%s'%(radec2str(tbrow['rasub'],tbrow['decsub'])) #imname+'_%04.0f_%04.0f'%(tbrow['xsub'],tbrow['ysub'])

    hdulist.flush()
    hdulist.close()
    
print "# of candidate:",ncand
print '# of score > 0:',ct

