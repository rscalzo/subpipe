#!/usr/bin/env python

from Utils.Constants import PipelinePath as CPP
import subprocess
from datetime import datetime,timedelta
import time

#wait for a bit after sending pause signal
pausetime=600

#pause the main pipeline
subprocess.call(['touch',CPP.pauselockfname])
print "Main pipeline paused with lockfile:",CPP.pauselockfname
print "Wait %d seconds."%pausetime
time.sleep(pausetime)

#get yesterday's date
date=(datetime.today()-timedelta(1)).strftime('%Y-%m-%d')

#make flats
cmd='$SUBPIPEHOME/subpipe_runflats.py %s'%date
print cmd
proc=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
stdoutstr=proc.communicate()

#for now, just use the latest flat

#in future, might check whether the flat has changed much from last


#unpause main pipeline
subprocess.call(['rm',CPP.pauselockfname])
