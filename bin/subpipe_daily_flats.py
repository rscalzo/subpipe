#!/usr/bin/env python

from Utils.Constants import PipelinePath as CPP
import subprocess
from datetime import datetime,timedelta
import time
from subpipe_runflats import main_runflats
import os

#wait for a bit after sending pause signal
pausetime=600

#pause the main pipeline
while os.path.exists(CPP.pauselockfname):
    print "pause file exists, waiting"
    time.sleep(pausetime)

subprocess.call(['touch',CPP.pauselockfname])
print "Main pipeline paused with lockfile:",CPP.pauselockfname
print "Wait %d seconds."%pausetime
time.sleep(pausetime)

#get yesterday's date
date=(datetime.today()-timedelta(1)).strftime('%Y-%m-%d')

#make flats
print "runflats:",date
main_runflats(date)
#cmd='$SUBPIPEHOME/subpipe_runflats.py %s'%date
#print cmd
#proc=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
#(stdoutstr,stderrstr)=proc.communicate()
#print stdoutstr
#print stderrstr

#for now, just use the latest flat

#in future, might check whether the flat has changed much from last


#unpause main pipeline
subprocess.call(['rm',CPP.pauselockfname])
