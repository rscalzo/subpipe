#! /usr/bin/env python
########################################
# FY - 02/08/2011                      #
# modified to copy data from Miner     #
#                                      #
# FY - 31/05/2011                      #
# this script is to ssh to Markle and  #
# call STAP_find_obsid_files.py        #
# it goes through find_obsid_new.list  #
# and scp each file that is not in the #
# processed file list                  #
########################################


#first write into to_copy_from_markle.list
#remove find_obsid_new.list
#remove from list after copying is successful
  
import os,  re, sys, time, glob
import pyfits
import numpy as np
from STAP_comm import STAP_callexternal

#define home, raw, new and cal path
homepath=os.environ["SUBPIPEHOME"]
rawpath=os.environ["SUBRAWPATH"]
newpath=os.environ["SUBNEWPATH"]
calpath=os.environ["SUBCALPATH"]
scratch=os.environ["SUBSCRATCH"]
scratchpath=scratch+'/raw_to_new/'
flatpath="%s/ff"%calpath

#set filename of new image list
newlist=homepath+"find_obsid_new.list"

#set filename of image queue to be copied from markle
#queue=homepath+"to_copy_from_markle.list"
queue=homepath+"to_copy.list"
#os.system("rm %s"%queue)

#set scp target for STAP_find_obsid_files.py
#target="yuanfang@miner:%s"%(newlist)
target=newlist

#set minimum available disk space, size of five images
minsize=585881280*5

#set scp cmd on markle
#scpmarkle="scp skymap@markle"
scpmarkle="cp"

#define copied list
copied=homepath+"copied.list"
#os.system("rm %s"%copied)

list='/priv/miner3/skymap/fang/subpipe/test_phase_1/SN-Fields.list'
datadir='/priv/miner3/skymap/tisseran/SkyMapperImages/SN-Fields/'
inlist=open(list)
flist=inlist.readlines()
flist=[item.strip() for item in flist if len(item.strip()) >0]
inlist.close()
output=open(newlist,'w')
newtime=time.mktime(time.gmtime())
output.write('%d\n'%newtime)
#for item in flist:
#    [fname,rastr,decstr,filtname,expstr]=re.match('\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S)\s+(\S+)',item).groups()
#    fullpath=datadir+fname
#    output.write(fullpath+'\n')
output.close()

def main():
    #the main loop
    ct=0
    while ct ==0:  
        #now check for the new filelist
        if os.path.exists(newlist):
            temp=open(newlist)
            newtimestamp=temp.readline()
            lines=[line.strip() for line in temp.readlines()]
            temp.close()
            #add files from the queue
            if os.path.exists(queue):
                tempq=open(queue)
                skip=tempq.readline()
                for line in tempq.readlines():
                    lines.append(line.strip())
                tempq.close()
            #update the queue
            tempq=open(queue,'w')
            tempq.write(newtimestamp)
            tempq.writelines(["%s\n"%line for line in lines])
            tempq.close()
            os.remove(newlist)
            #get nocopylist
            try:
                nocopy=open(copied,'r+')
                nocopylist=[line.strip() for line in nocopy.readlines() if
                            re.search('.fits',line)]
            except:
                nocopy=open(copied,'w')
                nocopylist=''
            
            for fullname in lines:
                if re.search("\S+.fits",fullname):
                    file=os.path.basename(fullname)
                    if not file in nocopylist:
                        #check free disk space
                        diskcheck=os.statvfs(rawpath)
                        diskfree=diskcheck.f_frsize*diskcheck.f_bavail
                        if diskfree > minsize:
                            cmd="%s %s %s/"%(scpmarkle,fullname,rawpath)
                            print "cp a new image from Miner:"
                            [status,msg]=STAP_callexternal(cmd,timeout=120,poll=5,combinestderr=True)
                            if status==0:
                                #apply fits64to32_cal
                                [newname,filtname]=fits64to32_cal("%s/%s"%(rawpath,file),scratchpath,flatpath)
                                if len(newname) >0 :
                                    [crval1,crval2]=apply_roughwcs(newname,scratchpath)
                                    if crval1 > -999 and crval2 > -999:
                                        fieldname=crval_to_fieldname(crval1,crval2)
                                        newfullpath='%s/%s/%s'%(newpath,fieldname,filtname)
                                        for ext in range(1,33):
                                            newfullpathext='%s/%02d'%(newfullpath,ext)
                                            if not os.path.exists(newfullpathext):
                                                os.system('mkdir -p %s'%newfullpathext)
                                            cmd='cp %s/%02d/%s_%d.fits %s/'%(scratchpath,ext,newname,ext,newfullpathext)
                                            print cmd
                                            os.system(cmd)
                                            
                                cmd='rm -fr %s/*'%scratchpath
                                print cmd
                                os.system(cmd)

                                #log copied and split files
                                nocopy.write("%s\n"%file)
                                nocopy.flush()
                                #update the queue
                                removefromqueue(queue,fullname)
                            else:
                                print "Failed copying %s from Markle"%fullname
                        else:
                            print "Data disk is nearly full. Will not download more data."
                            
                    else:
                        #update the queue
                        removefromqueue(queue,fullname)
                        print "%s already copied"%file
                else:
                    #update the queue
                    removefromqueue(queue,fullname)
                    print "Do not recognize %s as a fits file"%fullname
            nocopy.close()
        else:
            print "Failed querying for new images on Markle"
            print msg

        ct=ct+1
        cmd='rm -fr %s'%scratchpath
        #print cmd
        os.system(cmd)           
  
def removefromqueue(queue,linetodel):
    tempq=open(queue)
    qlines=tempq.readlines()
    tempq.close()
    tempq=open(queue,'w')
    tempq.writelines(["%s"%line for line in qlines if not line.strip() in linetodel.strip()])
    tempq.close()


def fits64to32_cal(rawname,scratchpath,flatpath):
    #write out to scratch

    newname=os.path.basename(rawname)
    newname=os.path.splitext(newname)[0]
    
    #find out filter name and field center
    temp=pyfits.open(rawname)
    filtname=temp[0].header['FILTNAME']
    if filtname in 'q': filtname='g'
    temp.close()
    #find relavant flat field
    flatbase="%s/%s"%(flatpath,filtname)
    #make sure the flats are there
    flatmiss=0
    for ind in range(1,33):
        if not os.path.exists("%s_%d.fits"%(flatbase,ind)):
            flatmiss=1
                
        if not flatmiss:
            ##call fits64to32 to
            ##remove overscan
            ##apply flat field
            ##align amps
            cmd="/export/miner3/skymap/fang/software/fits64to32_cal %s %s %s %s %s"%(rawname,scratchpath,"-flat",flatbase,"-align")

        else:
            print "Missing/incomplete flat field in %s"%flatbase
            print "Splitting without applying flatfield"
            cmd="/export/miner3/skymap/fang/software/fits64to32_cal %s %s %s"%(rawname,scratchpath,"-align")

        ##Start external call ##
        [status,msg]=STAP_callexternal(cmd,timeout=200,poll=10,combinestderr=True)
        if status !=0 :
            print "Failed to process raw image %s"%rawname
            print msg
            if status ==10:
                print "image might be empty (shutter problem?)"
                cmd='rm -f %s/*/%s_*.fit*'%(scratchpath,newname)
                print cmd
                os.system(cmd)
                return(['',''])
        else:
            return([newname,filtname])
    else:
        print "File %s found in processed directory %s"%(newname,newpath) 
        return(['',''])

def apply_roughwcs(newname,scratchpath):
    
    cmd=homepath+'/WCS/SM-ASTROMETRY.py %s -impath %s'%(newname,scratchpath)
    [status,msg]=STAP_callexternal(cmd,timeout=600,poll=10,combinestderr=True)
    
    if status !=0 :
        print "Failed to apply WCS on image %s"%newname
        print msg
        return([-999,-999])
    else:
        try:
            hdulist=pyfits.open('%s/17/%s_17.fits'%(scratchpath,newname))
            head=hdulist[0].header
            crval1=head['CRVAL1']
            crval2=head['CRVAL2']
            hdulist.close()
            return([crval1,crval2])
        except:
            print "CRVAL not found in the header of image %s"%newname
            return([-999,-999])

def crval_to_fieldname(crval1,crval2):
    crval1h=crval1/15.
    rahour=np.floor(crval1h)
    ramin=np.round((crval1h-rahour)*60.)
    ramin=np.round(ramin/10.)*10
    rastr='%02d%02d'%(rahour,ramin)
    rastr=rastr[0:3]
    if crval2 <0: decsign='-'
    else: decsign='+'
    crval2a=abs(crval2)
    decd=np.round(crval2a)
    decstr='%02d'%decd
    return rastr+decsign+decstr

def field_center_to_name(rastr,decstr):
    raparts=rastr.split(':')
    decparts=decstr.split(':')
    if not decparts[0][0] in ['-','+']:decparts[0]='+'+decparts[0]
    fieldname='SSS'+raparts[0]+raparts[1]+decparts[0]+decparts[1]
    return fieldname

    
if __name__ == "__main__":
        main()
