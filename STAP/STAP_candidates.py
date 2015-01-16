#!/usr/bin/env python

import os

#from mydjango.SMTP_scanner.scanner_models.models import Scanner_Date,Scanner_Candidate_Good,Scanner_ObjectType, Scanner_Candidate_Lookup
import datetime
import pyfits
import numpy as np
from STAP.STAP_tools.close_match_radec import close_match_radec


def STAP_banner(step,screenwidth=80):
    # Creates nicely formatted banner to place between steps.
    myname = "Executing step:  "+step
    mypad1 = " "*int((screenwidth-4-len(myname))/2)
    mypad2 = " "*(screenwidth-4-len(myname+mypad1))
    banner = ("+"+("-"*(screenwidth-4))+"+\n"
           +  "|" + mypad1 + myname + mypad2 + "|\n"
           +  "+"+("-"*(screenwidth-4))+"+\n")
    return banner


import matplotlib
#matplotlib.use('Agg')
from matplotlib import pyplot

from STAP.STAP_tools.datetime2mjd import mjd2datetime

RB_CUTOFF=40.
ZP={'v':22.0,'g':24.0,'r':23.0,'i':23.0}
PIXSCALE=0.5/3600.

def update_all_cand_tables(candfile,log=None):
    #update cand_nightid.fits
    #update cand_name.fits
    #update cand lookup

    banner=STAP_banner('update_cand_tables')
    print_msg(banner,log=log)
            
    dailydir=os.environ['CANDDAILY']
    lookupfile=os.environ['CANDLOOKUP']
    candpath=os.environ['CANDPATH']
    lcpath=os.environ['LCFIGPATH']

    if not os.path.exists(candfile):
        errmsg='ERROR: %s does not exist.'%candfile
        print_msg(errmsg,log=log)
        return -1
        
    hdulist=pyfits.open(candfile)
    candhead=hdulist[0].header
    canddata=hdulist[1].data
    hdulist.close()
    if len(canddata) ==0:
        errmsg='WARNING: no candidate found in %s.'%candfile
        print_msg(errmsg,log=log)
        return -1
            
    #real-bogus cut
    real=np.where(canddata['rbscore']>=RB_CUTOFF)[0]
    nnewcand=len(real)
    if nnewcand ==0:
        errmsg='NOTE: no real candidate (passed rbscore test) found in %s.'%candfile
        print_msg(errmsg,log=log)
        return -1
       
    #=======================
    #start updating.....
    canddata=canddata[real]
    newra=canddata['rasub']
    newdec=canddata['decsub']
    
    #get subname
    subname=os.path.basename(candfile).split('.')[0]
    #get nightid
    nightid=nightid_from_dateobs(candhead['DATE-OBS'])

    #update cand_night.fits
    dailyfile='%s/cand_%s.fits'%(dailydir,nightid)
    if not os.path.exists(dailyfile):
        #create
        print_msg('Creating %s'%dailyfile,log=log)
        newtbhdu=init_daily_table(nnewcand)
        newtbhdu,newindex=fill_daily_table(newtbhdu,canddata,candhead,nnewcand,first=0,subname=subname)
        newtbhdu.writeto(dailyfile,clobber=True)
    else:
        #update
        hdulist=pyfits.open(dailyfile,mode='update')
        m2,m1,miss=close_match_radec(newra,newdec,hdulist[1].data['ra'],hdulist[1].data['dec'],PIXSCALE*1.5,allow=1,silent=1)
        nrow1 = hdulist[1].data.shape[0]
        nmiss=len(miss)
        nmatch=len(m1)
        newtbhdu=pyfits.new_table(hdulist[1].columns,nrows=nrow1+nmiss)
        #update tb[m2]
        newindex=[]
        if nmatch >0:
            print_msg('Updating %s'%dailyfile,log=log)
            newtbhdu,newindex1=update_daily_table(newtbhdu,canddata[m2],candhead,nmatch,tbindex=m1,subname=subname)
            newindex=np.array(np.append(newindex,newindex1),dtype='int')
        #add canddata[miss]1
        if nmiss >0:
            print_msg('Filling %s'%dailyfile,log=log)
            newtbhdu,newindex2=fill_daily_table(newtbhdu,canddata[miss],candhead,nmiss,first=nrow1,subname=subname)
            newindex=np.array(np.append(newindex,newindex2),dtype='int')
        hdulist[1]=newtbhdu
        hdulist.flush()
        hdulist.close()

          
    #=======================
    #match lookup table
    #only new entries in the daily table
    dailyhdulist=pyfits.open(dailyfile,mode='update')
    dailytb=dailyhdulist[1].data
    if not os.path.exists(lookupfile):
        #start a lookupfile
        init_lookup(lookupfile)
        m1=[]
        miss=np.arange(len(newindex))
        lookup=pyfits.open(lookupfile,mode='update')
        lenlookup=0
    else:
        lookup=pyfits.open(lookupfile,mode='update')
        if lookup[1].data is None:
            #empty lookup, everything is new
            m1=[]
            miss=np.arange(len(newindex))
            lenlookup=0
        else:
            lenlookup=len(lookup[1].data)
            #match
            m2,m1,miss=close_match_radec(dailytb['ra'][newindex],dailytb['dec'][newindex],lookup[1].data['ra'],lookup[1].data['dec'],PIXSCALE*1.5,allow=1,silent=1)
            
    nmatch=len(m1)
    nmiss=len(miss)        
    if nmatch >0:
        for m2ind,mindex in enumerate(m1):
            dailyindex=[newindex[m2[m2ind]]]
            #update cand_name for matched ones
            candnamefile='%s/cand_%s.fits'%(candpath,lookup[1].data[mindex]['name'])
            print_msg('Updating %s'%candnamefile,log=log)
            candlcfile='%s/cand_%s_lc.png'%(lcpath,lookup[1].data[mindex]['name'])
            hdulist=pyfits.open(candnamefile,'update')
            candtbhdu=hdulist[1]
            candtbhdu=update_cand_name(candtbhdu,dailytb[dailyindex])
            plot_cand_lc(candtbhdu.data[0]['mjd'],candtbhdu.data[0]['mag'],candtbhdu.data[0]['mag_err'],candtbhdu.data[0]['mlim'],candtbhdu.data[0]['filter'],candlcfile)
            hdulist.flush()
            hdulist.close()
            #update cand_nightid.fits with name, type_auto, type_best, type_user, ra/dec_best
            print_msg('Updating %s'%dailyfile,log=log)
            dailytb[dailyindex]['name']=candtbhdu.data[0]['name']
            dailytb[dailyindex]['ra_best']=candtbhdu.data[0]['ra_best']
            dailytb[dailyindex]['dec_best']=candtbhdu.data[0]['dec_best']
            dailytb[dailyindex]['type_auto']=get_last_part(candtbhdu.data[0]['type_auto'])
            dailytb[dailyindex]['type_best']=get_last_part(candtbhdu.data[0]['type_best'])
            dailytb[dailyindex]['type_best_user']=get_last_part(candtbhdu.data[0]['type_best_user'])
            #update lookup tb with ra/dec_best
            #no extension of lookup table
            lookup[1].data[mindex]['ra']=candtbhdu.data[0]['ra_best']
            lookup[1].data[mindex]['dec']=candtbhdu.data[0]['dec_best']
    if nmiss >0:
        newlookuphdu=pyfits.new_table(lookup[1].columns,nrows=lenlookup+nmiss)
        for lind,mindex in enumerate(miss):
            dailyindex=[newindex[mindex]]
            #create cand_name for unmatched ones
            candtbhdu=init_cand_name(dailytb[dailyindex])
            candname=candtbhdu.data[0]['name']
            candnamefile='%s/cand_%s.fits'%(candpath,candname)
            print_msg('Creating %s'%candnamefile,log=log)
            candlcfile='%s/cand_%s_lc.png'%(lcpath,candname)
            plot_cand_lc(candtbhdu.data[0]['mjd'],candtbhdu.data[0]['mag'],candtbhdu.data[0]['mag_err'],candtbhdu.data[0]['mlim'],candtbhdu.data[0]['filter'],candlcfile)
            candtbhdu.writeto(candnamefile,clobber=True)            
            #update cand_nightid.fits
            print_msg('Updating %s'%dailyfile,log=log)
            dailytb[dailyindex]['name']=candtbhdu.data[0]['name']
            dailytb[dailyindex]['ra_best']=candtbhdu.data[0]['ra_best']
            dailytb[dailyindex]['dec_best']=candtbhdu.data[0]['dec_best']
            #update lookup tb with name, ra/dec_best
            #extend lookup table
            newlookuphdu.data[lenlookup+lind]['name']=candtbhdu.data[0]['name']
            newlookuphdu.data[lenlookup+lind]['ra']=candtbhdu.data[0]['ra_best']
            newlookuphdu.data[lenlookup+lind]['dec']=candtbhdu.data[0]['dec_best']
        lookup[1]=newlookuphdu
    lookup.flush()
    lookup.close()
    dailyhdulist.flush()
    dailyhdulist.close()
    
    print_msg('Done',log=log)
    

def print_msg(msg,log=None):
    if log is None:
        print msg
    else:
        mystdout=open(log,"a")
        mystdout.write('%s\n'%msg)
        mystdout.close()
    

def nightid_from_dateobs(dateobs):
    return ''.join(dateobs.split('T')[0].split('-'))

def init_daily_table(nrows):
    collist=[{"N":"name","F":"20A"},
             {"N":"ra_best","F":"D"},
             {"N":"dec_best","F":"D"},
             {"N":"mjd","F":"D"},
             {"N":"filter","F":"1A"},
             {"N":"mag","F":"E"},
             {"N":"mag_err","F":"E"},
             {"N":"ra","F":"D"},
             {"N":"dec","F":"D"},
             {"N":"mlim","F":"E"},
             {"N":"fwhm","F":"E"},
             {"N":"subname","F":"60A"},
             {"N":"field","F":"10A"},
             {"N":"ccd","F":"I"},
             {"N":"type_auto","F":"10A"},
             {"N":"type_best","F":"10A"},
             {"N":"type_best_user","F":"10A"},
             {"N":"rbscore","F":"E"},
             ]
    colarray=[]
    for c in collist:
        colarray.append(pyfits.Column(name=c["N"],format=c["F"]))

    return pyfits.new_table(pyfits.ColDefs(colarray),nrows=nrows)
    
def fill_daily_table(tbhdu,canddata,candhead,ncand,first=0,subname=None):
    mjd=candhead['MJD-OBS']
    filter=candhead['FILTNAME']
    if subname is None: subname=candhead['SUBNAME']
    field=subname.split('_')[-4]
    ccd=int(subname.split('_')[-1])                             
    try:
        zeropoint=candhead['ZP']
    except:
        zeropoint=ZP[filter]
    try:
        mlim=f2m(candhead['subflim'],zeropoint)
    except:
        mlim=-99.
        
    for ind in range(ncand):
        index=ind+first
        tbhdu.data['mjd'][index]=mjd
        tbhdu.data['filter'][index]=filter
        mag,mag_err=f2m(canddata['f4sub'][ind],zeropoint,flux_err=canddata['df4sub'][ind])
        tbhdu.data['mag'][index]=mag
        tbhdu.data['mag_err'][index]=mag_err
        tbhdu.data['ra'][index]=canddata['rasub'][ind]
        tbhdu.data['dec'][index]=canddata['decsub'][ind]
        tbhdu.data['mlim'][index]=mlim
        tbhdu.data['fwhm'][index]=canddata['fwhmsub'][ind]
        tbhdu.data['subname'][index]=subname
        tbhdu.data['field'][index]=field
        tbhdu.data['ccd'][index]=ccd
        tbhdu.data['rbscore'][index]=canddata['rbscore'][ind]
        tbhdu.data['type_auto'][index]=' '
        tbhdu.data['type_best'][index]=' '
        tbhdu.data['type_best_user'][index]=' '
    return tbhdu,first+np.arange(ncand)

def update_daily_table(tbhdu,canddata,candhead,ncand,tbindex=None,subname=None):
    mjd=candhead['MJD-OBS']
    filter=candhead['FILTNAME']
    if tbindex is None:tbindex=range(ncand)
    if subname is None: subname=candhead['SUBNAME']
    field=subname.split('_')[-4]
    ccd=int(subname.split('_')[-1])                             
    try:
        zeropoint=candhead['ZP']
    except:
        zeropoint=ZP[filter]
    try:
        mlim=f2m(candhead['subflim'],zeropoint)
    except:
        mlim=-99.
        
    for ind in range(ncand):
        index=tbindex[ind]
        tbhdu.data['mjd'][index]=mjd
        tbhdu.data['filter'][index]=filter
        mag,mag_err=f2m(canddata['f4sub'][ind],zeropoint,flux_err=canddata['df4sub'][ind])
        tbhdu.data['mag'][index]=mag
        tbhdu.data['mag_err'][index]=mag_err
        tbhdu.data['ra'][index]=canddata['rasub'][ind]
        tbhdu.data['dec'][index]=canddata['decsub'][ind]
        tbhdu.data['mlim'][index]=mlim
        tbhdu.data['fwhm'][index]=canddata['fwhmsub'][ind]
        tbhdu.data['subname'][index]=subname
        tbhdu.data['field'][index]=field
        tbhdu.data['ccd'][index]=ccd
        tbhdu.data['rbscore'][index]=canddata['rbscore'][ind]

    return tbhdu,tbindex

def init_lookup(lookupfile):
    c1=pyfits.Column(name='name',format='20A')
    c2=pyfits.Column(name='ra',format='D')
    c3=pyfits.Column(name='dec',format='D')
    lookuptb=pyfits.new_table([c1,c2,c3],nrows=0)
    lookuptb.writeto(lookupfile,clobber=True)
    return
 

def init_cand_name(tbdata):
    collist1=[{"N":"name","F":"20A","I":' '},
              {"N":"ra_best","F":"D","I":0.},
              {"N":"dec_best","F":"D","I":0.},
              {"N":"n_night","F":"I","I":0},
              ]
    collist2=[{"N":"mjd","F":"PD"},
              {"N":"filter","F":"PA"},
              {"N":"mag","F":"PE"},
              {"N":"mag_err","F":"PE"},
              {"N":"ra","F":"PD"},
              {"N":"dec","F":"PD"},
              {"N":"mlim","F":"PE"},
              {"N":"fwhm","F":"PE"},
              {"N":"rbscore","F":"PE"},
              {"N":"subname","F":"PA"}
              ]
    collist3=[{"N":"type_auto","F":"PA"},
              {"N":"type_auto_comment","F":"PA"},
              {"N":"type_best","F":"PA"},
              {"N":"type_best_user","F":"PA"},
              {"N":"type_best_user_comment","F":"PA"},
              ]
    colarray=[]
    for c in collist1:
        colarray.append(pyfits.Column(name=c["N"],format=c["F"]))
    for c in collist2:
        if c["F"][1] in 'A':
            colarray.append(pyfits.Column(name=c["N"],format=c["F"],array=[list(tbdata[c["N"]][0])]))
        else:
            colarray.append(pyfits.Column(name=c["N"],format=c["F"],array=[tbdata[c["N"]]]))
    for c in collist3:
        colarray.append(pyfits.Column(name=c["N"],format=c["F"],array=[' ']))
    
    tbhdu=pyfits.new_table(pyfits.ColDefs(colarray))
    tbhdu.data[0]['ra_best']=tbhdu.data[0]['ra'][0]
    tbhdu.data[0]['dec_best']=tbhdu.data[0]['dec'][0]
    tbhdu.data[0]['name']=name_from_radec(tbhdu.data[0]['ra_best'],tbhdu.data[0]['dec_best'])
    tbhdu.data[0]['n_night']=1

    #classify auto type
    
    return tbhdu

def update_cand_name(tbhdu,tbdata):
    collist2=[{"N":"mjd","F":"PD"},
              {"N":"filter","F":"PA"},
              {"N":"mag","F":"PE"},
              {"N":"mag_err","F":"PE"},
              {"N":"ra","F":"PD"},
              {"N":"dec","F":"PD"},
              {"N":"mlim","F":"PE"},
              {"N":"fwhm","F":"PE"},
              {"N":"rbscore","F":"PE"},
              {"N":"subname","F":"PA"}
              ]
    for c in collist2:
        if c["F"][1] in 'A':
            tbhdu.data[0][c["N"]]=np.append(tbhdu.data[0][c["N"]],list(tbdata[c["N"]][0]))
        else:
            tbhdu.data[0][c["N"]]=np.append(tbhdu.data[0][c["N"]],tbdata[c["N"]])

    #=====================================
    #check for repeated entries here
    #no same mjd?
    #=====================================
    
    #update ra_best, dec_best, use weighted mean
    weights=2.5/np.log(10)/tbhdu.data[0]['mag_err']
    tbhdu.data[0]['ra_best']=sum(tbhdu.data[0]['ra']*weights)/sum(weights)
    tbhdu.data[0]['dec_best']=sum(tbhdu.data[0]['dec']*weights)/sum(weights)
    #update n_night
    #same night always on same UT day
    lastmjd=np.floor(tbhdu.data[0]['mjd'][-1])
    allmjd=np.floor(tbhdu.data[0]['mjd'][0:-1])
    if any(abs(lastmjd-allmjd) > 0.5): tbhdu.data[0]['n_night']=tbhdu.data[0]['n_night']+1

    #classify auto type
    
    tbhdu.update()
    return tbhdu

def plot_cand_lc(mjd,mag,mag_err,mlim,filter,output):
    #mjd, mag, mag_err, mlim, filter should have same dimension
    nmjd=len(mjd)
    nmag=len(mag)
    if nmjd != nmag:
        return -1
    nmagerr=len(mag_err)
    if nmagerr !=nmag:
        return -1
    nfilter=len(filter)
    if nfilter !=nmag:
        return -1
    nmlim=len(mlim)
    if nmlim !=nmag:
        return -1
        
    filterset=['u','v','g','r','i','z'] #filters to go through
    filterclr=['c','b','g','r','m','k']
    filtersym=['p','s','*','o','d','x']


    mjd=mjd-50000
    
    fig=pyplot.figure(figsize=(8,5))
    ax1=pyplot.axes([0.1,0.15,0.75,0.75])
    
    #figure out mjd range
    mjd1=min(mjd)
    mjd2=max(mjd)
    mjdr=mjd2-mjd1
    mjdpad=mjdr*0.07
    if mjdpad <0.5 : mjdpad=0.5
    
    #figure out mag range
    maglb=min(mag-mag_err)
    maghb=max(mag+mag_err)
    magr=maghb-maglb
    magpad=magr*0.07
    mag1=maglb-magpad
    mag2=maghb+magpad
    for ind in range(nmlim):
        if mlim[ind]< mag1 and mlim[ind]>0:
            mlim[ind]=mag1
        elif mlim[ind]> mag2 and mlim[ind]>0:
            mag2=mlim[ind]
    xlim=(mjd1-mjdpad,mjd2+mjdpad)
    ylim=(mag2,mag1)

    for filtind,filt in enumerate(filterset):
        ind=[testind for testind,testfilt in enumerate(filter) if testfilt in filt]
        if len(ind)>0:
            plot=pyplot.errorbar(mjd[ind],mag[ind],yerr=mag_err[ind],fmt=filterclr[filtind]+filtersym[filtind],label=filt)
            plot=pyplot.plot(mjd[ind],mlim[ind],filterclr[filtind]+'v')
    
    xlim=ax1.set_xlim(xlim)
    ylim=ax1.set_ylim(ylim)
    xticks=ax1.get_xticks()
    xticks=list(set(np.floor(xticks*10.)/10.))
    xticks.sort()

    legend=pyplot.legend(loc=(1.05,0.3))

    ax2=pyplot.twiny()
    xlim=ax2.set_xlim(ax1.get_xlim())
    tmp=ax2.set_xticks(xticks)
    xticklabels=[]
    for xtick in xticks:
        xticklabels.append('%.1f'%xtick)
    tmp=ax2.set_xticklabels(xticklabels)
    tmp=ax2.set_xlabel('MJD-50000')

    #set ax1 label to date
    xintticks=list(set(np.floor(xticks)))
    xintticks.sort()
    datestr=[]
    for xinttick in xintticks:
        date=mjd2datetime(xinttick+50000)
        datestr.append('%s-%s-%s'%(date.day,date.month,date.year))
    tmp=ax1.set_xticks(xintticks)
    tmp=ax1.set_xticklabels(datestr,rotation=10)
    xlim=ax1.set_xlim(xlim)
    tmp=ax1.set_xlabel('UT Date')
    tmp=ax1.set_ylabel('Mag')

    pyplot.savefig(output)
    return


def f2m(flux,zp,flux_err=None):
    mag = -2.5*np.log10(flux)+zp
    if not flux_err is None:
        mag_err = 2.5*flux_err/flux/np.log(10)
        return mag,mag_err
    else:
        return mag

def name_from_radec(ra,dec):
    #ra, dec in degree
    name='smt_'+radec2str(ra,dec)
    return name

def radec2str(rad,decd):
    ra=rad/15.
    rastr=sixtystr(ra,decimal=2,delim='').replace('.','')
    if decd<0: decsign='-'
    else: decsign='+'
    dec=abs(decd)
    decstr=sixtystr(dec,decimal=1,delim='').replace('.','')
    return rastr+decsign+decstr

def sixtystr(x,decimal=0,delim=':'):
    xh=np.floor(x)
    xm=np.floor((x-xh)*60.)
    xs=(x-xh-xm/60.)*3600.
    xs=np.round(xs*10**decimal)/10**decimal
    
    if decimal==0:
        xsstr='%02f'%xs
    else:
        xsstr=eval('\'%0'+'%d.%df'%(3+decimal,decimal)+'\'%xs')
        
    return delim.join(('%02d'%xh,'%02d'%xm,xsstr))

def sixty(deg):
    adeg=abs(deg)
    dd=np.floor(adeg)
    mm=np.floor((adeg-dd)*60.)
    ss=((adeg-dd)*60.-mm)*60.
    if deg<0: dd= -1. *dd
    return (dd,mm,ss)

def get_last_part(alist):
    joined=''.join(alist)
    return joined[joined.rfind(' ')+1:]
    
