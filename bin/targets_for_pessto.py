#!/usr/bin/env python

import os
import argparse
import pyfits
import numpy as np
from Utils.Constants import FilenamesXsient, FilenamesSub, ThumbNames, PipelinePath
from Utils.Record import AsciiRecord, KeywordRecord
from Utils.Catalog import CatalogEntry

from mydjango.settings import MEDIA_URL
from django.core.urlresolvers import reverse
import mydjango.followup.models as fu

import datetime

def main():
    parser = argparse.ArgumentParser(
        description='Write SkyMapper Transient candidates csv file for PESSTO')
    parser.add_argument('--maxage', type=float, default=30,
                        help='maximum age of the earliest identification of a transient candidate')
    parser.add_argument('--output', default=PipelinePath.django_media+'/trigger/pessto.csv',
                        help='full path of the output csv file')
    

    args = parser.parse_args()
    maxage=args.maxage
    output=args.output
    
    now=datetime.datetime.now()
    cands=fu.Transient.objects.filter(type__type__contains='Cand')
    

    targets=[]
    for cand in cands:
        if cand.transient_change_set.filter(type__type__contains='Cand').order_by('time').values_list('time',flat=True)[0] < now -datetime.timedelta(maxage):
            #    if cand.transient_change_set.latest('time').time < now -datetime.timedelta(maxage):
            continue
        target=PesstoTarget(cand)
        targets.append(target)
        
    AsciiRecord.write_ascii_file(output,targets,header=True)


def getLC(trans):
    #get lightcurve
    xsientnames = FilenamesXsient(trans.name,trans.field.id,trans.ccd)
    lcfile=xsientnames.abssublcfname
    lctable=pyfits.getdata(lcfile)
    return lctable

def getDet(trans):
    lctable=getLC(trans)
    #get detections
    detind=np.where(lctable['MAG_CAL_ERR'] >=0)[0]
    if len(detind)==0:
        return None
    return lctable[detind]

def firstDetect(trans):
    det=getDet(trans)
    if len(det)==0:
        return None
    return det[0]

def lastDetect(trans):
    det=getDet(trans)
    if len(det)==0:
        return None
    return det[-1]

def lastLimit(trans):
    lctable=getLC(trans)
    detind=np.where(lctable['MAG_CAL_ERR'] >=0)[0]
    #get limits
    limind=np.where(lctable['MAG_CAL_ERR'] < 0)[0]
    if len(limind)==0:
        return None
    limind2=np.where(lctable[limind]['JD_OBS']<lctable[detind][0]['JD_OBS'])[0]
    if len(limind2)==0:
        return None
    return lctable[limind][limind2][-1]
    
def thumbURLs(trans):
    lcrow=lastDetect(trans)
    if not lcrow:
        return [None,None,None]
    subid=FilenamesSub.subid_from_fields(
        runtag=lcrow['RUNTAG'],field=int(lcrow['SM_FIELD']),filter=lcrow['FILTER'],
        obsseq=lcrow['OBSSEQ'], subfield=int(lcrow['SUBFIELD']))
    thumbs=ThumbNames(trans.name,subid)
    return [MEDIA_URL+'thumbs/%s/'%trans.name+thumbs.new,MEDIA_URL+'thumbs/%s/'%trans.name+thumbs.ref,MEDIA_URL+'thumbs/%s/'%trans.name+thumbs.diff]

def lastComment(trans):
    log=trans.transient_change_set.latest('time')
    comment='{0}: typed as {1}; {2}'.format(log.time,log.type.type,log.comment)
    return comment

def cleanurl(x):
    return os.path.normpath(x).replace('http:/','http://')
    

class PesstoTarget(AsciiRecord,KeywordRecord,CatalogEntry):
    """
    csv file for PESSTO marshall retrieval.
    """
    
    _ascii_separator = " | "
    _ascii_fields = [
        ['candidateID',     0, lambda x: x.strip(),'{0:22s}',   '-',  'lambda x: x.name', lambda x: x],
        ['RA',              1, lambda x:  float(x),'{0:10.6f}', -99., 'lambda x:   x.ra', lambda x: x],
        ['DEC',             2, lambda x:  float(x),'{0:9.5f}' , -99., 'lambda x:  x.dec', lambda x: x],       
        ['mag',             3, lambda x:  float(x),'{0:5.2f}' , -99., 'lambda x:  getDet(x)', lambda x: x[-1]['MAG_CAL']],
        ['magerr',          4, lambda x:  float(x),'{0:6.2f}' , -99., 'lambda x:  getDet(x)', lambda x: x[-1]['MAG_CAL_ERR']],
        ['mjd',             5, lambda x:  float(x),'{0:9.3f}' , -99., 'lambda x:  getDet(x)', lambda x: x[-1]['JD_OBS']-2400000.5],
        ['filt',            6, lambda x: x.strip(),'{0:4s}' ,   '-', 'lambda x:  getDet(x)', lambda x: x[-1]['FILTER']],
        ['discMag',         7, lambda x:  float(x),'{0:7.2f}' , -99., 'lambda x:  getDet(x)', lambda x: x[0]['MAG_CAL']],
        ['discMJD',         8, lambda x:  float(x),'{0:9.3f}' , -99., 'lambda x:  getDet(x)', lambda x: x[0]['JD_OBS']-2400000.5],
        ['discFilt',        9, lambda x: x.strip(),'{0:8s}' ,   '-',  'lambda x:  getDet(x)', lambda x: x[0]['FILTER']],
        ['bestType',       10, lambda x: x.strip(),'{0:10s}' ,  '-',  'lambda x:  x.type.type', lambda x: x],
#        ['catalogType',     7, lambda x:   float(x),'{0:9.5f}' , '', lambda x:  x],
#        ['hostZ',           8, lambda x:   float(x),'{0:5.2f}' , '', lambda x:  x],
        ['noneMag',        11, lambda x:  float(x),'{0:7.2f}' , -99., 'lambda x:  lastLimit(x)', lambda x: x['MAG_CAL']],
        ['noneMJD',        12, lambda x:  float(x),'{0:9.3f}' , -99., 'lambda x:  lastLimit(x)', lambda x: x['JD_OBS']-2400000.5],
        ['noneFilt',       13, lambda x: x.strip(),'{0:8s}' ,   '-',  'lambda x:  lastLimit(x)', lambda x: x['FILTER']],
        ['newThumbURL',    14, lambda x: x.strip(),'{0:160s}' , '-',  'lambda x:  thumbURLs(x)', lambda x: cleanurl(x[0])],
        ['refThumbURL',    15, lambda x: x.strip(),'{0:160s}' , '-',  'lambda x:  thumbURLs(x)', lambda x: cleanurl(x[1])],
        ['diffThumbURL',   16, lambda x: x.strip(),'{0:160s}' , '-',  'lambda x:  thumbURLs(x)', lambda x: cleanurl(x[2])],
        ['finderURL',      17, lambda x: x.strip(),'{0:100s}' , '-',  'lambda x:  FilenamesXsient(x.name,x.field.id,x.ccd).finderfname', lambda x: cleanurl(x.replace(PipelinePath.django_media,MEDIA_URL))],
        ['candidateURL',   18, lambda x: x.strip(),'{0:80s}' ,  '-',  'lambda x:  "http://www.mso.anu.edu.au/skymapper/smt/"+reverse("trans_info",args=(x.name,))', lambda x: cleanurl(x)],
        ['status',         19, lambda x: x.strip(),'{0:15s}' ,  '-',  'lambda x:  x.follow_up_status.status', lambda x: x],
        ['numDet',         20, lambda x:    int(x),'{0:6d}' ,    0,   'lambda x:  x.n_det', lambda x: x],
        ['comment',        21, lambda x: x.strip(),'{0:100s}' , '-',  'lambda x:  lastComment(x)', lambda x: x],
#        ['numSpectra',      14, lambda x:   float(),'{0:9.5f}' , '', lambda x:  x],
#        ['specType',        15, lambda x:   float(),'{0:9.5f}' , '', lambda x:  x],
#        ['discPhase',       16, lambda x:   float(),'{0:9.5f}' , '', lambda x:  x],
#        ['minObsDate',      17, lambda x:   float(),'{0:9.5f}' , '', lambda x:  x],
#        ['Z',               18, lambda x:   float(),'{0:9.5f}' , '', lambda x:  x],
#        ['subtractedLightcurveURL',19, lambda x:   float(),'{0:9.5f}' , '', lambda x:  x],
#        ['isFollowed',             20, lambda x:   float(),'{0:9.5f}' , 'No', lambda x:  x],
        ]
        
    _keyword_fields = [ (f[0],f[4]) for f in _ascii_fields ]

    
    def _init_from_transient(self,trans):
        #init from SMT transient record
        lastfval=None
        for fid,f in enumerate(self._ascii_fields):
            if f[5] == self._ascii_fields[fid-1][5]:
                try:
                    val=f[2](f[6](lastfval))
                except:
                    val=f[4]
            else:
                try:
                    fval=eval(f[5])(trans)
                    lastfval=fval
                    val=f[2](f[6](fval))
                except:
                    val=f[4]
            if not val: val=f[4]
            setattr(self,f[0],val)

    def __init__(self,trans):
        self._init_from_transient(trans)



if __name__ == "__main__":
    main()


