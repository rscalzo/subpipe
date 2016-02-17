#! /usr/bin/env python

import argparse
import os
import re
import ephem
import glob
import pyfits
from STAP_comm import print_cmd_line
from STAP.STAP_display import STAP_make_thumbs_xy, sky2xy
from Utils.Constants import FilenamesSub, Find_Image_Set, ThumbNames, follow_types
import subprocess

def makethumbs(lcfnames,xref=None, outputdir='.',localonly=False,redoall=False):
    #the initial call
    print_cmd_line("STAP_makethumbs.py", lcfnames, xref=xref,outputdir=outputdir,localonly=localonly)
    if not lcfnames:
        if not xref is None: 
            if not os.path.exists(xref): xref=None
        if xref:
            # FY - get the objectnames from xref file
            lcfnames=[]
            head=pyfits.getheader(xref)
            for (key, value) in head.items():
                namematch = re.match("NAME(\d{4})", key)     
                if namematch:
                    id=int(namematch.group(1))
                    objtype=head["TYPE{0:04d}".format(id)]
                    if objtype == "?" or objtype in follow_types:
                        if os.path.exists('./%s_sub_lc.fits'%value):
                            lcfnames.append('./%s_sub_lc.fits'%value)
        else:
            # RS 2013/10/05:  For now let's only get the subtracted light curves.
            # I think it would be useful to do the unsubtracted ones, but let's
            # implement that when we have a little more free time.
            lcfnames=glob.glob('./*_sub_lc.fits')
    if not lcfnames:
        print "No object (of unknown or followup type) that needs thumbnails."
        return

    # RS 2014/02/19:  Don't forget to print the command line so that we can
    # reproduce what was done in a debugging context.
    print_cmd_line("STAP_makethumbs.py", lcfnames, outputdir='.',
                   localonly=False, redoall=False)

    #make thumbs of _new _ref and _diff for all thumbs
    lc_registry = extract_lcfiles(lcfnames,localonly=localonly,redoall=redoall)

    for subid in lc_registry:
        allxs, allys, allnames=merge_xy_radec(lc_registry[subid]["imnames"]["newname"],lc_registry[subid]["xs"],lc_registry[subid]["ys"],
                                              lc_registry[subid]["ras"],lc_registry[subid]["decs"],lc_registry[subid]["thumbnew"],lc_registry[subid]["limnew"])
        STAP_make_thumbs_xy(lc_registry[subid]["imnames"]["newname"],allxs,allys,outdir=outputdir,outputnames=allnames)

        allxs, allys, allnames=merge_xy_radec(lc_registry[subid]["imnames"]["refname"],lc_registry[subid]["xs"],lc_registry[subid]["ys"],
                                              lc_registry[subid]["ras"],lc_registry[subid]["decs"],lc_registry[subid]["thumbref"],lc_registry[subid]["limref"])
        STAP_make_thumbs_xy(lc_registry[subid]["imnames"]["refname"],allxs,allys,outdir=outputdir,outputnames=allnames)

        allxs, allys, allnames=merge_xy_radec(lc_registry[subid]["imnames"]["subname"],lc_registry[subid]["xs"],lc_registry[subid]["ys"],
                                              lc_registry[subid]["ras"],lc_registry[subid]["decs"],lc_registry[subid]["thumbdiff"],lc_registry[subid]["limdiff"])
        STAP_make_thumbs_xy(lc_registry[subid]["imnames"]["subname"],allxs,allys,outdir=outputdir,outputnames=allnames)

    return 0

def extract_lcfiles(lcfnames,localonly,redoall):
    badsubid = []
    lc_registry = {}
    for lcfile in lcfnames:
        objname=pyfits.getheader(lcfile)['CANDNAME']
        tb=pyfits.getdata(lcfile)
        hasobsid='OBSID' in tb.dtype.names
        for rowid, row in enumerate(tb):
            if hasobsid:
                subid=FilenamesSub.subid_from_fields(
                    runtag=row['RUNTAG'],field=int(row['SM_FIELD']),filter=row['FILTER'],
                    obsid=row['OBSID'],obsseq=row['OBSSEQ'], subfield=int(row['SUBFIELD']))
            else:
                subid=FilenamesSub.subid_from_fields(
                    runtag=row['RUNTAG'],field=int(row['SM_FIELD']),filter=row['FILTER'],
                    obsseq=row['OBSSEQ'], subfield=int(row['SUBFIELD']))
            if (not subid in lc_registry) and (not subid in badsubid):
                imnames = Find_Image_Set(subid,localonly=localonly)
                if imnames:
                    lc_registry[subid]={"imnames":imnames,"objnames":[],"xs":[],"ys":[],"limnames":[],"ras":[],"decs":[],
                                        "thumbnew":[],"thumbref":[],"thumbdiff":[],
                                        "limnew":[],"limref":[],"limdiff":[]}
                else:
                    badsubid.append(subid)
            if subid in lc_registry:
                allthumbnames=ThumbNames(objname,subid,find=(not redoall))
                if not allthumbnames.exists:
                    if row['X_IMAGE']>0 and row['Y_IMAGE']>0:
                        lc_registry[subid]["thumbnew"].append(allthumbnames.new)
                        lc_registry[subid]["thumbref"].append(allthumbnames.ref)
                        lc_registry[subid]["thumbdiff"].append(allthumbnames.diff)
                        lc_registry[subid]["xs"].append(row['X_IMAGE'])
                        lc_registry[subid]["ys"].append(row['Y_IMAGE'])
                    else:
                        lc_registry[subid]["limnew"].append(allthumbnames.new)
                        lc_registry[subid]["limref"].append(allthumbnames.ref)
                        lc_registry[subid]["limdiff"].append(allthumbnames.diff)
                        lc_registry[subid]["ras"].append(row['ALPHA_J2000'])
                        lc_registry[subid]["decs"].append(row['DELTA_J2000'])
                else:
                    pass
                    #print "Not re-making thumbnails for object: {0}, subid: {1}.".format(objname,subid)

    return lc_registry


def merge_xy_radec(imname,xs,ys,ras,decs,objnames,limnames):
    if not ras:
        return xs,ys,objnames
    
    isgzip = False
    if '.gz' in imname[-3:]:
        isgzip = True
    
    newimname=imname
    if isgzip:
        subprocess.call(['gunzip',imname])
        newimname=imname[:-3]
        
    [lxs,lys]=sky2xy(newimname,ras,decs)
    if isgzip:
        subprocess.call(['gzip',newimname])

    return xs+lxs, ys+lys, objnames+limnames


def main():
    parser = argparse.ArgumentParser(description="Make thumbnails for the candidates")
    parser.add_argument("--lcfnames", nargs="+",
                        help='lightcurve files for candidates',default=None)
#    parser.add_argument("--newname",help='new image name')
#    parser.add_argument("--refname",help='ref image name')
#    parser.add_argument("--subname",help='sub image name')
    parser.add_argument("--xref", default=None,help="optional input xref file")
    parser.add_argument("--localonly",action='store_true', default=False)
    parser.add_argument("--output",help='output directory name')
    args = parser.parse_args()
    
    
    makethumbs(args.lcfnames,xref=args.xref,outputdir=args.output,localonly=args.localonly)

if __name__ == "__main__":
    main()

