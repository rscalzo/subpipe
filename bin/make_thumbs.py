#!/usr/bin/env python

# ============================================================================
# FY 2011/09/28: Make thumbnails in scratch disk
# ============================================================================

import os, argparse
import pyfits
from STAP.STAP_display import STAP_make_thumbs_xy, STAP_make_thumbs
from Utils.RealBogus import RealBogus

# Parse the arguments
# RS 2012/02/23:  Disabled thumbnail generation for Bogus detections.
#     Added --training flag to re-enable for training.
parser = argparse.ArgumentParser \
         (description='Make png thumbnails for all candidates.')
parser.add_argument('diffname', help='Difference image filename')
parser.add_argument('--training', action='store_true', default=False,
                    help='Write Bogus thumbs out to disk')
args = parser.parse_args()

subname=args.diffname
candfile=subname+'.cands'
if not os.path.exists(candfile):
    exit(0)
imname=os.path.basename(subname).replace('.fits','')
subpath=os.path.dirname(subname)
if subpath == '':  subpath = "."

#find new and reference files for the subtraction
subhdr=pyfits.getheader(subname)
inim=subhdr['TARGET']
tmplim=subhdr['TEMPLATE']
inim=subpath+'/'+os.path.basename(inim)
tmplim=subpath+'/'+os.path.basename(tmplim)

#get the list of candidates
candtb=pyfits.getdata(candfile)
if candtb is None:
    exit(0)

# RS 2012/02/20:  We always have x and y for our candidates, and we don't
# need RA and DEC to know where to take out the thumbnails.  So I changed
# this to useradec = 0.
useradec=0

# RS:  Check which things are Bogus and which are Real.
if args.training:  extractidx = (candtb['rbscore'] > -1)
else:  extractidx = (candtb['rbscore'] > RealBogus.rbscore_thresh)

if useradec ==1:
    ras=candtb['rasub'][extractidx]
    decs=candtb['decsub'][extractidx]
    #make thumbs out of sub image
    STAP_make_thumbs(subname,ras,decs,box=1.0,outdir=subpath,usesubid=imname,suffix='_diff')
    #make thumbs out of  new image
    STAP_make_thumbs(inim,ras,decs,box=1.0,outdir=subpath,usesubid=imname,suffix='_new')
    #make thumbs out of  ref image
    STAP_make_thumbs(tmplim,ras,decs,box=1.0,outdir=subpath,usesubid=imname,suffix='_ref')
    
else:    
    xs=candtb['xsub'][extractidx]
    ys=candtb['ysub'][extractidx]
    #make thumbs out of sub image
    STAP_make_thumbs_xy(subname,xs,ys,box=120,outdir=subpath,usesubid=imname,suffix='_diff')
    #make thumbs out of  new image
    STAP_make_thumbs_xy(inim,xs,ys,box=120,outdir=subpath,usesubid=imname,suffix='_new')
    #make thumbs out of  ref image
    STAP_make_thumbs_xy(tmplim,xs,ys,box=120,outdir=subpath,usesubid=imname,suffix='_ref')
