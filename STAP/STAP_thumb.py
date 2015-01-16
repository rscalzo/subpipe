#! /usr/bin/env python
"""
Solve WCS of input image and update the header

Syntax: STAP_WCS.py imname [-timeout TIMEOUT] [-poll POLL]
    
Inputs:
    imname: filename of input image including full path

Optional inputs:
    timeout: maximum running time (in seconds) allowed for external call
             default is 300s
    poll: polling interval (in seconds) for external call
             default is 10s

Description:
    Solve for the WCS of input image using either solve_field or other method
    Option to change timeout and polling interval

Exit status:
    0, successful
    1, see specific returned error message
    2, syntax error

Specifications:
    External system call:
        /usr/local/astrometry/bin/solve-field
        or else
    Memory requriements: image + sourcelist + catalog 
    Scratch disk space requirements: none  
    Typical wall clock time needed: 45 s for solve-field
    
    Possible requirments: 
        config for source extraction
    
"""

import argparse
from STAP_comm import STAP_callexternal
import sys

##Parse input arguments##
parser = argparse.ArgumentParser(description='Solve WCS of input image and update the header')
parser.add_argument('imname', 
                    help='filename of input image')
parser.add_argument('-timeout',type=float, default=600,
                    help='maximum running time allowed for external call (default: %(default)ss)')
parser.add_argument('-poll',type=float, default=10,
                    help='polling interval for external call (default: %(default)ss)')
# RS 2011/04/28:  added log file option to pass to STAP_callexternal
parser.add_argument('-log',default=None,
                    help='optional log file (default: write to stdout)')
# RS 2012/02/23:  Disabled thumbnail generation for Bogus detections.
# Added -training flag to re-enable for training set generation only.
parser.add_argument('-training',action='store_true',default=False,
                    help='generate thumbnails for Bogus detections as well as Reals')
args = parser.parse_args()

##-----------!!------------##
##Define other variables here##

##DEFINE EXTERNAL COMMAND STRING HERE##
cmd_name="make_thumbs.py"
cmd_ext="%s %s"%(cmd_name,args.imname)
if args.training:  cmd_ext += " --training"
##-------------------------##

##Start external call ##
[status,message,stdoutstr]=STAP_callexternal(cmd_ext,timeout=args.timeout,poll=args.poll,log=args.log)
if status !=0 :
    sys.exit(message)
