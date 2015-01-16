#! /usr/bin/env python
"""
Apply flat fielding and/or fringe correction 

Syntax: STAP_imcorr.py imname outname [-ff FFNAME] [-fr FRNAME] [-noampjoin]
[-timeout TIMEOUT]
    
Inputs:
    imname: filename of input image including full path
    outname: output filename of corrected image

Optional inputs:
    ffname: filename of flat field including full path
    frname: filename of fringe pattern
    noampjoin: if flag set, no alignment between left and right halves of CCD
    timeout: maximum running time (in seconds) allowed for each step
             default is 20s, with polling interval fixed to 1 s    

Description:
    Apply flat fielding on input image if a flat field is supplied.
    Apply fringe correction on input image if a fringe pattern is supplied.
    Align the right half of CCD to the left, unless noampjoin flag is set.
    Write the corrected image to outname.
    Option to change timeout.

Exit status:
    0, successful
    1, see specific returned error message
    2, syntax error

Specifications:
    External system call:
        imarith
    Python function:
        ampjoin        
    Memory requriements: 32M x 2 + bits
    Scratch disk space requirements: none  
    Typical wall clock time needed: 3 s
    Config files needed: none
    Enviroment variables: none if imarith is in the path
"""

import argparse
import sys
import pyfits

parser = argparse.ArgumentParser(description='Apply flat fielding and/or fringe correction.')
parser.add_argument('imname', 
                    help='filename of input image')
parser.add_argument('outname',
                    help='filename of output image')
parser.add_argument('-ff',dest='ffname',
                    help='filename of flat field')
parser.add_argument('-fr',dest='frname',
                    help='filename of fring pattern')
parser.add_argument('-noampjoin',action="store_const",const=True,default=False,
                    help='if set, two halves of CCD are not aligned')
parser.add_argument('-timeout',type=int, default=60,
                    help='maximum running time allowed for each correction (default: %(default)s)')
# RS 2011/04/28:  added log file option to pass to STAP_callexternal
parser.add_argument('-log',default=None,
                    help='optional log file (default: write to stdout)')
args = parser.parse_args()

if args.noampjoin is True and args.ffname is None and args.frname is None:
    sys.exit("No action to be done on the input image")


def ampjoin():
    pass
