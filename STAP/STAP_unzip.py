#! /usr/bin/env python

import argparse
import pyfits
import numpy as np
import subprocess
import shlex
import os

def unzip(imgname, timeout=None):
    """Unzip a gzipped file

    """
    if not imgname[-3:] == '.gz':
        print imgname + ' is not a gzipped file'
        return(1)
    else:
        subprocess.call(['gunzip',imgname])
        return

    
def main():
    """Wrapper allowing STAP.unzip() to be run from the command line"""
    parser = argparse.ArgumentParser(description=
                                     'Unzip a gzipped file')
    parser.add_argument('imgname', help='filename of input image')
    parser.add_argument('--timeout',type=int, default=15,
                        help='maximum running time allowed (default: %(default)s)')
    args = parser.parse_args()
    unzip(args.imgname, timeout=args.timeout)

if __name__ == "__main__":
    main()
