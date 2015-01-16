#!/usr/bin/env python

#SM-WCS - put WCS on a SkyMapper MOSAIC image  
#uses Astronet on amp 33 to find pointing of telescop
#puts a rough WCS on each chip based on ASTROnet solution
#matches stars on each chip against ucac3
#fits a ZPN + distortion matrix to a grid
#then puts a TNX-5order onto each image
#this provides the best astrometry possible within
#the confines of the WCS system that can be
#easily interpreted by Sextractor, SWARP, DS9 etc.
#Based on CASU Code written by Mike Irwin and Jim Lewis        
#Brian Schmidt, Australian National University April 2006 in PERL     
#moved to Python May 2011
#
#21Jul11 - added TNX to be calculated from same CRVAL at center of array
#21Jul11 - made fit go from 2mass rather than ucac2 - makes it slow, but better
#

import argparse
import numpy
import pyfits
import os
import sys
import stat
import subprocess
import ephem
#
#subroutines
#

# RS 2011/11/02:  We already had to create an environment variable $BRIANWCS
# in which to store global configuration files (in $BRIANWCS/etc/).  We might
# as well not clutter up $PATH with the executables in $BRIANWCS/bin/.

brianwcspath = os.environ["BRIANWCS"]
if not os.path.exists(brianwcspath):
    print "SM-WCS.py:  Whoops!  Please setenv BRIANWCS and try again."
    exit(-1)

#
# write header of MEF with rough WCS based on RA-CENT,DEC-CENT or 
#

def roughwcs(raval,decval,header,tan,projp3):
    import sys
    import math

    
    if header.get('IMAGEID')==None:
        sys.stderr.write('SM-WCS.py:  NO IMAGEID in image header\n')
        exit(-1)
    else:
        imageid=header.get('IMAGEID')


    if header.get('ROTSKYPA')==None:
        sys.stderr.write('SM-WCS.py:  NO ROTSKYPA in image header - setting to 0\n')
        rotskypa=0
    else:
        rotskypa=header.get('ROTSKYPA')    

    if header.get('NAXIS1')==None:
        sys.stderr.write('SM-WCS.py:  NO NAXIS1 in image header\n')
        exit(-1)
    else:
        naxis1=header.get('NAXIS1')


    #which amp do we really want to use...important if images are joined
        
    if naxis1 == 2048: #image is wrong
        chip=imageid*2
        if chip > 32:
            chip=chip-1 # way images are joined right now
                        # top half of array, 1,1 is in lower number chip
    else: 
        sys.stderr.write('SM-WCS.py:  image is not valid size, expecting naxis1=2048\n')
        exit(-1)

        
    scale=0.497/3600 # SkyMapper scale - hardwired, could put in configuration

    
    rotskypa+=90
    rotskypa=rotskypa*math.pi/180.
    flipx=1
    flipy=-1
    ra=180/3.1415926535*ephem.hours(raval) # converts to degrees
    dec=180/3.1415926535*ephem.degrees(decval) #converts to degrees

    CRPIX1=crpix1[chip]
    CRPIX2=crpix2[chip]


    tpa=ra
    tpd=dec

    print tpa,tpd,raval,decval

    aap=0*scale
    bbp=1*scale*flipx
    ddp=1*scale*flipy
    eep=0*scale

    #correct for rotation
    ca=math.cos(rotskypa)
    sa=math.sin(rotskypa)

    cd1_1 = aap*ca + bbp*sa
    cd1_2 = -aap*sa + bbp*ca
    cd2_1 = ddp*ca + eep*sa
    cd2_2 = -ddp*sa + eep*ca
    

    header.update('CRVAL1',tpa)
    header.update('CRVAL2',tpd)
    header.update('CRPIX1',CRPIX1)
    header.update('CRPIX2',CRPIX2)
    header.update('CD1_1',cd1_1)
    header.update('CD2_1',cd2_1)
    header.update('CD1_2',cd1_2)
    header.update('CD2_2',cd2_2)
    

    if tan == True:
        header.update('CTYPE1','RA---TAN')
        header.update('CTYPE2','DEC--TAN')
    else:
        header.update('CTYPE1','RA---ZPN')
        header.update('CTYPE2','DEC--ZPN')
        header.update('PROJP1',1.0)
        header.update('PROJP3',projp3)
        header.update('PV2_0',0.0)
        header.update('PV2_1',1.0)
        header.update('PV2_2',0.0)
        header.update('PV2_3',projp3)        
    
    #NOTE these need to be committed when done using flush

    
#
# return crval for an image for a given amp using WCS on it.
#
def getcrval(amp,image):
# getcrval takes an image with a WCS on it, use crpix values to find CRVAL
# for the array. Useful to find central RA and DEC of array

# make system call to xy2sky
# xy2sky image.fits CRPIX1 CRPIX2

  system= 'xy2sky '+ image + ' ' + str(crpix1[amp]) + ' ' + str(crpix2[amp]) 
  print system
  sky2xy = os.popen(system)
  line=sky2xy.readline()
  racent=line.split()[0]
  deccent=line.split()[1]
  return(racent,deccent)

############################################################
# MAIN PROGRAM
#

startingdir=os.getcwd()

configurefile = brianwcspath + '/etc/configure.wcs.skymapper'
if not os.path.exists (configurefile):
    print "SM-WCS.py:  Can't find configuration file ", configurefile, "!"
    exit(-1)
crval1 = -99999
crval2 = -99999
tan = False
p3 = 3.0
spliced = 0
mosaic = 0
astronet_init = True
noresid = False
ucac2 = False

parser = argparse.ArgumentParser(description='Take SkyMapper FITS files and put rough WCS on it')
parser.add_argument('imname', help='input filename (files assumed to be in N/filename_N.fits)')
parser.add_argument('-p3', help='value to put for ZPN 3rd order coef')
parser.add_argument('-tan', help='do a TAN projection, not ZPN',action='store_true')
parser.add_argument('-ucac2', help='use ucac2 rather than 2mass',action='store_true')
parser.add_argument('-noresid', help='straight zpn rather than zpn and distortions',action='store_true')
parser.add_argument('-crval', help='directly provide crval1 and crval2 in degrees',nargs=2)
parser.add_argument('-astronet', help="use astrometry.net for each chip instead of Brian's code",nargs=2)

args=parser.parse_args()
filename=args.imname
tan=args.tan
noresid=args.noresid
ucac2=args.ucac2

#
# readin CRPIX from configure.wcs.skymapper
#
crpix1,crpix2=numpy.loadtxt(configurefile, usecols=(1,2), unpack=True)

  
#
# Do we want to run astronet to get the position of the chip
# or the center of the array if a MEF
# If so, make a temporary copy of file
#


if astronet_init==True:   # copy correct extension to tmp_filename
    os.chdir('17')
    useext=os.path.splitext(filename)[0] + '_17.fits'
    
    hdulist = pyfits.open(useext)
    prihdr = hdulist[0].header
    if prihdr.get('NEXTEND')!= None:
        print >> sys.stderr, \
                'SM-WCS.py:  I don\'t do multi-extension FITS files, bye'
        exit(-1)
        
    ra=prihdr['RA']
    dec=prihdr['DEC']
    hdulist.close()

    system='solve-field -T --parity neg --no-plots --no-fits2fits '+  useext + ' --scale-low .493 --scale-high 0.501 --scale-units app --ra ' + ra + ' --dec ' + dec + ' --radius 2 --overwrite --use-sextractor --sextractor-path /usr/local/bin/sex -N ' + os.path.splitext(useext)[0] + '_solved.fits'
    print system

    os.system(system + '>> /dev/null')

    if os.path.exists(os.path.splitext(useext)[0] + '.solved')==True:
        (RACENT,DECCENT)=getcrval(33,os.path.splitext(useext)[0] + '_solved.fits')
    else:
        sys.stderr.write('SM-WCS.py:  ASTRONET failed on image '+ filename +'\n')
        exit(-1)
    os.chdir(startingdir)    


# RS 2012/05/09:  Override chip-by-chip with flat astronet solution.
# This is really just a debugging tool, not meant for production use.

if args.astronet:
    print "Overriding Brian's WCS code with astrometry.net for each chip"
    for chip in range(1,33):
        os.chdir("%02d" % chip)
        file = os.path.splitext(filename)[0] + '_' + str(chip) + '.fits'
        wcsfile = file.replace('.fits','_wcs.fits')
        os.system("solve-mosaic_single.py {0} {1}".format(file, wcsfile))
        if not os.path.exists(wcsfile):
            sys.stderr.write('SM-WCS.py:  ASTRONET failed on image '+ file+'\n')
    exit(0)

# RS 2012/05/09:  Override chip-by-chip with flat astronet solution.

#
# We know where we are pointing, now start puting WCSs in
#
#close my first images header


for chip in range(1,33):
    # os.chdir(str(chip))
    os.chdir("%02d" % chip)
    file=filename + '_' + str(chip) + '.fits'
    file=os.path.splitext(filename)[0] + '_' + str(chip) + '.fits'
    
    # get header for this image
    hdulist = pyfits.open(file,mode='update')
    hdr = hdulist[0].header
    #print hdr

    #put rough WCS into header
    roughwcs(RACENT,DECCENT,hdr,tan,3.5)
    #write this header 
    hdulist.flush()
    hdulist.close()    


    #now run WCS-ZPN software

    if ucac2 == True:
        astrocat = "-ucac2 "
    else:
        astrocat = "-2mass "

    # first do a WCS match, but do not apply
    system=brianwcspath+'/bin/WCS.pl ' + file + ' -usewcs ' + astrocat + '-donotapply -force'
    os.system(system)

    if noresid == False:
        fitres="-fitresid "
    else:
        fitres="";
        

    
    if os.path.exists(os.path.splitext(file)[0] + '.fits.wcsmatch')==False:
        sys.stderr.write('WCS.pl failed on image '+ file +'\n')
    else: 
       #this fits a grid
       system=brianwcspath+'/bin/brian_fitwcs ' + file + ' ' + file + '.wcsmatch ' + fitres + ' -outputgrid 2048 4096 50 50' # -outfile'
       print system
       # need to add writing this info to each header, reading in what the program outputs
       os.system(system)

       if os.path.exists(os.path.splitext(file)[0] + '.fits.wcsmatch.grid')==False:
           sys.stderr.write('brian_fitwcs failed on image '+ file +'\n')
       else:
           # finally use the grid produced by above to put a 5 TNX into the header via IRAF Yuk!

           system=brianwcspath+'/bin/updateWCS.py ' + file + ' --matchlist ' + file \
               + '.wcsmatch.grid --order 5 --refpt \'' + RACENT + ' ' + DECCENT + '\''
           print system
           os.system(system)    

    #done here, go up to first directory
    
    os.chdir(startingdir)
