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
#FY - modified to just apply astrometry.net solution

import argparse
import numpy
import pyfits
import math
import os
import sys
import ephem
from STAP_comm import STAP_callexternal
from STAP_SEx import SEx
from Utils.Catalog import SExtractorDetection
from Utils.Photometry import calc_seeing, calc_elong


#
#subroutines
#

#
# write header of MEF with rough WCS based on RA-CENT,DEC-CENT or 
#

def roughwcs(raval,decval,header,tan,projp3,rotskypa_kwval):
    import math
    
    if header.get('IMAGEID')==None:
        print 'NO IMAGEID in image header'
        exit(-1)
    else:
        imageid=header.get('IMAGEID')

    # RS 2012/07/05:  Disabling this block.  We now get the rotator angle
    # directly from the astrometry.net initial guess.
    # if header.get('ROTSKYPA')==None:
    #     print 'NO ROTSKYPA in image header - setting to 0'
    #     rotskypa=0
    # else:
    #     rotskypa=header.get('ROTSKYPA')    
    rotskypa = rotskypa_kwval

    if header.get('NAXIS1')==None:
        print 'NO NAXIS1 in image header'
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
        print 'image is not valid size, expecting naxis1=2048'
        exit(-1)

    scale=0.497/3600 # SkyMapper scale - hardwired, could put in configuration

    rotskypa+=90
    rotskypa=rotskypa*math.pi/180.
    flipx=1
    flipy=-1
    coo = ephem.Equatorial(ephem.hours(raval), ephem.degrees(decval), epoch=2000)

    CRPIX1=crpix1[chip]
    CRPIX2=crpix2[chip]

    tpa=coo.ra*180/math.pi
    tpd=coo.dec*180/math.pi

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
    header.update('ROTSKYPA',rotskypa_kwval,"Rotator angle")

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
    
    # RS 2012/11/29:  what the hell, put AZ and EL in there too.
    # We'll need the MJD; try to extract it from the header
    if 'MJD-OBS' in header:
        # Siding Spring Observatory information; should be in Utils.Constants
        sso = ephem.Observer()
        sso.lon, sso.lat = '149.066086', '-31.277039'
        # ephem runs on this stupid "Dublin Julian Date"; convert from MJD
        sso.date = ephem.Date(float(header['MJD-OBS']) + 2400000.5 - 2415020)
        # Make a PyEphem "body" instance and compute it
        body = ephem.readdb("{id},f,{ra},{dec},0.0,2000.0".format(
                             id=0, ra=str(coo.ra), dec=str(coo.dec)))
        body.compute(sso)
        header.update('AZ', body.az*180/math.pi)
        header.update('EL', body.alt*180/math.pi)

    #NOTE these need to be committed when done using flush
    
    
#
# return crval for an image for a given amp using WCS on it.
#
def getcrval(amp,image):

# make system call to xy2sky
# xy2sky image.fits CRPIX1 CRPIX2

  # RS 2012/06/03:  Now assumes WCStools are in the user's $PATH.
  system=  'xy2sky '+ image + ' ' + str(crpix1[amp]) + ' ' + str(crpix2[amp]) 
  print system
  sky2xy = os.popen(system)
  line=sky2xy.readline()
  racent=line.split()[0]
  deccent=line.split()[1]
  return(racent,deccent)

############################################################
# MAIN PROGRAM
#

homepath=os.environ["SUBPIPEHOME"]

configurefile= homepath + '/WCS/etc/configure.wcs.skymapper'

parser = argparse.ArgumentParser(description='Take SkyMapper FITS files and put rough WCS on it')
parser.add_argument('imname', help='input filename (files assumed to be in (path/)N/filename_N.fits)')
parser.add_argument('-impath', help='input file path',default='./')
                    
args=parser.parse_args()
filename=args.imname
impath=args.impath

#
# readin CRPIX from configure.wcs.skymapper
#
crpix1,crpix2=numpy.loadtxt(configurefile, usecols=(1,2), unpack=True)

  
#
# Do we want to run astronet to get the position of the chip
# or the center of the array if a MEF
# If so, make a temporary copy of file
#

useext=impath+'/17/'+os.path.splitext(filename)[0] + '_17.fits'
solvedfits=os.path.splitext(useext)[0] + '_solved.fits'

hdulist = pyfits.open(useext)
prihdr = hdulist[0].header
if prihdr.get('NEXTEND')!= None:
    print  'SM-ASTROMETRY.py:  Not a single FITS file, bye'
    exit(-1)

# RS 2012/03/18:  don't bother trying to align flat field images
if 'IMAGETYP' in prihdr and prihdr['IMAGETYP'] != 'object':
    print "SM-ASTROMETRY.py:  Sorry, image type is '{0}' "\
          "and not a science image".format(prihdr['IMAGETYP'])
    print "Exiting normally, but not trying to solve for WCS either."
    exit(0)

# RS 2012/03/18:  even so, check to see if there are any objects.
# Yucky hardcoded stuff here, beware!
try:
    SEx(useext, do_seeing=True)
except:
    print "SM-ASTROMETRY.py:  SExtractor failed on image", useext
    exit(-2)
tmpcat = os.path.basename(useext + ".stars")
detections = SExtractorDetection.read_fits_file(tmpcat)[0]
n_objs = len(detections)
if n_objs < 20:
    print "SM-ASTROMETRY.py:  Only", n_objs, "objects found in image", useext
    print "Not enough to solve field; bailing."
    exit(-3)
# RS 2012/06/28:  Calculate a rough seeing estimate from chip 17.
qq_seeing = calc_seeing(detections, verbose=True, quartiles=True)
qq_elong = calc_elong(detections, verbose=True, quartiles=True)
seeing, seewid = qq_seeing[1], (qq_seeing[2]-qq_seeing[0])/qq_seeing[1]
elong, elongwid = qq_elong[1], (qq_elong[2]-qq_elong[0])/qq_elong[1]


# RS 2012/03/18:  getting on with it
ra=prihdr['RA']
dec=prihdr['DEC']
hdulist.close()


cmd='solve-field -T --parity neg --no-plots --no-fits2fits '+  useext + ' --scale-low .493 --scale-high 0.501 --scale-units app --ra ' + ra + ' --dec ' + dec + ' --radius 15 --overwrite --use-sextractor --sextractor-path /usr/local/bin/sex -N ' + solvedfits
[status,stdoutstr]=STAP_callexternal(cmd,timeout=600,getstdout=True,combinestderr=True)
print stdoutstr
if status !=0 :
    print "SM-ASTROMETRY.py:  astrometry.net terminated abnormally on image %s"%filename
    exit(-4)
else:
    if os.path.exists(os.path.splitext(useext)[0] + '.solved')==True:
        (RACENT,DECCENT)=getcrval(33,solvedfits)
    else:
        print "SM-ASTROMETRY.py:  astrometry.net failed to solve image %s"%filename
        exit(-5)

#
# We know where we are pointing, now start puting WCSs in
#
#close my first images header

# RS 2012/07/05:  Extract ROTSKYPA information from the header of chip 17.
# This should generally be more reliable than the header information, and if
# astrometry.net died we didn't want to work with that image anyways.
hdr = pyfits.getheader(solvedfits)
rotskypa = math.atan2(hdr['CD1_2'], hdr['CD1_1']) * 180/math.pi
if rotskypa < 0:  rotskypa += 360.0

for chip in range(1,33):
    chippath='%s/%02d/'%(impath,chip)
    file=chippath+os.path.splitext(filename)[0] + '_' + str(chip) + '.fits'
    
    # get header for this image
    hdulist = pyfits.open(file,mode='update')
    hdr = hdulist[0].header
    #print hdr

    #put rough WCS into header
    roughwcs(RACENT,DECCENT,hdr,False,3.5,rotskypa)
    # RS 2012/06/28:  Put rough seeing estimate into header.
    hdr.update("SEEING", seeing, "SExtractor-calculated median seeing in pixels")
    hdr.update("SEEWID", seewid, "75%-25% fractional seeing distribution width")
    hdr.update("ELONG", elong, "SExtractor-calculated median elongation")
    hdr.update("ELONGWID", elongwid, "75%-25% fractional elongation distrib. width")
    #write this header 
    hdulist.flush()
    hdulist.close()    
