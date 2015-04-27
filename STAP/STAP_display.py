#!/usr/bin/env python

#make thumbnails centered at list of ra/dec's


#imname='../data/new/111-18/g/15/Skymapper_805310385_00000_2011-04-12T19:16:27_15.fits'

#ras=[168.98898]
#decs=[-18.014603]
#outdir='.'

import pyfits
import numpy as np
from math import pi
import subprocess
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
from matplotlib import cm
import os
import sys
import ephem

def STAP_make_thumbs_xy(imname,xs,ys,box=120,outdir=None,usesubid=None,objnames=None,suffix='',outputnames=None):
    if len(xs) != len(ys):
        print "input xs/ys arrays don't match"
        return(1)
    if outputnames:
        if len(xs) != len(outputnames):
            print "length of output stamp names do not match length of coordinates arrays"
            return(1)
    if objnames:
        if len(xs) != len(objnames):
            print "length of object names do not match length of coordinates arrays"
            return(1)
    if not os.path.exists(imname):
        imname=imname+'.gz'
        if not os.path.exists(imname):
            print "can't find image %s"%imname
            return(1)
    try:
        fullimdata=pyfits.getdata(imname)
        imheader=pyfits.getheader(imname)
    except:
        print "can't read image %s"%imname
        return(1)

    if not outputnames:
        if usesubid is None:
            imnamebase=os.path.splitext(os.path.basename(imname))[0]
            imnamebase='_'.join(imnamebase.split('_')[0:4])
        else:
            imnamebase=usesubid
        
    if outdir==None:
        outdir=os.path.dirname(imname)
        if outdir in '':
            outdir='.'
            
    [imy,imx]=fullimdata.shape
    slow,shigh=zscale(fullimdata)
    nstamps=len(xs)

    for ind in range(nstamps):
        xct=xs[ind]-1.
        yct=ys[ind]-1.

        xll=xct-box/2
        xhh=xct+box/2
        yll=yct-box/2
        yhh=yct+box/2
           
        #not really necessary, always xll <xhh, yll<yhh
        xlow=np.floor(np.min([xll,xhh]))
        xhigh=np.ceil(np.max([xll,xhh]))
        ylow=np.floor(np.min([yll,yhh]))
        yhigh=np.ceil(np.max([yll,yhh]))
        
        if xlow <0:xlow=0
        if ylow <0:ylow=0
        if xhigh >imx-1: xhigh=imx-1
        if yhigh >imy-1: yhigh=imy-1
        
        if yhigh<=ylow or xhigh<=xlow:
            if objnames:
                print "Stamp size less than 1 for %s"%objnames[ind]
            else:
                print "Stamp size less than 1 for %f/%f"%(xs[ind],ys[ind])
            continue
        
        imdata=fullimdata[ylow:yhigh+1,xlow:xhigh+1]
        xrel=xct-xlow
        yrel=yct-ylow
        
        #figure out rotation angle
        #from header
        if ('CD1_1' in imheader) and ('CD2_2' in imheader):
            cd1_1=imheader['CD1_1']
            cd2_2=imheader['CD2_2']
            if cd1_1<0 and cd2_2>0: rotated=imdata
            if cd1_1<0 and cd2_2<0: 
                rotated=np.fliplr(np.rot90(imdata,k=2))
                yrel=yhigh-yct
            if cd1_1>0 and cd2_2>0: 
                rotated=np.fliplr(imdata)
                xrel=xhigh-xct
            if cd1_1>0 and cd2_2<0: 
                rotated=np.rot90(imdata,k=2)
                xrel=xhigh-xct
                yrel=yhigh-yct
        else: 
            print "CD1_1 and/or CD2_2 are not found in header, thumb orientation might be incorrect"
            #flip left/right by default
            rotated=np.fliplr(imdata)
       
        #scale the data
        #slow,shigh=zscale(rotated)
        rotated=rescale(rotated,slow,shigh)
        
        #output file name
        if outputnames:
            stampname=outdir + '/%s'%outputnames[ind]
        else:
            if objnames:    
                stampname=outdir+'/%s_%s%s.png'%(objnames[ind],imnamebase,suffix)
            else:
                stampname=outdir+'/'+imnamebase+'_%04.0f_%04.0f%s.png'%(xs[ind],ys[ind],suffix)
        
        #plot data
        pyplot.figure(figsize=(4,4))
        ax=pyplot.axes([0,0,1,1],frameon=False)
        ax.set_axis_off()
        pyplot.imshow(rotated,cmap=cm.gray,origin='lower',interpolation='nearest')

        #add cross hair
        xsize=xhigh-xlow
        ysize=yhigh-ylow
        ssize=np.min([xsize,ysize])
        hline=np.array([xrel-ssize/3.,xrel-ssize/10.,xrel+ssize/10.,xrel+ssize/3.])/xsize
        vline=np.array([yrel-ssize/3.,yrel-ssize/10.,yrel+ssize/10.,yrel+ssize/3.])/ysize
        pyplot.axhline(y=yrel,xmin=hline[0],xmax=hline[1],c='1')
        pyplot.axhline(y=yrel,xmin=hline[2],xmax=hline[3],c='1')
        pyplot.axvline(x=xrel,ymin=vline[0],ymax=vline[1],c='1')
        pyplot.axvline(x=xrel,ymin=vline[2],ymax=vline[3],c='1')
        
        pyplot.savefig(stampname)

def STAP_make_thumbs(imname,ras,decs,box=1.,pbox=None,outdir=None,usesubid=None,objnames=None,suffix=''):
    if len(ras) != len(decs):
        print "input ra/dec arrays don't match"
        return(1)
    if not objnames is None:
        if len(ras) != len(objnames):
            print "length of object names do not match length of coordinates arrays"
            return(1)
    if not os.path.exists(imname):
        imname=imname+'.gz'
        if not os.path.exists(imname):
            print "can't find image %s"%imname
            return(1)

    nstamps=len(ras)
    [xs,ys]=sky2xy(imname,ras,decs)
    
    if not pbox:
        #figure out box size in pixel
        header=pyfits.gethead(imname)
        if "CD1_1" in header: 
            scale=header["CD1_1"]
            pbox=round(box/60./scale/2.)*2.
        else:
            ras=np.array(ras)
            decs=np.array(decs)
            hboxra=box/2./60./np.cos(decs*pi/180.)
            hboxdec=box/2./60.
            ralimls=ras-hboxra
            ralimhs=ras+hboxra
            declimls=decs-hboxdec
            declimhs=decs+hboxdec
            [xlls,ylls]=sky2xy(imname,ralimls,declimls)
            [xhhs,yhhs]=sky2xy(imname,ralimhs,declimhs)
            pbox=np.mean(xhhs-xlls)
            pbox2=np.mean(yhhs-ylls)
            if pbox2 >pbox: pbox=pbox2
            pbox=round(pbox/2.)*2.

    if not objnames:
        objnames=[]
        for ind in range(nstamps):
            radecstr=radec2str(ras[ind],decs[ind])
            objnames.append('J'+radecstr)
    
    STAP_make_thumbs_xy(imname,xs,ys,box=pbox,outdir=outdir,usesubid=usesubid,objnames=objnames,suffix=suffix)    

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
    
def sky2xy(imname,ras,decs):
    cmd="sky2xy %s"%(imname)
    nobj=len(ras)
    for ind in range(nobj):
        cmd="%s %f %f"%(cmd,ras[ind],decs[ind])
    cmd_orig = cmd
    cmd = cmd.split(" ")
    process=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    (stdoutstr,stderrstr) = process.communicate()
    status = process.returncode
    if status != 0:
        print cmd_orig, "failed with the following message:"
        print stderrstr
        print "Couldn't convert ra/dec to image coordinates"
        return [[-999]*nobj,[-999]*nobj]
    else:
        lines=stdoutstr.strip().split('\n')
        if len(lines) !=nobj:
            print "Output does not match input number of objects"
            return [[-999]*nobj,[-999]*nobj]
        else:
            x=[]
            y=[]
            for ind in range(nobj): 
                x.append(float(lines[ind].split()[4]))
                y.append(float(lines[ind].split()[5]))
            return [x,y]

def xy2sky(imname,xs,ys):
    cmd="xy2sky %s"%(imname)
    nobj=len(xs)
    for ind in range(nobj):
        cmd="%s %f %f"%(cmd,xs[ind]+1.,ys[ind]+1.)
    cmd = cmd.split(" ")
    process=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    (stdoutstr,stderrstr) = process.communicate()
    status = process.returncode
    if status != 0:
        print "Couldn't convert x/y to sky coordinates"
        return [[-999]*nobj,[-999]*nobj]
    else:
        lines=stdoutstr.strip().split('\n')
        if len(lines) !=nobj:
            print "Output does not match input number of objects"
            return [[-999]*nobj,[-999]*nobj]
        else:
            ras=[]
            decs=[]
            for ind in range(nobj):
                ras.append(ephem.hours(stdoutstr.split()[0])*180./pi)
                decs.append(ephem.degrees(stdoutstr.split()[1])*180./pi)
            return [ras,decs]

def rescale(im,z1,z2):
    z21=z2-z1
    if z21==0: z21=1. 
    imdata=(im-z1)/z21
    imdata=np.clip(imdata,0,1)
    return imdata

def zscale (image, nsamples=1000, contrast=0.25, bpmask=None, zmask=None): 
    """Implement IRAF zscale algorithm 
    nsamples=1000 and contrast=0.25 are the IRAF display task defaults 
    bpmask and zmask not implemented yet 
    image is a 2-d numpy array 
    returns (z1, z2) 
    """ 
    MAX_REJECT = 0.5 
    MIN_NPIXELS = 5 
    KREJ = 2.5 
    MAX_ITERATIONS = 5 
    
    # Sample the image 
    samples = zsc_sample (image, nsamples, bpmask, zmask) 
    npix = len(samples) 
    samples.sort() 
    zmin = samples[0] 
    zmax = samples[-1] 
    # For a zero-indexed array 
    center_pixel = (npix - 1) // 2 
    if npix%2 == 1: 
        median = samples[center_pixel] 
    else: 
        median = 0.5 * (samples[center_pixel] + samples[center_pixel + 1]) 
        
    # 
    # Fit a line to the sorted array of samples 
    minpix = max(MIN_NPIXELS, int(npix * MAX_REJECT)) 
    ngrow = max (1, int (npix * 0.01)) 
    ngoodpix, zstart, zslope = zsc_fit_line (samples, npix, KREJ, ngrow, 
                                             MAX_ITERATIONS) 
  
    if ngoodpix < minpix: 
        z1 = zmin 
        z2 = zmax 
    else: 
        if contrast > 0: zslope = zslope / contrast 
        z1 = max (zmin, median - (center_pixel - 1) * zslope) 
        z2 = min (zmax, median + (npix - center_pixel) * zslope) 
    return z1, z2
    
def zsc_sample (image, maxpix, bpmask=None, zmask=None): 
     
    # Figure out which pixels to use for the zscale algorithm 
    # Returns the 1-d array samples 
    # Don't worry about the bad pixel mask or zmask for the moment 
    # Sample in a square grid, and return the first maxpix in the sample 
    nc = image.shape[0] 
    nl = image.shape[1] 
    stride = max (1.0, np.sqrt((nc - 1) * (nl - 1) / float(maxpix))) 
    stride = int (stride) 
    samples = image[::stride,::stride].flatten() 
    return samples[:maxpix] 
     
def zsc_fit_line (samples, npix, krej, ngrow, maxiter): 
 
    MAX_REJECT = 0.5 
    MIN_NPIXELS = 5 
    GOOD_PIXEL = 0 
    BAD_PIXEL = 1 
    # 
    # First re-map indices from -1.0 to 1.0 
    xscale = 2.0 / (npix - 1) 
    xnorm = np.arange(npix) 
    xnorm = xnorm * xscale - 1.0 
 
    ngoodpix = npix 
    minpix = max (MIN_NPIXELS, int (npix*MAX_REJECT)) 
    last_ngoodpix = npix + 1 
 
    # This is the mask used in k-sigma clipping.  0 is good, 1 is bad 
    badpix = np.zeros(npix, dtype="int32") 
 
    # 
    #  Iterate 
 
    for niter in range(maxiter): 
 
        if (ngoodpix >= last_ngoodpix) or (ngoodpix < minpix): 
            break 
         
        # Accumulate sums to calculate straight line fit 
        goodpixels = np.where(badpix == GOOD_PIXEL) 
        sumx = xnorm[goodpixels].sum() 
        sumxx = (xnorm[goodpixels]*xnorm[goodpixels]).sum() 
        sumxy = (xnorm[goodpixels]*samples[goodpixels]).sum() 
        sumy = samples[goodpixels].sum() 
        sum = len(goodpixels[0]) 
 
        delta = sum * sumxx - sumx * sumx 
        # Slope and intercept 
        intercept = (sumxx * sumy - sumx * sumxy) / delta 
        slope = (sum * sumxy - sumx * sumy) / delta 
         
        # Subtract fitted line from the data array 
        fitted = xnorm*slope + intercept 
        flat = samples - fitted 
 
        # Compute the k-sigma rejection threshold 
        ngoodpix, mean, sigma = zsc_compute_sigma (flat, badpix, npix) 
 
        threshold = sigma * krej 
 
        # Detect and reject pixels further than k*sigma from the fitted line 
        lcut = -threshold 
        hcut = threshold 
        below = np.where(flat < lcut) 
        above = np.where(flat > hcut) 
 
        badpix[below] = BAD_PIXEL 
        badpix[above] = BAD_PIXEL 
         
        # Convolve with a kernel of length ngrow 
        kernel = np.ones(ngrow,dtype="int32") 
        badpix = np.convolve(badpix, kernel, mode='same') 
 
        ngoodpix = len(np.where(badpix == GOOD_PIXEL)[0]) 
         
        niter += 1 
 
    # Transform the line coefficients back to the X range [0:npix-1] 
    zstart = intercept - slope 
    zslope = slope * xscale 
 
    return ngoodpix, zstart, zslope 
 
def zsc_compute_sigma (flat, badpix, npix): 

    GOOD_PIXEL = 0 
 
    # Compute the rms deviation from the mean of a flattened array. 
    # Ignore rejected pixels 
 
    # Accumulate sum and sum of squares 
    goodpixels = np.where(badpix == GOOD_PIXEL) 
    sumz = flat[goodpixels].sum() 
    sumsq = (flat[goodpixels]*flat[goodpixels]).sum() 
    ngoodpix = len(goodpixels[0]) 
    if ngoodpix == 0: 
        mean = None 
        sigma = None 
    elif ngoodpix == 1: 
        mean = sumz 
        sigma = None 
    else: 
        mean = sumz / ngoodpix 
        temp = sumsq / (ngoodpix - 1) - sumz*sumz / (ngoodpix * (ngoodpix - 1)) 
        if temp < 0: 
            sigma = 0.0 
        else: 
            sigma = np.sqrt (temp) 
 
    return ngoodpix, mean, sigma 

