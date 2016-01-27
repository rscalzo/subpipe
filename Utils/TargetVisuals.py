#!/usr/bin/env python

"""
RS 2014/02/01:  Target Visuals Module

Useful for corralling subs, making finding charts etc.
"""

import os
import sys
import pyfits
import subprocess
import numpy as np
from Utils.Constants import FilenamesSub, FilenamesXsient
from Utils.Constants import PipelinePath as CPP
import mydjango.followup.models as fu
import mydjango.jobstats.models as js
from STAP.STAP_display import zscale, rescale, sky2xy
import matplotlib.pyplot as pypl
from matplotlib import cm


def get_sublist(xsient_name, filter='g', only_detections=True):
    """Look up all subtractions in the database involving a given transient"""

    # Look up transient in database by name
    try:
        xsient = fu.Transient.objects.get(name=xsient_name)
    except fu.Transient.DoesNotExist:
        print "run_ds9:  Can't find {0} in database".format(xsient_name)
        return
    # pointinglist = xsient.pointing_det.get_query_set().filter(filter=filter)
    if only_detections:
        pointingmgr = xsient.pointing_det
    else:
        pointingmgr = xsient.pointing_obs
    pointinglist = pointingmgr.get_query_set().filter(filter=filter)

    # Get a list of all the subtraction jobs in this field and whether they
    # actually overlap with the transient's location.
    sublist = js.SubtractionJob.objects.select_related(depth=1).filter(
            ccd=xsient.ccd, pointing__in=pointinglist,
            minra__lt=xsient.ra, mindec__lt=xsient.dec,
            maxra__gt=xsient.ra, maxdec__gt=xsient.dec)
    print "{0} was observed in {1} pointings".format(xsient_name, len(sublist))
    return xsient, sublist


def cutout_bounds(x, y, nx, ny, naxis1, naxis2, verbose=False):
    """
    Grabs the coordinate limits of a cutout

    Inputs
        x, y:  Pixel coordinates of cutout center within image
        nx, ny:  Extent of cutout in pixels
        naxis1, naxis2:  Size of image in pixels
    """
    # Default extent of cutout
    xc0, xc1, yc0, yc1 = 0, nx, 0, ny
    # Boundaries of cutout within image if cutout is wholly contained 
    xm0, ym0 = int(np.round(x)) - nx/2, int(np.round(y)) - ny/2
    xm1, ym1 = xm0 + nx, ym0 + ny
    # Deal with edge cases involving partial overlap
    if xm0 < 0:
        xc0, xm0 = -xm0, 0
    if xm1 > naxis1:
        xc1, xm1 = xc1 - (xm1 - naxis1), naxis1
    if ym0 < 0:
        yc0, ym0 = -ym0, 0
    if ym1 > naxis2:
        yc1, ym1 = yc1 - (ym1 - naxis2), naxis2
    # Return final array index bounds in cutout and in larger image
    returnvals = ((xc0, xc1, yc0, yc1), (xm0, xm1, ym0, ym1))
    if verbose:
        print "Utils.TargetVisuals.cutout_bounds():", returnvals
    return returnvals


def makefinder(xsient_name):
    """Make a finding chart using matplotlib"""

    # Get a list of all the subs in which stuff was *detected*
    xsient, sublist = get_sublist(xsient_name, only_detections=True)

    # Find discovery image
    newfn = None
    for sub in sublist:
        submeta = FilenamesSub(subid=sub.jobname)
        subfitspath = submeta.absolute_dir + '/'
        subfitsname = subfitspath + submeta.sub_fits_name + '.gz'
        if os.path.exists(subfitsname):
            subhdr = pyfits.getheader(subfitsname)
            trial_newfn = subfitspath + subhdr['TARGET'] + ".gz"
            if os.path.exists(trial_newfn):
                newfn = trial_newfn
                break
    if not newfn:
        print "Error:  couldn't find discovery image for", xsient.name
        return None

    # Plot the image
    # This is all adapted from STAP.STAP_display.STAP_make_thumbs_xy
    with pyfits.open(newfn) as hdulist:
        hdr, data = hdulist[0].header, hdulist[0].data
    scratchdir = CPP.scratch + "/subinspect/"
    subprocess.call(['mkdir', '-p', scratchdir])
    # Find x and y coordinates of transient in image
    #hdrfname = os.path.basename(newfn).replace(".fits.gz", ".hdr")
    #hdrfname = scratchdir + hdrfname
    #with open(hdrfname, 'w') as hdrfile:
    #    hdrfile.write(str(hdr))
    fname = scratchdir + os.path.basename(newfn)
    fitsname= fname.replace(".fits.gz", ".fits")
    os.system('cp %s %s'%(newfn,scratchdir))
    os.system('gunzip %s'%(fname))
    idx, idy = sky2xy(fitsname, [xsient.ra], [xsient.dec])
    idx, idy = int(idx[0]), int(idy[0])
    print "{0} appears at ({1}, {2}) in original coords".format(
        xsient.name, idx, idy)
    os.system('rm -fr %s'%fitsname)

    # Reorient according to WCS
    # This incantation should get the data facing with N = up, E = left.
    data, idx, idy = np.rot90(data.T), idx, data.shape[0] - idy
    if hdr['CD2_2'] > 0:
        print "Rotating 180 degrees"
        idx, idy = data.shape[1] - idx, data.shape[0] - idy
        data = np.flipud(np.fliplr(data))
    print "new x, y =", idx, idy
    print "data.shape =", data.shape
    # Remove subimage roughly corresponding to finding chart,
    # respecting image boundaries
    lenx, leny = 1024, 1024
    xlo, xhi = idx - lenx/2, idx + lenx/2
    ylo, yhi = idy - leny/2, idy + leny/2
    if xlo < 0:  xlo, xhi = 0, lenx
    if ylo < 0:  ylo, yhi = 0, leny
    if xhi > data.shape[1]:
        xlo, xhi = data.shape[1] - lenx - 1, data.shape[1] - 1
    if yhi > data.shape[0]:
        ylo, yhi = data.shape[0] - leny - 1, data.shape[0] - 1
    print "xlo, xhi, ylo, yhi = [{0}:{1},{2}:{3}]".format(
            xlo, xhi, ylo, yhi)
    data = data[ylo:yhi, xlo:xhi]
    idx, idy = idx - xlo, idy - ylo
    # Rescale using IRAF zscale algorithm from stsci.numdisplay
    slo, shi = zscale(data)
    data = rescale(data, slo, shi)
    print "data.shape =", data.shape
    print "new idx, idy =", idx, idy
    sys.stdout.flush()

    # Display
    pypl.figure() # (figsize=(8,8))
    ax=pypl.axes([0,0,1,1],frameon=False)
    ax.set_axis_off()
    # Rotate 270 degrees...
    # for i in range(3): data = np.rot90(data)
    pypl.imshow(data, cmap=cm.gray, origin='lower', interpolation='nearest')
    # Crosshairs
    # pypl.plot(idx, idy, marker='s', ms=12, mew=3, mec='lightgreen', mfc='None')
    for i in range(4):
        vec = [np.cos(i*np.pi/2), np.sin(i*np.pi/2)]
        x = idx + np.array([30*vec[0], 70*vec[0]])
        y = idy + np.array([30*vec[1], 70*vec[1]])
        pypl.plot(x, y, lw=1.5, color='lightgreen')
    # Object name
    pypl.text(70, leny-100, xsient.name,
              ha='left', va='center', color='lightgreen')
    # Compass rose: first E arrow, then N arrow
    pypl.arrow(200, 100, -100, 0, width=1.5, head_width=15, color='lightgreen')
    pypl.text( 70, 100, 'E', ha='right', va='center', color='lightgreen')
    pypl.arrow(200, 100, 0, +100, width=1.5, head_width=15, color='lightgreen')
    pypl.text(200, 230, 'N', ha='center', va='bottom', color='lightgreen')
    # Angular scale
    pypl.plot([lenx-220, lenx-100], [100, 100], lw=1.5, color='lightgreen')
    pypl.text(lenx-160, 120, '1 arcmin', ha='center', color='lightgreen')

    # Finally, show/save
    subprocess.call(['mkdir', '-p', CPP.finderpath])
    xmeta = FilenamesXsient(xsient.name, xsient.field.id, xsient.ccd)
    fcfname = xmeta.finderfname
    print "Saving finding chart in", fcfname
    pypl.savefig(fcfname)
    pypl.show()

def ds9history(xsient_name):
    """Run ds9 on all existing subtractions"""

    # Get a list of all the subs in which the area was *observed*
    xsient, sublist = get_sublist(xsient_name, only_detections=False)

    # Make a ds9 regions file for the SN
    scratchdir = CPP.scratch + "/subinspect/"
    subprocess.call(['mkdir', '-p', scratchdir])
    subfnregSN = scratchdir + "{0}.reg".format(xsient.name)
    with open(subfnregSN, 'w') as snregfile:
        snregfile.write("global color=blue width=2\nj2000\n")
        snregfile.write("j2000; box({0},{1},0.004,0.004,0) "
                        "# color=cyan,width=2".format(xsient.ra, xsient.dec))

    # Now construct a list of files to send to ds9
    ds9pars, reffn = [ ], None
    for sub in sublist:
        submeta = FilenamesSub(subid=sub.jobname)
        subfitspath = submeta.absolute_dir + '/'
        subfitsname = subfitspath + submeta.sub_fits_name + '.gz'
        subfnreg = subfitspath + submeta.sub_reg_name
        if os.path.exists(subfitsname) and os.path.exists(subfnreg):
            ds9pars.append(subfitsname)
            ds9pars.extend(['-region', subfnreg, '-region', subfnregSN])
            # Find a NEW or REF as well.  The filename will be listed in
            # the FITS header as 'TARGET' for NEW, 'TEMPLATE' for REF.
            if not reffn:
                subhdr = pyfits.getheader(subfitsname)
                trial_reffn = subfitspath + subhdr['TARGET'] + ".gz"
                if os.path.exists(trial_reffn):
                    reffn = trial_reffn
        else:
            print "Couldn't find", subfitsname, "or its region file"
    if reffn:
        ds9pars.extend([reffn, '-region', subfnregSN])

    # Now run ds9!
    print "Launching ds9"
    pars = ['ds9', '-rotate', '90', '-zoom', 'to', '0.25', '-zscale'] + ds9pars
    print " ".join(pars)
    sys.stdout.flush()
    subprocess.call(pars)
    return
