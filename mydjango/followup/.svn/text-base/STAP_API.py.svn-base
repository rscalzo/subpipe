#!/usr/bin/env python

# ============================================================================
# RS 2012/06/14:  Django API for pipeline
# ----------------------------------------------------------------------------
# This allows the pipeline to interact with followup.models.
# ============================================================================

import ephem
import datetime
import numpy as np
from Utils import Constants
import mydjango.jobstats.models as js
import mydjango.followup.models as fu
from django.contrib.auth.models import User
from django.db.transaction import commit_manually, commit


def register_xsient(xdict):
    """Registers a new transient in Django.  Called by subpipe_master.py.

    xdict:  A list of keywords for constructing the Transient record.
    """

    try:
        xsient = fu.Transient.objects.get(name=xdict['name'],field__id=xdict['field'],ccd=xdict['subfield'])
        # RS 2013/09/25:  Fixed bug re: overwriting autotypes.  Don't update
        # the type unless it is explicitly overridden in the hash we passed in.
        if 'autotype' in xdict:
            xsient.type = fu.TransientType.objects.get_or_create(
                    type=xdict['autotype'])[0]
    except fu.Transient.DoesNotExist:
        # Fill some fields with boilerplate if we had to create a new entry.
        if 'autotype' not in xdict:  xdict['autotype'] = '?'
        newtype = fu.TransientType.objects.get_or_create(type=xdict['autotype'])[0]
        status = fu.FollowUpStatus.objects.get_or_create(status='new')[0]
        field = js.SkymapperField.objects.get(id=xdict['field'])
        xsient = fu.Transient(name=xdict['name'],
                              ra=xdict['ra'], dec=xdict['dec'],
                              ccd=xdict['subfield'], field=field,
                              type=newtype,
                              follow_up_status=status)

    # Save comments to database
    if (not xsient.id) or 'comment' in xdict: 
        if 'comment' in xdict:
            comment=xdict['comment']
        else:
            comment=''
        user = User.objects.get_or_create(username='skybot')[0]
        xsient.savelog(user=user,comment=comment)

    for ix,jd in enumerate(xdict['jd']):
        dateobs=ephem.Date(jd-2415020)
        daterange=[dateobs.datetime()-datetime.timedelta(seconds=2),dateobs.datetime()+datetime.timedelta(seconds=2)]
        
        smpointing = js.SciencePointing.objects.filter(dateobs__range=daterange)
        if len(smpointing)==0:
            print "WARNING: Cannot find pointing for field: {0} at dateobs: {1}".format(xsient.field.id,dateobs)
            continue
        if len(smpointing)>1:
            print "WARNING: Found more than one pointing for field: {0} at dateobs: {1}".format(xsient.field.id,dateobs)
            continue
        xsient.pointing_obs.add(smpointing[0])
        if not xdict['islim'][ix]:
            xsient.pointing_det.add(smpointing[0])
    xsient.save()
        
    # Update transient history based on comments in hash
    #xsientfnames = Constants.FilenamesXsient(
    #        xdict['name'], xdict['field'], xdict['subfield'])
    #if 'comment' in xdict:
    #    with open(xsientfnames.abshistfname, "a") as commentfile:
    #        commentfile.write(xdict['comment'])


def register_xdict(xdict,smpointings,smdates):
    """Register a transient, with required input of smpointings and
    corresponding dates.
    To be called by register_xlist.
    """
    try:
        xsient = fu.Transient.objects.get(name=xdict['name'])
        if 'autotype' in xdict:
            xsient.type = fu.TransientType.objects.get_or_create(
                type=xdict['autotype'])[0]
    except fu.Transient.DoesNotExist:
        # Fill some fields with boilerplate if we had to create a new entry.
        if 'autotype' not in xdict:  xdict['autotype'] = '?'
        newtype = fu.TransientType.objects.get_or_create(type=xdict['autotype'])[0]
        status = fu.FollowUpStatus.objects.get_or_create(status='new')[0]
        field = js.SkymapperField.objects.get(id=xdict['field'])
        xsient = fu.Transient(name=xdict['name'],
                              ra=xdict['ra'], dec=xdict['dec'],
                              ccd=xdict['subfield'], field=field,
                              type=newtype,
                              follow_up_status=status)

    # Save comments to database
    if (not xsient.id) or 'comment' in xdict: 
        if 'comment' in xdict:
            comment=xdict['comment']
        else:
            comment=''
        user = User.objects.get_or_create(username='skybot')[0]
        xsient.savelog(user=user,comment=comment)

    for ix,jd in enumerate(xdict['jd']):
        dateobs=ephem.Date(jd-2415020).datetime()
        tdiff=np.abs(smdates-dateobs)
        tmin=tdiff.argsort()[0]
        if tdiff[tmin].seconds<2.:
            xsient.pointing_obs.add(smpointings[tmin])
            if not xdict['islim'][ix]:
                xsient.pointing_det.add(smpointings[tmin])
        else:
            print "WARNING: Cannot find pointing for field: {0} at dateobs: {1}".format(xsient.field.id,dateobs)
            continue
    xsient.save()

@commit_manually
def batch_create(xnew):
    for xdict in xnew:
        if 'autotype' not in xdict:  xdict['autotype'] = '?'
        newtype = fu.TransientType.objects.get_or_create(type=xdict['autotype'])[0]
        status = fu.FollowUpStatus.objects.get_or_create(status='new')[0]
        field = js.SkymapperField.objects.get(id=xdict['field'])
        xsient = fu.Transient.objects.create(name=xdict['name'],
                                             ra=xdict['ra'], dec=xdict['dec'],
                                             ccd=xdict['subfield'], field=field,
                                             type=newtype,
                                             follow_up_status=status)
    commt()

@commit_manually
def batch_update(xlist,xsients,smpointings,smdates):
    for xsient in xsients:
        xupdate=[xdict for xdict in xlist if xdict['name']==xsient.name]
        for xdict in xupdate:
            if 'autotype' in xdict:
                xsient.type = fu.TransientType.objects.get_or_create(
                    type=xdict['autotype'])[0]

            # Save comments to database
            if (not xsient.id) or 'comment' in xdict: 
                if 'comment' in xdict:
                    comment=xdict['comment']
                else:
                    comment=''
                user = User.objects.get_or_create(username='skybot')[0]
                xsient.savelog(user=user,comment=comment)

            for ix,jd in enumerate(xdict['jd']):
                dateobs=ephem.Date(jd-2415020).datetime()
                tdiff=np.abs(smdates-dateobs)
                tmin=tdiff.argsort()[0]
                if tdiff[tmin].seconds<2.:
                    xsient.pointing_obs.add(smpointings[tmin])
                    if not xdict['islim'][ix]:
                        xsient.pointing_det.add(smpointings[tmin])
                else:
                    print "WARNING: Cannot find pointing for field: {0} at dateobs: {1}".format(xsient.field.id,dateobs)
                    continue
            xsient.save()
    commit()


def register_xlist(xlist):
    """Register a list of new transients. Called by subpipe_master.py.
    This function is used instead of register_xsient to limit database queries
    and save time.

    xlist: result of xref, a dictionary of id, xdict
    """
    #the list will be broken into sub groups for each skymapper field
    fnames=[xdict['field'] for id, xdict in xlist.items()]
    fnames=list(set(fnames))
    for fname in fnames:
        smpointings=list(js.SciencePointing.objects.filter(field__id=fname))
        smdates=np.array([sm.dateobs for sm in smpointings])
        for id, xdict in xlist.items():
            if xdict['field']!=fname:
                continue
            register_xdict(xdict,smpointings,smdates)


