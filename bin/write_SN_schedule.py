#!/usr/bin/env python

"""
RS 2012/08/17:  SkyMapper SN Search Schedule Writing Script

This is a more advanced version of the hacky little script I ran during
the search commissioning.  Useful features:

    -- Now contains an objective function which it tries to optimize.
       It doesn't do a full nonlinear optimization a la plane schedules
       but at least there is now a way to incorporate observing
       constraints using a universal prescription.

    -- Incorporates standard stars, search fields and follow-up fields.

    -- Can provide single-line answers to interface with the scheduler.

FY 2015/08: 
    -- Use SkymapperField database table to find the observing constraints

    -- Get lastobs from database and schedule records.

FY 2015/09/29:
    -- Schedule each filter seperately
    
    -- Rely on database for lastobs and quality (won't be in db if WCS fails)

"""

import os
import re
import sys
import glob
import ephem
import argparse
import subprocess
import numpy as np
import cPickle as pickle
from Utils import Constants
from Utils.Catalog import match_coords
from Utils.KS91Moon import KS91Moon
import mydjango.jobstats.models as js
import mydjango.followup.models as fu

# ----------------------------------------------------------------------------
#                 Important constants global to the script
# ----------------------------------------------------------------------------


# format string for writing start/stop times in scheduler lines
timefmt = "%H:%M:%S.0"

# Siding Spring Observatory information
sso = Constants.sso

# Today's date
tomorrow_utc = ephem.now()
today_utc = ephem.Date(ephem.now() - 1.0)
tomorrow_str = tomorrow_utc.datetime().strftime("%Y-%m-%d")
today_str = today_utc.datetime().strftime("%Y-%m-%d")

# Scheduler files
schedpath = os.path.join(Constants.PipelinePath.raw, "search_schedules")
schlogpath = os.path.join(Constants.PipelinePath.raw, "scheduler_logs")
if not os.path.exists(schedpath):
    os.makedirs(schedpath)
search_schedfile = os.path.join(schedpath, "bad_schedule.txt")
follow_schedfile = os.path.join(schedpath, "sn_schedule.txt")


# ----------------------------------------------------------------------------
#                              The main routine
# ----------------------------------------------------------------------------


def read_schedule_fields(qname, qfname, cadence, pmult, minalt):
    """
    Read in a list of fields to observe from a schedule file.
    This effectively initializes a queue for observations with a common
    cadence, priority multipler, and so forth.

    Inputs
        qname:  name of the queue (any string, e.g. "Follow-up")
        qfname:  name of text file on disk to read field centers from
        cadence:  desired cadence of observations, in days
        pmult:  priority multiplier to apply; a crude way to ensure
            follow-up fields are done preferentially, for example
    """

    fields = Constants.SkymapperFieldCenter.read_ascii_file(qfname)
    for f in fields:
        f.set_cadence(cadence, 0.1*cadence)
        f.qname = qname
        f.pmult = pmult
        f.minalt=minalt
    return fields


def read_all_fields():
    """
    Reads in a standard suite of field centers.  Absolute pathnames,
    cadences, and relative priorities of different classes of field centers
    are all hard-wired here until we find a better way to represent them.
    """

    # Read in field centers and assign cadences and priorities.
    fields = [ ]
    etcpath = Constants.PipelinePath.etc
    is_new_field = lambda f: f.id not in [fi.id for fi in fields]
    #add minimum elevation 
    fields += filter(is_new_field, read_schedule_fields(
            "Follow-up",  etcpath + '/skymapper_followup.txt', 2.0, 2.0, 28.0))
    fields += filter(is_new_field, read_schedule_fields(
            "K2 Field 6", etcpath + '/skymapper_kepler6.txt',  3.0, 1.1, 28.0))
    fields += filter(is_new_field, read_schedule_fields(
            "SN Search",  etcpath + '/skymapper_snfields.txt', 4.0, 1.0, 40.0))

    # Read in state information about when fields were last observed.
    state_fn = Constants.PipelinePath.lastobs_state_fn
    if os.path.exists(state_fn):
        field_id, lastobsdate = np.loadtxt(state_fn, unpack=True)
        lastobs = dict(zip(field_id, lastobsdate))
        for f in fields:
            if f.id in lastobs:
                f.set_obs_date(lastobs[f.id])
    return fields


def read_fields_fromdb():
    """
    Read in fields from database.
    One field item for each filter
    """
    #chek active supernova
    #now=datetime.datetime.now()
    #trans=fu.Transient.objects.filter(type__type__in=Constants.follow_types)
    #for tran in trans:
    #    det=tran.pointing_det.all().order_by('dateobs')[0].dateobs
    #    if (now-det).days<60:
    #        continue
        
    sfs=js.SkymapperField.objects.filter(rank__gte=1) #consider all fields being searched
    sfs=list(sfs) #cache the results
    #filters requested
    fields=[]
    for sf in sfs:
        filtexp={}
        for filtinfo in sf.exp.split(","):
            filt, xpt = filtinfo.split(":")
            if filt in filtexp:
                filtexp[filt]['filters'].append(filt)
                filtexp[filt]['exptime']['filt'].append(float(xpt))
            else:
                filtexp[filt]={}
                filtexp[filt]['filters']=[filt]
                filtexp[filt]['exptime']={filt:[float(xpt)]}
        
        for filt in ['u','v','g','r','i','z']: #going through Skymapper filters in order
            if not filt in filtexp:
                continue
            sf.filt=filt
            sf.filters=filtexp[filt]['filters']
            sf.exptime=filtexp[filt]['exptime']
            f=Constants.SkymapperFieldCenter(dict=sf)
            #last observation for this filter
            lastobs= list(sf.sciencepointing_set.filter(program='3PS',filter=filt,status=0).order_by('-dateobs'))
            if len(lastobs)>0:
                f.set_obs_date(ephem.Date(lastobs[0].dateobs))
            fields.append(f)
    
    # set qname
    # convert lastobs to ephem.Date style last_observed
    for f in fields:
        if f.rank==2:
            f.qname+=' Follow-up'
        if f.label==3:
            f.qname+=' K2 Field 6'
        if f.label==4:
            f.qname+=' DES SN Field'
        
    #    f.set_obs_date(ephem.Date(f.lastobs))
        
    # Read in state information about when fields were last observed.
    # in case the database is behind
    #state_fn = Constants.PipelinePath.lastobs_state_fn
    #if os.path.exists(state_fn):
    #    field_id, lastobsdate = np.loadtxt(state_fn, unpack=True)
    #    lastobs = dict(zip(field_id, lastobsdate))
    #    for f in fields:
    #        if f.id in lastobs:
    #            if f.last_observed<lastobs[f.id]-.5: #timezone problem
    #                f.set_obs_date(lastobs[f.id])

    return fields


def dump_lastobs_state(fields):
    """
    Updates information about when fields were last observed.
    The "last observation" information is kept separate from the scheduling
    configuration, because it shouldn't be in svn and it could in principle
    be drawn from a number of sources, including databases.  Here it is kept
    in a simple text file.
        fields:  list of ScheduleField objects to check off
    """

    with open(Constants.PipelinePath.lastobs_state_fn, "w") as outfile:
        for f in fields:
            outfile.write("{0:4d} {1:12.6f}\n".format(f.id, f.last_observed))
    


def update_lastobs_from_scheduler_log(logfname):
    """
    Uses a scheduler log file to mark certain fields as observed.
    """

    # Read in fields from disk
    #fields = read_all_fields()
    fields = read_fields_fromdb()
    # Now read the scheduler log
    manifest = [ ]
    with open(logfname) as logfile:
        for line in logfile:
            fsched = Constants.SchedulerLogEntry(line)
            if fsched.obsid != None:
                manifest.append(fsched)
    # We now have a list of SchedulerLogEntry instances corresponding to
    # what was observed, and a list of ScheduleFields corresponding to which
    # fields we're keeping state for.  Match the two.
    print "# Read {0} fields from scheduler log {1}".format(
            len(manifest), logfname)
    match_idx, best_dr = match_coords(fields, manifest, tol=90.0)
    print "# Matched {0} of {1} fields from scheduler log".format(
            len([i for i in match_idx if i != None]), len(fields))
    print "# fields matched:", [fields[i].id for i in range(len(match_idx))
                              if match_idx[i] is not None]
    # Update the "last observed" bit now.
    print "# Marking as done..."
    for i in range(len(match_idx)):
        if match_idx[i] is not None:
            if fields[i].last_observed<manifest[match_idx[i]].obsdate:
                fields[i].set_obs_date(manifest[match_idx[i]].obsdate)
    # Dump results to disk.
    dump_lastobs_state(fields)


def write_schedule_deadreckon(date, follow_schedfile, search_schedfile,
                              dump_state=False):
    """
    Schedules observations for the *next* available night after given date.
    Writes results to two schedule files, one for high priority to observe
    at any cost, and one for medium priority to observe in bad seeing.

    Inputs
        date:  ephem.Date object telling when to start clock
            (scheduled night will be next sunset after this date)
        follow_schedfile:  Name of text file to hold schedule for
            high-priority search fields to observe in any conditions
        search_schedfile:  Name of text file to hold schedule for
            average-priority search fields to observe in bad seeing
        dump_state:  Update "last observed" state for scheduled fields?
    """

    # Read in fields from disk.
    #fields = read_all_fields()
    fields = read_fields_fromdb()

    # Text to be written to schedule files
    search_schedstr, follow_schedstr = "", ""

    # First figure out when astronomical twilight is.  This will be when
    # the center of the sun is 18 degrees below the horizon.
    # (Our use of next_setting assumes we're running AFTER local sunrise!)
    sso.date = date
    sso.horizon = '-18'
    sun = ephem.Sun()
    sunset = sso.next_setting(sun, use_center=True)
    sunrise = sso.next_rising(sun, use_center=True)
    night_start, night_end = sunset, sunrise

    # RS 2015/07/16:  We now actually predict sky brightness due to moon
    # and use it to influence our scheduling.
    moon = KS91Moon()
    moon.compute(0.5*(sunset + sunrise))
    stat_str = ("# Scheduling on {0} -- moon phase = {1}\n"
                "# sunset = {2}, sunrise = {3}\n".format
                (sso.date, moon.moon_phase, sunset, sunrise))
    search_schedstr += stat_str
    follow_schedstr += stat_str
            
    # Start the clock at effective sunset.
    sso.date = night_start
    for f in fields:
        f.init_ephem_body()
    # mininum time step to forward when no field can be scheduled.
    obslengths=[f.obs_length for f in fields]
    mintimestep=np.min(obslengths)
    # Main event loop
    my_schedule, sched_time, n = [ ], night_start, 1
    while sched_time < night_end:
        schedline = ""
        sso.date = sched_time
        priority_thresh = 50.0
        # Exclude search fields badly affected by moonlight.  Moon flux of
        # 500 nL corresponds to a background level of about 10 counts/sec
        # or a limiting magnitude of about 19.5 in g-band.  For follow-up
        # fields we should, of course, get whatever data we can.
        moon.compute(sso)
        for f in fields:
            f.calc_priority(sso.date, min_alt=f.minalt)
            if (f.priority < 150.0
                    and moon.sky_brightness(f.az, f.alt, k=0.35) > 500):
                f.priority = 0.0
                #print "setting priority to zero due to moon."
        avail_fields = [f for f in fields if f not in my_schedule
                        and f.priority > priority_thresh]
        if len(avail_fields) < 1:
            sched_time = ephem.Date(sched_time + mintimestep)
            #log nothing to observe
            #check_fields=[f for f in fields if f not in my_schedule]
            #check_priority=np.array([f.priority for f in check_fields])
            #check_f=check_fields[check_priority.argmax()]
            #print ephem.Date(sched_time),"max priority:",check_f.id,check_f.priority,check_f.calc_priority(sso.date,verbose=True,min_alt=f.minalt)
            continue
        # Observe the highest-priority field that we haven't done yet.
        priority = np.array([f.priority for f in avail_fields])
        f = avail_fields[priority.argmax()]
        airmass = 1.0/np.sin(f.body.alt)

        # Figure out the best block of time to observe this field.
        # Do this by running a threshold on the observing priority.
        obs_start, obs_end = ephem.Date(sched_time), ephem.Date(sched_time)
        fmax_prio, dt = f.priority, f.obs_length/2. #5*ephem.minute
        while night_start < obs_start and f.priority > priority_thresh:
            obs_start -= dt
            sso.date = obs_start
            f.calc_priority(sso.date,min_alt=f.minalt)
        f.priority = fmax_prio
        while obs_end < night_end and f.priority > priority_thresh:
            obs_end += dt
            sso.date = obs_end
            f.calc_priority(sso.date,min_alt=f.minalt)
        f.priority = fmax_prio
        obs_start, obs_end = ephem.Date(obs_start+dt), ephem.Date(obs_end-dt)
        # observation window can not be smaller than obs_length
        if obs_end-obs_start<f.obs_length:
            sched_time = ephem.Date(sched_time + mintimestep)
            continue
        # Which queue does it go into?
        if f.priority > 150.0:
            prior_id = '3RD_PARTY'
        else:
            prior_id = 'BAD_SEEING'
        # Print a comment about this field.
        obs_length = (obs_end-obs_start)/ephem.minute
        schedline += (
                "# {0} field {1} filter {8} w/priority = {2:.1f} (thresh={3:.1f})\n"
                "# sched_time = {4}, last observed {5}\n"
                "# length of window = {6:.1f} min, airmass = {7:.2f}\n"
                .format(f.qname, f.id, f.priority, priority_thresh,
                        ephem.Date(sched_time), ephem.Date(f.last_observed),
                        obs_length, airmass, f.filt)
            )
        # Mark the field as observed and advance the time.
        f.set_obs_date(sched_time)
        my_schedule.append(f)
        sched_time += f.obs_length
        # Print the schedule file lines.
        for filt in f.filters:
            for exptime in f.exptime[filt]:
                schedline += (
                    "{t0} {t1} {T} {ra},{dec},{f},{rot},{exp},{n},{prog}\n"
                    .format(t0=obs_start.datetime().strftime(timefmt),
                            t1=obs_end.datetime().strftime(timefmt),
                            ra=str(f.body.a_ra),
                            dec=str(f.body.a_dec),
                            exp=exptime,
                            T=prior_id, f=filt, rot=0, n=n, prog=42))
                n += 1
        # Append it to the correct schedule.
        if prior_id == '3RD_PARTY':
            follow_schedstr += schedline
        else:
            search_schedstr += schedline

    # Write schedules to appropriate output files.
    with open(follow_schedfile, 'w') as output:
        output.write(follow_schedstr)
    with open(search_schedfile, 'w') as output:
        output.write(search_schedstr)

    # If we need to dump state, do that last.
    if dump_state:  dump_lastobs_state(fields)


# ----------------------------------------------------------------------------
#                           Command line interface
# ----------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
            description='Write SkyMapper scheduler schedule files')
    parser.add_argument('--dump_state', action='store_true', default=False,
            help='dump state on when fields were last observed')
    parser.add_argument('--mark_done', nargs='+',
            help='names of scheduler logs from which to mark fields as DONE')
    parser.add_argument('--schedule_today_utc', action='store_true', default=False,
            help='automatically generate schedule for next UTC night')
    parser.add_argument('--test', action='store_true', default=False,
            help='perform a simple test run for the scheduler')
    args = parser.parse_args()
    if args.mark_done:
        for schedlog in args.mark_done:
            update_lastobs_from_scheduler_log(schedlog)
    elif args.schedule_today_utc:
        #today_glob_str = "".join(['*.', today_str, '*'])
        #tomorrow_glob_str = "".join(['*.', tomorrow_str, '*'])
        #for schedlog in (
        #        glob.glob(os.path.join(schlogpath, today_glob_str)) +
        #        glob.glob(os.path.join(schlogpath, tomorrow_glob_str))):
        #    update_lastobs_from_scheduler_log(schedlog)
        write_schedule_deadreckon(
                ephem.now(), follow_schedfile, search_schedfile,
                dump_state=args.dump_state)
        search_schedfile_next = os.path.join(schedpath,
                "bad_schedule_{0}.txt".format(tomorrow_str).replace('-',''))
        follow_schedfile_next = os.path.join(schedpath,
                "sn_schedule_{0}.txt".format(tomorrow_str).replace('-',''))
        subprocess.call(['cp', search_schedfile, search_schedfile_next])
        subprocess.call(['cp', follow_schedfile, follow_schedfile_next])
    elif args.test:
        ephem.now = lambda: ephem.Date("2015/06/20 06:00:00")
        Constants.PipelinePath.lastobs_state_fn = os.path.join(
                Constants.PipelinePath.etc, "skymapper_lastobs_test.txt")
        write_schedule_deadreckon(
                ephem.now(), follow_schedfile, search_schedfile,
                dump_state=False)
        search_schedfile_next = os.path.join(schedpath, "bad_schedule_test.txt")
        follow_schedfile_next = os.path.join(schedpath, "sn_schedule_test.txt")
        subprocess.call(['cp', search_schedfile, search_schedfile_next])
        subprocess.call(['cp', follow_schedfile, follow_schedfile_next])

if __name__ == '__main__':
    main()
