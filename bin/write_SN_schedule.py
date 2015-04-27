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
"""

import os
import re
import sys
import glob
import ephem
import argparse
import numpy as np
import cPickle as pickle
from Utils import Constants
from Utils.Catalog import match_coords

# ----------------------------------------------------------------------------
#                 Important constants global to the script
# ----------------------------------------------------------------------------


# format string for writing start/stop times in scheduler lines
timefmt = "%H:%M:%S.0"

# Siding Spring Observatory information
sso = Constants.sso

# Global sets of fields for useful commentary purposes
fieldidsets = { }


# ----------------------------------------------------------------------------
#                              The main routine
# ----------------------------------------------------------------------------


def read_schedule_fields():
    """Reads in fields, and their state, from the standard configuration"""
    # Read in fields to observe from disk.  Read straight from text files,
    # as these may be modified at any time (particularly follow-up fields).
    fields, fieldsets = [ ], { }
    for fsetname, fname, cadence, cwidth in [
            ("SkyMapper",   "skymapper_snfields.txt",   0.5, 0.2),
            # ("SkyMapper",   "skymapper_snfields_snzoo.txt",   0.5, 0.1),
            # ("SkyMapper",   "skymapper_20150201_hacked.txt",   4.0, 0.2),
            # ("SkyMapper",   "skymapper_snfields_cosmology.txt",   3.0, 1.0),
            ("Follow-up",   "skymapper_followup.txt",   1.0, 0.1),
            ("LSQ overlap", "skymapper_lsqoverlap.txt", 3.0, 1.0),
            ("K1 overlap",  "skymapper_kepler1.txt",    2.0, 0.5),
            ("K3 overlap",  "skymapper_kepler3.txt",    3.0, 0.2),
            ]:
        absfn = Constants.PipelinePath.etc + '/' + fname
        ftmp = Constants.SkymapperFieldCenter.read_ascii_file(absfn)
        for f in ftmp:  f.set_cadence(cadence, cwidth)
        # HACK:  just observe SN Zoo fields and follow-up fields for now
        if len(fields) < 2:  fields = fields + ftmp
        # fields = fields + ftmp
        fieldidsets[fsetname] = [f.id for f in ftmp]
    SKMfields, LSQfields, TRGfields, K1fields, K3fields = fieldidsets
    # Read in state information about when fields were last observed.
    state_fn = Constants.PipelinePath.lastobs_state_fn
    if os.path.exists(state_fn):
        field_id, lastobsdate = np.loadtxt(state_fn, unpack=True)
        lastobs = dict(zip(field_id, lastobsdate))
        for f in fields:
            if f.id in lastobs:
                f.set_obs_date(lastobs[f.id])
    return fields


def dump_lastobs_state(fields):
    """Updates information about when fields were last observed
    
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
    """Uses a scheduler log file to mark certain fields as observed."""

    # Read in fields from disk
    fields = read_schedule_fields()
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
            fields[i].set_obs_date(manifest[match_idx[i]].obsdate)
    # Dump results to fisk.
    dump_lastobs_state(fields)


def write_schedule_singleobs(verbose=False, dump_state=False):
    """Generates a single OB request for the SkyMapper scheduler."""

    # Read in fields from disk
    fields = read_schedule_fields()
    sso.date = ephem.now()
    # Check what's up at this time that we haven't observed.
    # Pick the highest-priority field that's up and observe it.
    for f in fields:
        f.init_ephem_body()
        f.calc_priority(sso.date)
    priority = np.array([f.priority for f in fields])
    f = fields[priority.argmax()]
    # In "single-obs" mode, we'll need to save the state of the program
    # because Patrick only wants to make one call at a time.
    filt_idx = 0
    schfname = Constants.PipelinePath.etc + "/schedule_single_obs.txt"
    statepkl = Constants.PipelinePath.etc + "/schedule_single_obs.pkl"
    # If a state file exists, load the state.
    # What we need to know:  field, n, next filter.
    if os.path.exists(statepkl):
        with open(statepkl) as infile:
            tmpf, n, filt_idx, obdate = pickle.load(infile)
        # If we have any filters left to do, observe the same field.
        # Otherwise, start the first filter of whatever field we
        # picked in the shenanigans above.
        if filt_idx < len(f.filters):  f = tmpf
        else:  filt_idx = 0
        # See if this observation is for a new night
        if obdate < sso.date - 0.33:  n = 0
    else:
        n, filt_idx = 0, 0
    # Make a few conversions.
    filt = f.filters[filt_idx]
    ra = ephem.hours(f.ra*ephem.pi/180.0)
    dec = ephem.degrees(f.dec*ephem.pi/180.0)
    obdate = ephem.now()
    # Now that the state has been reset, write the schedule line.
    with open(schfname, "w") as outfile:
        outstr = "{ra},{dec},{f},{rot},{exp},{n},{prog}".format(
            ra=str(ra), dec=str(dec), f=filt, exp=f.exptime[filt],
            rot=0, n=n, prog=42)
        print outstr
        outfile.write(outstr + "\n")
    # Finally, dump the pickle state back out.
    filt_idx += 1
    f.set_obs_date(obdate)
    f.body = f.rise = f.set = None
    n += 1
    with open(statepkl, 'w') as outfile:
        pickle.dump([f, n, filt_idx, 1.0*obdate], outfile)
    # If we need to dump overall state, do that last.
    if dump_state:  dump_lastobs_state(fields)


def write_schedule_deadreckon(verbose=False, dump_state=False,
                              start_night=None, stop_night=None,
                              exclude_moon=True):
    """Writes a schedule file to stdout."""

    # Read in fields from disk
    fields = read_schedule_fields()

    # First figure out when astronomical twilight is.  This will be when
    # the center of the sun is 18 degrees below the horizon.
    sso.date = ephem.now()
    sso.horizon = '-18'
    # RS 2013/09/13:  Added optional switches for start and stop of night.
    # Defaults to the entire night.
    sun = ephem.Sun()
    if stop_night:
        sunrise = ephem.Date(str(sso.date).split()[0] + " " + stop_night)
    else:
        sunrise = sso.next_rising(sun, use_center=True)
    if start_night:
        sunset = ephem.Date(str(sso.date).split()[0] + " " + start_night)
    else:
        sunset = sso.next_setting(sun, use_center=True)
    # If we're running this in the middle of a night, "sunset" may already
    # have happened; but we always mean the sunset at the start of the
    # current night of observing.  So make sure sunset < sunrise.
    if sunset > sunrise:
        sunset = sso.previous_setting(sun, use_center=True)
    print "# SN scheduler:  fixed sun stuff"
    # Default:  start the night when the sun goes down
    night_start, night_end = sunset, sunrise

    # RS 2013/11/25:  Add moon calculations.  The goal here is to make sure
    # that if the Moon is more full than X% (we use 50% here), then scheduled
    # observations take place when the Moon is not in the sky.  To do that,
    # we need to calculate the moon rise and set times *most pertinent* for
    # observations on the given UT date.  First see if the Moon is more than
    # half full in the middle of the night.
    moon = ephem.Moon()
    moon.compute(0.5*(sunset + sunrise))
    sso.date = sunset
    moonrise = sso.next_rising(moon, use_center=True)
    moonset = sso.next_setting(moon, use_center=True)
    # If moon is more than half full, we need to see when it's up.
    if moon.moon_phase > 0.5 and exclude_moon:
        # If the moon is *already* up at the start of the night, then the next
        # moon rising will be *after* the next moon setting.  In that case,
        # we should start observations after it sets.
        if sunset < moonset < sunrise < moonrise:
            night_start = moonset
        # If the moon rises during the night, then the next moon setting will
        # be after the next moon rising.  In that case, we should observe
        # *until* the moon comes up.
        elif sunset < moonrise < sunrise < moonset:
            night_end = moonrise
        # If the moon is up all night, don't even observe.
        elif sunset < sunrise < moonset < moonrise:
            night_end = night_start
        # If none of these three situations cover things, then, ummm...
        else:
            print "# WARNING:  moon calculations screwed up, aborting"
            return
    # If the moon is less than half full, we don't care when it's up
    print "# SN scheduler:  fixed moon stuff"
    print "# HACK:  only observing Kepler Field 1 right now,",
    print "see read_schedule_fields() to fix"
    print "# DOUBLE HACK:  disabling moon phase lockout temporarily"
    print "# Scheduling on", sso.date, "-- moon phase =", moon.moon_phase
    print "# sunset   =", sunset,   "sunrise =", sunrise
    if moon.moon_phase > 0.5:
        print "# moonrise =", moonrise, "moonset =", moonset

    # Start the clock at effective sunset.
    sso.date = night_start
    if verbose:  print "# obstart  = {0}, obend  = {1}".format(night_start, night_end)
    if verbose:
        print "# Initializing pyephem FixedBody elements...";
        sys.stdout.flush()
    for f in fields:  f.init_ephem_body()
    print "# Attempting to schedule", len(fields), "fields"
    # Main event loop
    my_schedule, sched_time, n = [ ], night_start, 1
    while sched_time < night_end:
        sso.date = sched_time
        priority_thresh = 0.35*np.max([f.calc_priority(sso.date) for f in fields])
        priority_thresh = 50.0
        # Check what's up at this time that we haven't observed
        avail_fields = [f for f in fields if f not in my_schedule
                        and f.priority > priority_thresh]
        if len(avail_fields) < 1:
            sched_time = ephem.Date(sched_time + fields[0].obs_length)
            continue
        # Do follow-up fields if they're above the priority threshold.
        # Else, pick the highest-priority field that's up and observe it.
        priority = np.array([f.priority for f in avail_fields])
        f = avail_fields[priority.argmax()]

        # Figure out the best block of time to observe this field.
        # Do this by running a threshold on the observing priority.
        obs_start, obs_end = ephem.Date(sched_time), ephem.Date(sched_time)
        fmax_prio, dt = f.priority, ephem.minute
        while night_start < obs_start and f.priority > priority_thresh:
            obs_start -= dt
            sso.date = obs_start
            f.calc_priority(sso.date)
        f.priority = fmax_prio
        while obs_end < night_end and f.priority > priority_thresh:
            obs_end += dt
            sso.date = obs_end
            f.calc_priority(sso.date)
        f.priority = fmax_prio
        obs_start, obs_end = ephem.Date(obs_start+dt), ephem.Date(obs_end-dt)
        # Print a comment if the user has "verbose" set.
        if verbose:
            for fsetname in ("Follow-up", "LSQ overlap",
                             "K1 overlap", "K3 overlap", "SkyMapper"):
                if f.id in fieldidsets[fsetname]:  break
            print "#", fsetname,
            print "SN search field", f.id,
            print "scheduled with priority = {0:.1f}".format(f.priority),
            print "(thresh = {0:.1f})".format(priority_thresh)
            print "# sched_time =", ephem.Date(sched_time),
            print "last observed", ephem.Date(f.last_observed).datetime()
            print "# length of window = {0:.1f} min,".format(
                    (obs_end-obs_start)/ephem.minute),
            print "airmass = {0:.2f}".format(1.0/np.sin(f.body.alt))
        # Mark the field as observed and advance the time.
        f.set_obs_date(sched_time)
        my_schedule.append(f)
        sched_time += f.obs_length
        # Print the schedule file lines.
        for filt in f.filters:
            print "{t0} {t1} {T} {ra},{dec},{f},{rot},{exp},{n},{prog}"\
                  .format(t0=obs_start.datetime().strftime(timefmt),
                          t1=obs_end.datetime().strftime(timefmt),
                          ra=str(f.body.a_ra),
                          dec=str(f.body.a_dec),
                          exp=f.exptime[filt],
                          T="BAD_SEEING", # T="3RD_PARTY",
                          f=filt, rot=0, n=n, prog=42)
            n += 1

    # If we need to dump state, do that last.
    if dump_state:  dump_lastobs_state(fields)


# ----------------------------------------------------------------------------
#                           Command line interface
# ----------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
            description='Write SkyMapper scheduler schedule files')
    parser.add_argument('--verbose', action='store_true', default=False,
            help='print hash comments to standard output')
    parser.add_argument('--single_obs', action='store_true', default=False,
            help='only schedule a single observation')
    parser.add_argument('--dump_state', action='store_true', default=False,
            help='dump state on when fields were last observed')
    parser.add_argument('--mark_done', nargs='+',
            help='names of scheduler logs from which to mark fields as DONE')
    parser.add_argument('--start_night', default=None, # default="13:30:00",
            help='start of SN observing block, UT time (hh[:mm[:ss]])')
    parser.add_argument('--stop_night', default=None, # default="15:30:00",
            help='end of SN observing block, UT time (hh[:mm[:ss]])')
    parser.add_argument('--schedule_today_utc', action='store_true', default=False,
            help='automatically generate schedule for next UTC night')
    args = parser.parse_args()
    if args.mark_done:
        for schedlog in args.mark_done:
            update_lastobs_from_scheduler_log(schedlog)
    elif args.single_obs:
        write_schedule_singleobs(verbose=args.verbose,
                                 dump_state=args.dump_state)
    elif args.schedule_today_utc:
        # today_utc, tomorrow_utc = ephem.now(), ephem.Date(ephem.now() + 1.0)
        today_utc, tomorrow_utc = ephem.Date(ephem.now() - 1.0), ephem.now()
        today_glob_str = today_utc.datetime().strftime("*.%Y-%m-%d*")
        tomorrow_glob_str = tomorrow_utc.datetime().strftime("*.%Y-%m-%d*")
        tomorrow_str = tomorrow_utc.datetime().strftime("%Y%m%d")
        CPP_base_glob = Constants.PipelinePath.raw + "/scheduler_logs/"
        for schedlog in (glob.glob(CPP_base_glob + today_glob_str) +
                         glob.glob(CPP_base_glob + tomorrow_glob_str)):
            update_lastobs_from_scheduler_log(schedlog)
        sched_path = Constants.PipelinePath.raw + "/search_schedules/"
        sched_file = sched_path + "sn_schedule_" + tomorrow_str + ".txt"
        tmpstdout = sys.stdout
        with open(sched_file, 'w') as output:
            sys.stdout = output
            write_schedule_deadreckon(verbose=args.verbose,
                                      dump_state=args.dump_state,
                                      start_night=args.start_night,
                                      stop_night=args.stop_night,
                                      exclude_moon=False)
        sys.stdout = tmpstdout
        with open(sched_file) as whatwejustdid:
            for line in whatwejustdid:
                if line[0] != '#':
                    print line.strip()
    else:
        write_schedule_deadreckon(verbose=args.verbose,
                                  dump_state=args.dump_state,
                                  start_night=args.start_night,
                                  stop_night=args.stop_night,
                                  exclude_moon=False)

if __name__ == '__main__':
    main()
