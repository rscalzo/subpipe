#!/bin/tcsh

# ============================================================================
# RS 2011/12/05:  Pipeline Wrapper for Production Running
# ----------------------------------------------------------------------------
# This is just a wrapper around subpipe_master.py which clears the environment
# before running.  This will allow us to run cleanly under group accounts like
# the skymapper account, without worrying about running in an environment with
# a cluttered $PATH or other confusing things.
# ============================================================================

# Source all the pipeline environment definitions, then run the pipeline.
# Save a few environment settings which we know we'll need.
/usr/bin/env - SUBPIPEHOME=/export/maipenrai/skymap/skymap/subpipe HOME=$HOME SHELL=$SHELL \
    tcsh -c 'source $SUBPIPEHOME/subpipe_env.csh; $SUBPIPEHOME/subpipe_master.py'
    # tcsh -c 'source $SUBPIPEHOME/subpipe_env.csh; setenv; which swarp; swarp -v'
