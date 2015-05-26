#! /bin/bash

# ============================================================================
# RS 2011/12/05:  Clean Pipeline Environment
# ----------------------------------------------------------------------------
# These settings should be sourced so that the pipeline knows where to find
# all its components.  It's automatically sourced within a clean environment
# by running run_subpipe.csh, in this directory.
# ============================================================================

# If we're starting with a clean environment, we'll need to set these
# variables to keep run_subpipe.csh from grinding to a halt.

# All external packages are installed under $SUBPIPEHOME/ext/.  On maipenrai
# this is a soft link to:  /export/maipenrai/skymap/skymap/external_packages

export SUBPIPEEXT=${SUBPIPEHOME}/ext

# ----------------------------------------------------------------------------
#                The baseline production pipeline environment
# ----------------------------------------------------------------------------

# By default, this is a production install and we should be using the real
# Apache webserver, the real Postgres backend to Django, etc.
export SUBPIPE_PRODUCTION_INSTALL=1

# Paths that should always be in $PATH.  Make sure $LOCAL/bin is first,
# we may occasionally override system defaults there.
export PATH=${SUBPIPEEXT}/bin:${PATH}
# Home path for pipeline code
export SUBETCPATH=$SUBPIPEHOME/etc
export PATH=$SUBPIPEHOME/bin:${PATH}
# Subpath(s) for storing data associated with a particular pipeline run
export SUBPIPEDATA=$SUBPIPEHOME/data
export SUBRAWPATH=$SUBPIPEDATA/raw
export SUBNEWPATH=$SUBPIPEDATA/new
export SUBREDPATH=$SUBPIPEDATA/red
export SUBSUBPATH=$SUBPIPEDATA/sub
export SUBCALPATH=$SUBPIPEDATA/cal
export SUBREFPATH=$SUBPIPEDATA/ref
# Subpath(s) for storing metadata for the Django-powered data browsers
export MEDIAPATH=$SUBPIPEHOME/media
export SUBTHUMBS=$MEDIAPATH/thumbs
export LCFIGPATH=$MEDIAPATH/lc_figs
export SUBDBPATH=$SUBPIPEHOME/db
# FY:  path to candidate files
export CANDDAILY=$SUBPIPEDATA/cand_daily
export CANDPATH=$SUBPIPEDATA/candidates
export CANDLOOKUP=$SUBPIPEDATA/cand_lookup.fits
# Python environment -- include access to system packages, but override with
# anything we've installed locally that's important for us.
export PYTHONPATH=${PYTHONPATH}:${SUBPIPEHOME}:${SUBPIPEHOME}/STAP
# setenv PYTHONPATH ${PYTHONPATH}:/usr/local/python-2.7.1/lib/python2.7/site-packages/
# setenv PYTHONPATH ${PYTHONPATH}:${SUBPIPEEXT}/lib/python2.7/site-packages
# Source virtual python environment (if needed)
#source $SUBPIPEEXT/bin/activate.sh
# ... Brian's astrometry module ...........................................
# cdsclient binaries, which Brian's stuff needs to run
export PATH=${PATH}:${SUBPIPEEXT}/bin/cdsclient
# Brian's actual stuff
export BRIANWCS=${SUBPIPEHOME}/WCS
export brianwcs=${SUBPIPEHOME}/WCS

source $BRIANWCS/etc/brianwcs.sh
# ... external packages ...................................................
# cfitsio
export CFITSIO=${SUBPIPEEXT}
# DSS image archive interface
export DSS_ROOT=${SUBPIPEEXT}/src/dss/
# django database access
export DJANGO_SETTINGS_MODULE="mydjango.settings"

# ----------------------------------------------------------------------------
#                Load user-specific settings incl. overrides
# ----------------------------------------------------------------------------
# If there are environment settings specifically for this user, source them
# here.  This way one can, for example, override the production pipeline
# settings (always in this file) with his/her own test settings as desired.
# ----------------------------------------------------------------------------

export USER_SUBENV=$SUBPIPEHOME/subpipe_env_$USER.sh
if [ -e $USER_SUBENV ]
then
    echo Sourcing user settings in $USER_SUBENV ;
    source $USER_SUBENV ;
fi
