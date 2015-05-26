# ============================================================================
# RS 2011/12/05:  Clean Pipeline Environment
# ----------------------------------------------------------------------------
# These settings should be sourced so that the pipeline knows where to find
# all its components.  It's automatically sourced within a clean environment
# by running run_subpipe.csh, in this directory.
# ============================================================================

if ( ! $?SUBPIPEHOME ) then
   setenv SUBPIPEHOME `pwd`
   echo Using SUBPIPEHOME=${SUBPIPEHOME}
endif

# If we're starting with a clean environment, we'll need to set these
# variables to keep run_subpipe.csh from grinding to a halt.

if ( ! $?PATH ) then
   setenv PATH .
endif
if ( ! $?PERL5LIB ) then
   setenv PERL5LIB .
endif
if ( ! $?PYTHONPATH ) then
   setenv PYTHONPATH .
endif

# All external packages are installed under $SUBPIPEHOME/ext/.  On maipenrai
# this is a soft link to:  /export/maipenrai/skymap/skymap/external_packages

setenv SUBPIPEEXT ${SUBPIPEHOME}/ext

# ----------------------------------------------------------------------------
#                The baseline production pipeline environment
# ----------------------------------------------------------------------------

# By default, this is a production install and we should be using the real
# Apache webserver, the real Postgres backend to Django, etc.
setenv SUBPIPE_PRODUCTION_INSTALL 1

# Paths that should always be in $PATH.  Make sure $LOCAL/bin is first,
# we may occasionally override system defaults there.
setenv PATH ${SUBPIPEEXT}/bin:${PATH}:/bin:/usr/bin:/usr/local/bin:.
# Home path for pipeline code
setenv SUBETCPATH  $SUBPIPEHOME/etc
setenv PATH $SUBPIPEHOME/bin:${PATH}
# Subpath(s) for storing data associated with a particular pipeline run
setenv SUBPIPEDATA $SUBPIPEHOME/data
setenv SUBRAWPATH  $SUBPIPEDATA/raw
setenv SUBNEWPATH  $SUBPIPEDATA/new
setenv SUBREDPATH  $SUBPIPEDATA/red
setenv SUBSUBPATH  $SUBPIPEDATA/sub
setenv SUBCALPATH  $SUBPIPEDATA/cal
setenv SUBREFPATH  $SUBPIPEDATA/ref
# Subpath(s) for storing metadata for the Django-powered data browsers
setenv MEDIAPATH   $SUBPIPEHOME/media
setenv SUBTHUMBS   $MEDIAPATH/thumbs
setenv LCFIGPATH   $MEDIAPATH/lc_figs
setenv SUBDBPATH   $SUBPIPEHOME/db
# FY:  path to candidate files
setenv CANDDAILY $SUBPIPEDATA/cand_daily
setenv CANDPATH $SUBPIPEDATA/candidates
setenv CANDLOOKUP $SUBPIPEDATA/cand_lookup.fits
# Local scratch disk for temporary calculations
setenv SUBSCRATCH  /ramdisk/subpipe_production
# Python environment -- include access to system packages, but override with
# anything we've installed locally that's important for us.
setenv PYTHONPATH ${PYTHONPATH}:${SUBPIPEHOME}:${SUBPIPEHOME}/STAP
# setenv PYTHONPATH ${PYTHONPATH}:/usr/local/python-2.7.1/lib/python2.7/site-packages/
# setenv PYTHONPATH ${PYTHONPATH}:${SUBPIPEEXT}/lib/python2.7/site-packages
source $SUBPIPEEXT/bin/activate.csh
# ... Brian's astrometry module ...........................................
# cdsclient binaries, which Brian's stuff needs to run
setenv PATH ${PATH}:${SUBPIPEEXT}/bin/cdsclient
# Brian's actual stuff
setenv BRIANWCS ${SUBPIPEHOME}/WCS
source $BRIANWCS/etc/brianwcs.csh
# ... external packages ...................................................
# cfitsio
setenv CFITSIO ${SUBPIPEEXT}
# astrometry.net global install
setenv PATH ${PATH}:/usr/local/astrometry/bin/
setenv TMPDIR $SUBSCRATCH/astronet
# DSS image archive interface
setenv DSS_ROOT ${SUBPIPEEXT}/src/dss/
# django database access
setenv DJANGO_SETTINGS_MODULE mydjango.settings

# ----------------------------------------------------------------------------
#                Load user-specific settings incl. overrides
# ----------------------------------------------------------------------------
# If there are environment settings specifically for this user, source them
# here.  This way one can, for example, override the production pipeline
# settings (always in this file) with his/her own test settings as desired.
# ----------------------------------------------------------------------------

setenv USER_SUBENV $SUBPIPEHOME/subpipe_env_$USER.csh
if ( -e $USER_SUBENV ) then
    echo Sourcing user settings in $USER_SUBENV
    source $USER_SUBENV
endif
