#!/usr/bin/env csh
# Environment configuration for Brian's astrometry routines
# Required:  set environment variable BRIANWCS
if ( ! $?BRIANWCS ) then
    echo "PROBLEM:  environment variable BRIANWCS not set."
    echo "Please go and set it, I can't function without it!"
else
    if ( ! $?PERL5LIB ) then
        setenv PERL5LIB ${BRIANWCS}/lib/perl5
    else
        setenv PERL5LIB ${PERL5LIB}:${BRIANWCS}/lib/perl5
    endif
    setenv PERL5LIB ${PERL5LIB}:${BRIANWCS}/lib/perl5
    setenv PATH ${PATH}:${BRIANWCS}/bin
endif

# Required:  cdsclient should be in the user's $PATH
if ! { set tmpvar = `which vizquery` } then
    echo "PROBLEM:  Can't find cdsclient binaries; please add to your PATH."
endif
