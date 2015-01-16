#!/usr/bin/env python

# ============================================================================
# RS 2012/02/26:  Exceptions for externally called pipeline routines
# ----------------------------------------------------------------------------
# These allow us to report failure conditions in the subtraction pipeline in
# what I hope is a sane and uniform way.  We'll see how we go.
# ============================================================================

import os


class TrackableException(Exception):
    """A base exception class for failures to be tracked in the pipeline."""
    def __init__ (self, msg="", **kwargs):
        self.msg = msg
        for kw in kwargs:  setattr(self,kw,kwargs[kw])
    def __str__ (self):
        return self.msg


class ExternalFailure(TrackableException):
    """A class specifically meant for failures of external programs."""
    def __init__ (self, cmd="???", exit_code=0):
        msg = "{0} failed with exit code {1}".format(
               os.path.basename(cmd.split()[0]), exit_code)
        super(ExternalFailure, self).__init__(
              msg, cmd=cmd, exit_code=exit_code)
