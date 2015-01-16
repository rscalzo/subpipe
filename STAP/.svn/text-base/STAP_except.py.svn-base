#!/usr/bin/env python2.6

# ============================================================================
# RS 2012/02/26:  Exceptions for STAP_* subtraction pipeline routines
# ----------------------------------------------------------------------------
# These allow us to report failure conditions in the subtraction pipeline in
# what I hope is a sane and uniform way.  We'll see how we go.
# ============================================================================


class STAP_Error(Exception):
    def __init__ (self, msg="", **kwargs):
        self.msg = msg
        for kw in kwargs:  setattr(self,kw,kwargs[kw])
    def __str__ (self):
        return self.msg

class BadInputError(STAP_Error):
    pass

class ExternalFailure(STAP_Error):
    def __init__ (self, cmd="", exit_code=0):
        # TODO:  Figure out why normal positional arguments don't work here.
        cmdbin = cmd.split(" ")[0]
        msg = "{0} failed with exit code {1}".format(cmdbin, exit_code)
        super(ExternalFailure,
              self).__init__(msg, cmd=cmd, exit_code=exit_code)
