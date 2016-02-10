#!/usr/bin/env python

import time

class ProfilerStopwatch(object):
    """Time counter class to help us optimize performance
    
    A quick walltime performance counter.  Can run multiple clocks at the
    same time for different classes of things.
    """

    _default_tag = 'DEFAULT'

    def __init__(self):
        """Initialize internal state"""
        self.reset_all()

    def start(self, tag=_default_tag):
        """Start the stopwatch for a given chunk"""
        now = time.time()
        if tag not in self._start:
            self._start[tag] = now

    def stop(self, tag=_default_tag):
        """Stop the stopwatch for a given chunk"""
        now = time.time()
        if tag in self._start:
            if tag not in self._chunks:
                self._chunks[tag] = [ ]
            self._chunks[tag].append((self._start[tag], now))
            del self._start[tag]

    def reset_all(self):
        """Resets all stopwatches"""
        self._start, self._chunks = { }, { }

    def profile_call(self, func, *args, **kwargs):
        """Wraps and times a function call"""
        tag = "call(s) to {}".format(func.__name__)
        self.start(tag)
        results = func(*args, **kwargs)
        self.stop(tag)
        return results

    def report(self):
        """Report on whatever's still running"""
        total = dict([(tag, sum([c[1] - c[0] for c in self._chunks[tag]]))
                       for tag in self._chunks])
        maxlen = max([len(tag) for tag in self._chunks])
        print "Time used so far:"
        for key in sorted(total.keys(), key=lambda k: -total[k]):
            print "    {:.3f} sec  ".format(total[key]),
            if key[:10] == "call(s) to":
                print "{:d}".format(len(self._chunks[key])),
            print key
