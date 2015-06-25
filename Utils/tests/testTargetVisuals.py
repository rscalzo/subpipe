#!/usr/bin/env python

"""
RS 2015/05/07:  Target Visuals Module
"""

from Utils.TargetVisuals import cutout_bounds


# Exercise the routine that locates cutouts on an image

def cutout_bounds_harness(func, params, result):
    assert func(*params, verbose=True) == result

def test_cutout_bounds():
    paramset = [
            (1307, 2475, 50, 50, 2048, 4096), # wholly contained
            (57.2, 99.6, 50, 50, 2048, 4096), # wholly contained
            (   7, 2475, 50, 50, 2048, 4096), # cutoff on lower x edge
            (1307,    5, 50, 50, 2048, 4096), # cutoff on lower y edge
            (   7,    5, 50, 50, 2048, 4096), # cutoff on lower xy corner
            (2036, 2475, 50, 50, 2048, 4096), # cutoff on upper x edge
            (1307, 4079, 50, 50, 2048, 4096), # cutoff on upper y edge
            (2036, 4079, 50, 50, 2048, 4096), # cutoff on upper xy corner
            (2036,    5, 50, 50, 2048, 4096), # cutoff on mixed xy corner
            (   7, 4079, 50, 50, 2048, 4096), # cutoff on mixed xy corner
        ]
    resultset = [
            (( 0, 50,  0, 50), (1282, 1332, 2450, 2500)),
            (( 0, 50,  0, 50), (  32,   82,   75,  125)),
            ((18, 50,  0, 50), (   0,   32, 2450, 2500)),
            (( 0, 50, 20, 50), (1282, 1332,    0,   30)),
            ((18, 50, 20, 50), (   0,   32,    0,   30)),
            (( 0, 37,  0, 50), (2011, 2048, 2450, 2500)),
            (( 0, 50,  0, 42), (1282, 1332, 4054, 4096)),
            (( 0, 37,  0, 42), (2011, 2048, 4054, 4096)),
            (( 0, 37, 20, 50), (2011, 2048,    0,   30)),
            ((18, 50,  0, 42), (   0,   32, 4054, 4096)),
        ]
    for params, result in zip(paramset, resultset):
        yield cutout_bounds_harness, cutout_bounds, params, result
