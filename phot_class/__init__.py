#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``phot_class`` package applies the photometric fitting technique from
González-Gaitan et al. 2014 to identify peculiar Type Ia Supernovae (SNe Ia).
It provides multiple light curve fitting routines, various utilities for
parsing and manipulating data, and functions for tabulating fit/classification
results for an entire photometric survey.
"""

# Todo: Include an example(s)

from pathlib import Path as _Path

dust_path = _Path(__file__).resolve().parent / 'schlegel98_dust_map'
