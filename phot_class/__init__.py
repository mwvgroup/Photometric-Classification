#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``phot_class`` package applies the photometric fitting technique from
González-Gaitan et al. 2014 to identify peculiar Type Ia Supernovae (SNe Ia).
It provides multiple light curve fitting routines, various utilities for
parsing and manipulating data, and functions for tabulating fitting and
classification results for an entire photometric survey. For usage examples,
please see the individual documentation of each module.

Modules
-------

+--------------------+--------------------------------------------------------+
| Module Name        | Description                                            |
+====================+========================================================+
| ``classification`` | Tabulates light-curve fits for a given survey and      |
|                    | determines the corresponding classification            |
|                    | coordinates.                                           |
+--------------------+--------------------------------------------------------+
| ``fit_func_wraps`` | Provides wrapped versions of ``sncosmo`` minimization  |
|                    | routines to address bugs.                              |
+--------------------+--------------------------------------------------------+
| ``fitting``        | Runs a series of fits on individual light-curves and   |
|                    | tabulates the results.                                 |
+--------------------+--------------------------------------------------------+
| ``fom``            | Calculates the Figure of Merit (FOM) parameter.        |
+--------------------+--------------------------------------------------------+
| ``models``         | Defines custom ``Source`` classes for modeling         |
|                    | supernovae with sncosmo.                               |
+--------------------+--------------------------------------------------------+
| ``simulation``     | Simulates light-curves using ``sncosmo``.              |
+--------------------+--------------------------------------------------------+
| ``spectra``        | tabulates the properties of spectral features.         |
+--------------------+--------------------------------------------------------+
| ``utils``          | A collection of general utilities.                     |
+--------------------+--------------------------------------------------------+
"""

from pathlib import Path as _Path

dust_path = _Path(__file__).resolve().parent / 'schlegel98_dust_map'
