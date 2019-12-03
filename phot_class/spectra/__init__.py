#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module tabulates the properties of spectral features such as area,
velocity, and  equivalent width. All functions in this module are built to
support ``uarray`` objects from the ``uncertainties`` package as inputs.
"""

from .calc_properties import (bin_spectrum, correct_extinction, dust_map,
                              feature_area, feature_pew, feature_velocity,
                              find_peak_wavelength, guess_feature_bounds,
                              line_locations)
from .graphical_interface import (SpectrumInspector,
                                  tabulate_spectral_properties)

__all__ = (
    'bin_spectrum',
    'correct_extinction',
    'dust_map',
    'feature_area',
    'feature_pew',
    'feature_velocity',
    'find_peak_wavelength',
    'guess_feature_bounds',
    'line_locations',
    'SpectrumInspector',
    'tabulate_spectral_properties'
)
