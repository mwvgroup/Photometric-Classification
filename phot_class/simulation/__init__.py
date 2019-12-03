# !/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Generates light-curves using sncosmo and write models results to file.

Usage Example
-------------



Function Documentation
----------------------
"""

from .sncosmo_sims import (AVG_COLOR, AVG_STRETCH, COVARIANCE,
                           bg_stretch_color, generate_lc, sim_bg_params)

__all__ = ['AVG_COLOR', 'AVG_STRETCH', 'COVARIANCE',
           'bg_stretch_color', 'sim_bg_params', 'generate_lc']
