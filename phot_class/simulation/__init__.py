# !/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Generates light-curves using sncosmo and write models results to file.

Usage Example
-------------

>>> from phot_class import simulation
>>>
>>> # Average 91bg stretch and color determined using SIFTO
>>> print(simulation.avg_stretch)
>>> print(simulation.avg_color)
>>>
>>> # Stretch and color covariance determined using SIFTO
>>> print(simulation.covariance)
>>>
>>> # Return random array 5 stretch and color values for 91bg SNe
>>> stretch_array, color_array = bg_stretch_color(5)
>>>
>>> # Simulate parameters for 91bg SNe
>>> param_dicts = simulation.sim_bg_params(zmin=0, zmax=1, tmin=100, tmax=1000)
>>>
>>> # Generate light curves in SDSS band passes
>>> import sncosmo
>>> from phot_class import models
>>>
>>> models.register_sources()
>>> model = sncosmo.Model('sn91bg')
>>> phase_range = (-50, 50)
>>> light_curves = generate_lc(model, phase_range, param_dicts)

Function Documentation
----------------------
"""

from .sncosmo_sims import (AVG_COLOR as avg_color,
                           AVG_STRETCH as avg_stretch,
                           COVARIANCE as covariance,
                           bg_stretch_color, generate_lc, sim_bg_params)

__all__ = ['avg_color', 'avg_stretch', 'covariance',
           'bg_stretch_color', 'sim_bg_params', 'generate_lc']
