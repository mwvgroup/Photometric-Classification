#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module defines the SN91bgSource class, which acts as a custom SNCosmo
source object for 91bg like supernovae.
"""

import os

import numpy as np
import sncosmo
from scipy.interpolate import interpn

FILE_DIR = os.path.abspath(os.path.dirname(__file__))
COMPILED_MODEL_PATH = os.path.join(FILE_DIR, 'complete_template.npy')


class SN91bgSource(sncosmo.Source):
    _param_names = ['stretch', 'color', 'amplitude']
    param_names_latex = ['x1', 'c', 'a']

    def __init__(self):
        super(SN91bgSource, self).__init__()

        self.name = '91bg model'
        self.version = 'None'

        # Model spectra and initial guess for stretch, color, and amplitude
        self._flux_values = self._get_91bg_model()
        self._parameters = np.array([1, .5, 1])

        # Create an array of xi points for self._flux_values
        self._stretch_vals = np.arange(0.65, 1.26, .1)
        self._color_vals = np.array([0.0, 0.25, 0.5, 0.75, 1])
        self._phase = np.arange(-13, 101, 1)
        self._wave = np.arange(1000, 12001, 10)
        self._flux_points = [self._stretch_vals,
                             self._color_vals,
                             self._phase,
                             self._wave]

    @staticmethod
    def _get_91bg_model():
        """Load template spectra for SN 1991bg

        Template spans 7 stretch values, 5 color values, 114 dates,
        and 1101 wavelengths.

        Returns:
            An array of shape (7, 5, 114, 1101)
        """

        return np.load(COMPILED_MODEL_PATH)

    def _flux(self, phase, wave):
        """Interpolate template flux for stretch, color, time, and wavelength

        Time and wavelength are given by arguments phase and wave. Stretch,
        color, and amplitude are given by self._parameters/

        Args:
            phase (list): A list of days till maximum for the desired flux
            wave  (list): A list of wavelengths for the desired flux

        Returns:
            A 2d array of flux with shape (<len(phase)>, <len(wave)>)
        """

        stretch, color, amplitude = self._parameters
        requested_xi = np.empty((len(phase), len(wave), 4))
        for i, p in enumerate(phase):
            for j, w in enumerate(wave):
                requested_xi[i, j] = [stretch, color, p, w]

        print('calling interpolation')
        flux = interpn(points=self._flux_points,
                       values=self._flux_values,
                       xi=requested_xi)

        return amplitude * flux
