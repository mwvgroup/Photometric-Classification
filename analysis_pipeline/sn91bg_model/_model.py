#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module defines the SN91bgSource class, which acts as a custom SNCosmo
source object for 91bg like supernovae.
"""

import os

import numpy as np
import sncosmo
from scipy.interpolate import RectBivariateSpline

FILE_DIR = os.path.abspath(os.path.dirname(__file__))
COMPILED_MODEL_PATH = os.path.join(FILE_DIR, 'template.npy')


def bi_search(a, x):
    """
    Binary search

    Args:
        a (list): The sorted list in which the number x will be searched
        x  (num): The number to be searched

    Returns:
        The position of nearest neighbors of the number x
    """

    if x in a:
        try:
            return a.index(x)
        except AttributeError:
            return a.tolist().index(x)

    if x < a[0] or x > a[-1]:
        raise ValueError('x is out of range')

    left, right = 0, len(a)
    while abs(right - left) > 1:
        m = (left + right) // 2
        if a[m] > x:
            right = m

        if a[m] < x:
            left = m

    return left, right


def linear_interp(x0, x1, f, x):
    """
    Linear interpolation

    Args:
        x0, x1  (num): Grid points
        x       (num): Coordinates of interpolation point
        f      (list): The list contains values at x0 and x1

    Returns:
        Interpolation value
    """
    y0 = f[0]
    y1 = f[1]
    y = y0 + ((x - x0) * y1 - (x - x0) * y0) / (x1 - x0)
    return y


class SN91bgSource(sncosmo.Source):
    _param_names = ['x0', 'x1', 'c']
    param_names_latex = ['x0', 'x1', 'c']  # used in plotting display

    def __init__(self):
        super(SN91bgSource, self).__init__()

        self.name = '91bg model'
        self.version = 'color interpolation'

        # Model spectra and initial guess for stretch, color, and amplitude
        self._flux_values = self._get_91bg_model()[0]
        self._parameters = np.array([1., 1., 0.55])

        # 4-dimension grid points in the model
        self._stretch = np.array([0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25])
        self._color = np.array([0.0, 0.25, 0.5, 0.75, 1])
        self._phase = np.arange(-18, 101, 1)
        self._wave = np.arange(1000, 12001, 10)

        # Creat bi-cubic spline for phase and wavelength using 5 templates
        # (stretch: 0.65, color: [0,0.25,0.5,0.75,1.])
        self._model_flux = np.full(len(self._color), None)
        for i in range(len(self._color)):
            self._model_flux[i] = RectBivariateSpline(
                self._phase,
                self._wave,
                self._flux_values[i],
                kx=3, ky=3)

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
        """
        Return the flux at given phase and wave

        Time and wavelength are given by arguments phase and wave.
        Stretch, color, and amplitude are given by self._parameters

        Args:
            phase (list): A list of days till maximum for the desired flux
            wave  (list): A list of wavelengths for the desired flux

        Returns:
            A 2d array of flux with shape (<len(phase)>, <len(wave)>)
        """

        A, st, c = self._parameters

        # Linearly interpolate template flux by color
        if c in self._color:
            f = self._model_flux[self._color.tolist().index(c)](phase / (st / 0.65), wave)

        else:
            c1, c2 = bi_search(self._color, c)
            y = [self._model_flux[c1](phase / (st / 0.65), wave),
                 self._model_flux[c2](phase / (st / 0.65), wave)]

            f = linear_interp(self._color[c1], self._color[c2], y, c)

        return A * f
