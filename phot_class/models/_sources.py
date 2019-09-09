#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module defines the custom SNCosmo Sources for 91bg like supernovae."""

import os
from bisect import bisect

import numpy as np
import sncosmo
from scipy.interpolate import RectBivariateSpline, interpn

FILE_DIR = os.path.abspath(os.path.dirname(__file__))
COMPILED_MODEL_PATH = os.path.join(FILE_DIR, 'template.npy')


def load_template(phase_min=None, phase_max=None):
    """Load template spectra for SN 1991bg

    The full template spans 5 color values, 119 phases, and 1101 wavelengths.
    The phase range can be limited by specifying the phase_min and
    phase_max arguments. Returned grid coordinates are ordered as color,
    phase, and wavelength.

    Args:
         phase_min (float): Minimum phase of model
         phase_max (float): Maximum phase of model

    Returns:
        A tuple of grid coordinates for the modeled flux
        An array of modeled flux values
    """

    phase_min = phase_min or -float('inf')
    phase_max = phase_max or float('inf')

    # 4-dimension grid points in the model
    stretch = np.array([0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25])
    color = np.array([0.0, 0.25, 0.5, 0.75, 1])
    phase = np.arange(-18, 101, 1)
    wave = np.arange(1000, 12001, 10)

    # Read model and get flux
    modeled_flux = np.load(COMPILED_MODEL_PATH)

    # Limit model phase range
    phase_indices = ((phase_min <= phase) & (phase <= phase_max))
    modeled_flux = modeled_flux[:, :, phase_indices, :]
    phase = phase[phase_indices]

    return (stretch, color, phase, wave), modeled_flux


# noinspection PyPep8Naming, PyUnusedLocal
def SN91bg(name=None, version='phase_limited'):
    """Return a SN 1991bg-like source class for SNCosmo

    Versions include: 'phase_limited', 'color_interpolation'

    Args:
         name   (None): A dummy argument for compatibility with SNCosmo
         version (str): The version of the template to load
    """

    if version == 'phase_limited':
        return PhaseLimited()

    if version == 'color_interpolation':
        return ColorInterpolation()

    else:
        raise ValueError(f"Unidentified version: '{version}'.")


class PhaseLimited(sncosmo.Source):
    """An SNCosmo Source for SN 1991bg-like supernovae"""

    _param_names = param_names_latex = ['x0', 'x1', 'c']

    def __init__(self, min_phase=None, max_phase=None):
        """An SNCosmo Source for SN 1991bg-like supernovae.

        Parameters for this source include 'x0', 'x1', and 'c'

        Flux for this Source is determined by linearly interpolating for
        stretch and color and then using a 2d spline for phase and wavelength.

        By default, the phase of this source is reduced to match the phase
        range of the SNCosmo Salt2.4 model. Alternative phase ranges can
        be specified, but the spectroscopic templates for this source only
        extend from -18 to 101.

        Args:
            min_phase (float): The minimum phase
        """

        super(PhaseLimited, self).__init__()

        salt2_phase = sncosmo.get_source('salt2', version='2.4')._phase
        min_phase = min_phase or min(salt2_phase)
        max_phase = max_phase or max(salt2_phase)

        # Define version info and initial parameter values
        self.name = 'sn91bg'
        self.version = 'phase_limited'
        self._parameters = np.array([1., 1., 0.55])

        # Get phase limited 91bg model
        (self._stretch, self._color, self._phase, self._wave), self._template \
            = load_template(min_phase, max_phase)

    def _flux(self, phase, wave):
        """Return the flux for a given phase and wavelength

        Flux is determined by linearly interpolating for stretch and color
        and then using a 2d spline for phase and wavelength.

        Args:
            phase (ndarray): A list of days till maximum for the desired flux
            wave  (ndarray): A list of wavelengths for the desired flux

        Returns:
            A 2d array of flux with shape (<len(phase)>, <len(wave)>)
        """

        amplitude, stretch, color = self._parameters

        # Linearly interpolate template for current stretch and color
        interp_flux = interpn(
            points=[self._stretch, self._color],
            values=self._template,
            xi=[stretch, color])

        # Fit a spline in phase and wavelength space
        spline = RectBivariateSpline(
            self._phase,
            self._wave,
            interp_flux[0],
            kx=3, ky=3)

        return amplitude * spline(phase, wave)


class ColorInterpolation(sncosmo.Source):
    """An SNCosmo Source for SN 1991bg-like supernovae"""

    _param_names = param_names_latex = ['x0', 'x1', 'c']

    def __init__(self):
        """An SNCosmo Source for SN 1991bg-like supernovae.

        Parameters for this source include 'x0', 'x1', and 'c'

        """

        super(ColorInterpolation, self).__init__()

        self.name = 'sn91bg'
        self.version = 'color_interpolation'

        # Model spectra and initial guess for stretch, color, and amplitude
        coords, flux = load_template()
        self._template = flux[0]
        self._parameters = np.array([1., 1., 0.55])

        # 4-dimension grid points in the model
        self._stretch, self._color, self._phase, self._wave = coords

        # Creat bi-cubic spline for phase and wavelength using 5 templates
        # (stretch: 0.65, color: [0,0.25,0.5,0.75,1.])
        self._model_flux = np.full(len(self._color), None)
        for i in range(len(self._color)):
            self._model_flux[i] = RectBivariateSpline(
                self._phase,
                self._wave,
                self._template[i],
                kx=3, ky=3)

    @staticmethod
    def bi_search(a, x):
        """Binary search for value ``x`` in array ``a``

        Args:
            a (ndarray): The sorted list in which the number x will be searched
            x     (num): The number to be searched

        Returns:
            The position of nearest left neighbor of ``x``
            The position of nearest right neighbor of ``x``
        """

        index = bisect(a, x)
        return index - 1, index

    @staticmethod
    def linear_interp(x0, x1, f, x):
        """Linear interpolation

        Args:
            x0, x1  (num): Grid points
            x       (num): Coordinates of interpolation point
            f      (list): The list contains values at x0 and x1

        Returns:
            Interpolated value
        """

        return f[0] + ((x - x0) * f[1] - (x - x0) * f[0]) / (x1 - x0)

    def _flux(self, phase, wave):
        """Return the flux for a given phase and wavelength

        Flux is determined by adjusting the phase for the given stretch value
        and then using a spline to determine the template for the adjusted
         phase and wavelength. The resulting values are then linearly
         interpolated for color.

        Args:
            phase (ndarray): A list of days till maximum for the desired flux
            wave  (ndarray): A list of wavelengths for the desired flux

        Returns:
            A 2d array of flux with shape (<len(phase)>, <len(wave)>)
        """

        A, st, c = self._parameters

        # Linearly interpolate template flux by color
        if c in self._color:
            f = self._model_flux[self._color.tolist().index(c)](
                phase / (st / 0.65), wave)

        else:
            c1, c2 = self.bi_search(self._color, c)
            y = [self._model_flux[c1](phase / (st / 0.65), wave),
                 self._model_flux[c2](phase / (st / 0.65), wave)]

            f = self.linear_interp(self._color[c1], self._color[c2], y, c)

        return A * f
