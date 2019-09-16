#!/usr/bin/env python3.7
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


def _bi_search(a, x):
    """Binary search for value ``x`` in array ``a``

    Args:
        a (array): The sorted list in which the number x will be searched
        x   (num): The number to be searched

    Returns:
        The position of nearest left neighbor of ``x``
        The position of nearest right neighbor of ``x``
    """

    if any(np.diff(a) < 0):
        raise RuntimeError('Array is not sorted.')

    if not (min(a) <= x <= max(a)):
        raise RuntimeError('Given element outside array range')

    if x in a:
        return list(a).index(x)

    else:
        index = bisect(a, x)
        return [index - 1, index]


class SN91bg(sncosmo.Source):
    """An SNCosmo Source for SN 1991bg-like supernovae"""

    _param_names = param_names_latex = ['x0', 'x1', 'c']

    def __init__(self, min_phase=-18, max_phase=50):
        """An SNCosmo Source for SN 1991bg-like supernovae.

        Parameters for this source include 'x0', 'x1', and 'c'

        Flux for this Source is determined by linearly interpolating for
        stretch and color and then using a 2d spline for phase and wavelength.

        By default, the phase of this source is reduced to match the phase
        range of the SNCosmo Salt2.4 model. Alternative phase ranges can
        be specified, but be aware that the full spectroscopic templates only
        extends from -18 to 100.

        Args:
            min_phase (float): The minimum phase
        """

        super(SN91bg, self).__init__()
        self.name = 'sn91bg'

        # Define initial parameter values
        self._parameters = np.array([1., 1., 0.55])

        # Load 91bg model
        coords, self._template = load_template(min_phase, max_phase)
        (self._stretch, self._color, self._phase, self._wave) = coords
        self._splines = self.get_splines(*coords, self._template)

    @staticmethod
    def get_splines(stretch, color, phase, wave, template):
        """
        Fit splines to a flux template for phase and wavelength

        Args:
            stretch  (ndarray): Stretch coordinates for the flux template
            color    (ndarray): Color coordinates for the flux template
            phase    (ndarray): Phase coordinates for the flux template
            wave     (ndarray): Wave coordinates for the flux template
            template (ndarray): Template of flux values to fit to

        Returns:
            A two dimensional array of RectBivariateSpline with shape
            len(stretch) X len(color)
        """

        splines = np.empty((len(stretch), len(color)), RectBivariateSpline)
        for i, x1 in enumerate(stretch):
            for j, c in enumerate(color):
                splines[i][j] = RectBivariateSpline(
                    phase,
                    wave,
                    template[i][j],
                    kx=3, ky=3)

        return splines

    def _flux(self, phase, wave):
        """Return the flux for a given phase and wavelength

        Flux is determined by linearly interpolating for stretch and color
        and then using a 2d spline for phase and wavelength.

        Args:
            phase (ndarray): A list of days till maximum for the desired flux
            wave  (ndarray): A list of wavelengths for the desired flux

        Returns:
            A 2d array of flux with shape len(phase) X len(wave)
        """

        amplitude, x1, c = self._parameters
        x1_index = _bi_search(self._stretch, x1)
        c_index = _bi_search(self._color, c)
        if isinstance(x1_index, int) and isinstance(c_index, int):
            flux = amplitude * self._splines[x1_index, c_index](phase, wave)

        elif isinstance(x1_index, int):
            flux = [self._splines[x1_index, i](phase, wave) for i in c_index]
            flux = amplitude * np.interp(c, self._stretch[c_index], flux)

        elif isinstance(c_index, int):
            flux = [self._splines[i, c_index](phase, wave) for i in x1_index]
            flux = amplitude * np.interp(x1, self._stretch[c_index], flux)

        else:
            # Linearly interpolate template for current stretch and color
            spline_flux = np.zeros((2, 2, len(phase), len(wave)))
            for i in range(2):
                for j in range(2):
                    spline_flux[i][j] = self._splines[i][j](phase, wave)

            flux = interpn(
                points=[self._stretch[x1_index], self._color[c_index]],
                values=spline_flux,
                xi=[x1, c])[0]

        # Since the spline will extrapolate we enforce bounds
        flux[phase < min(self._phase)] = 0
        flux[phase > max(self._phase)] = 0
        return flux
