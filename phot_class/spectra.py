#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module calculates the properties of spectral features such as area,
velocity, and  equivalent width.
"""

import numpy as np
from astropy import units
from astropy.constants import c as speed_of_light
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values, std_devs, uarray

# In Angstroms
line_locations = {
    'CaHK': 3945.02,
    '4130': 4129.78,
    'MgII': 4481.00,
    'FeII': 5169.00,
    'SIIW1': 5449.20,
    'SIIW2': 5622.46,
    '5972': 5971.89,
    '6355': 6356.08,
    'CaII': 8578.79,  # Also 8498, 8542 and 8662
    'NaD': 5895.0,
    'unk1': 7000.0,
    'unk2': 8000.0,
    'FeI': 7000.0,
    'OI': 7774.0,
    'CaNIR': (11784. + 11839. + 11849.) / 3.}


def feature_area(wave, flux):
    """Calculate the area of a feature

    Args:
        wave         (ndarray): A sorted array of wavelengths for the feature
        flux (ndarray, uarray): An array of flux values for each wavelength

    Returns:
        The area of the feature
    """

    # Feature area = area under continuum - area under spectrum
    continuum_area = (wave[-1] - wave[0]) * (flux[0] + flux[-1]) / 2
    spectrum_area = np.trapz(y=flux, x=wave)
    return continuum_area - spectrum_area


def feature_pew(wave, flux):
    """Calculate the pseudo equivalent-width of a feature

    Args:
        wave         (ndarray): A sorted array of wavelengths for the feature
        flux (ndarray, uarray): An array of flux values for each wavelength

    Returns:
        The normalized flux
        The pseudo equivalent-width of the feature
    """

    # Fit a line to the end points
    x0, x1 = wave[0], wave[-1]
    y0, y1 = flux[0], flux[-1]
    m = (y0 - y1) / (x0 - x1)
    b = - m * x0 + y0

    continuum = m * wave + b
    norm_flux = flux / continuum
    pew = (x1 - x0) - np.trapz(y=norm_flux, x=wave)
    return norm_flux, pew


def feature_velocity(rest_frame, wave, flux, eflux=None, unit=None):
    """Calculate the velocity of a feature

    Args:
        rest_frame     (float): The rest frame wavelength of the feature
        wave         (ndarray): A sorted array of wavelengths for the feature
        flux (ndarray, uarray): An array of flux values for each wavelength
        eflux        (ndarray): An array of error values for each flux value
        unit      (PrefixUnit): Astropy unit for returned velocity (default km/s)

    Returns:
        The velocity of the feature
    """

    unit = units.km / units.s if unit is None else unit
    eflux = std_devs(flux) if eflux is None else eflux
    flux = nominal_values(flux)

    # Fit feature with a gaussian
    gaussian = lambda x, amplitude, avg, stddev, offset: \
        amplitude * np.exp(-((x - avg) ** 2) / (2 * stddev ** 2)) + offset

    (amplitude, avg, stddev, offset), cov = curve_fit(
        f=gaussian,
        xdata=wave,
        ydata=-(flux - 1),
        p0=[0.5, np.median(wave), 50., 0.],
        sigma=eflux)

    c = speed_of_light.to(unit).value
    avg = ufloat(avg, np.sqrt(cov[1][1]))
    return c * (
            ((((rest_frame - avg) / rest_frame) + 1) ** 2 - 1) /
            ((((rest_frame - avg) / rest_frame) + 1) ** 2 + 1)
    )


def calc_feature_properties(
        feat_name, wave, flux, feat_start, feat_end, eflux=None):
    """Calculate the properties of a single feature in a spectrum

    Velocity values are returned in km / s. Error values are determined
    both formally (summed in quadrature) and by re-sampling the feature
    boundaries five steps in either direction.

    Args:
        wave     (ndarray): An array of wavelengths
        flux     (ndarray): An array of flux for each wavelength
        eflux    (ndarray): The error for each flux value
        feat_name    (str): The name of the feature
        feat_start (float): Starting wavelength of the feature
        feat_end   (float): Ending wavelength of the feature

    Returns:
        - (The line velocity, its formal error, and its sampling error)
        - (The equivalent width, its formal error, and its sampling error)
        - (The feature area, its formal error, and its sampling error)
    """

    if eflux is not None:
        flux = uarray(flux, eflux)

    # Get rest frame location of the specified feature
    rest_frame = line_locations[feat_name]

    # Get indices for beginning and end of the feature
    idx_start = np.where(wave == feat_start)[0][0]
    idx_end = np.where(wave == feat_end)[0][0]
    if idx_end - idx_start <= 10:
        raise ValueError('Range too small. Please select a wider range')

    # We vary the beginning and end of the feature to estimate the error
    velocity, pequiv_width, area = [], [], []
    for i in np.arange(-5, 6):
        for j in np.arange(-5, 6):
            # Get sub-sampled wavelength/flux
            sample_start_idx = idx_start + i
            sample_end_idx = idx_end + j

            nw = wave[sample_start_idx: sample_end_idx]
            nf = flux[sample_start_idx: sample_end_idx]

            # Determine feature properties
            area.append(feature_area(nw, nf))
            norm_flux, pew = feature_pew(nw, nf)
            pequiv_width.append(pew)
            velocity.append(feature_velocity(rest_frame, nw, norm_flux))

    avg_velocity = np.mean(velocity)
    avg_ew = np.mean(pequiv_width)
    avg_area = np.mean(area)

    return (
        (
            avg_velocity.nominal_value,
            avg_velocity.std_dev,
            np.std(nominal_values(avg_velocity))
        ),
        (
            avg_ew.nominal_value,
            avg_ew.std_dev,
            np.std(nominal_values(pequiv_width))
        ),
        (
            avg_area.nominal_value,
            avg_area.std_dev,
            np.std(nominal_values(area))
        )
    )
