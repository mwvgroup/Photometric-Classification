#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module calculates the properties of spectral features such as area,
velocity, and  equivalent width.
"""

import numpy as np
from astropy import units
from astropy.constants import c as speed_of_light
from scipy.optimize import curve_fit

# Todo: Compare outputs to original IDL code
# Todo: Add flux error propagation

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
        wave (ndarray): A sorted array of wavelengths for the feature
        flux (ndarray): An array of flux values for each wavelength

    Returns:
        The area of the feature
    """

    x0 = wave[0]
    x1 = wave[-1]
    y0 = min(flux[0], flux[-1])
    y1 = max(flux[0], flux[-1])

    # Feature area = area under continuum - area under spectrum
    continuum_area = (x1 - x0) * y0 + (x1 - x0) * (y1 - y0) / 2.
    spectrum_area = np.trapz(y=flux, x=wave)
    return continuum_area - spectrum_area


def feature_pew(wave, flux):
    """Calculate the pseudo equivalent-width of a feature

    Args:
        wave (ndarray): A sorted array of wavelengths for the feature
        flux (ndarray): An array of flux values for each wavelength

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


def feature_velocity(wave, flux, eflux, rest_frame, unit=None):
    """Calculate the velocity of a feature

    Args:
        wave     (ndarray): A sorted array of wavelengths for the feature
        flux     (ndarray): An array of flux values for each wavelength
        eflux    (ndarray): An array of error values for each flux value
        rest_frame (float): The rest frame wavelength of the feature
        unit  (PrefixUnit): Astropy unit to return velocity in (default km / s)

    Returns:
        The velocity of the feature
    """

    unit = units.km / units.s if unit is None else unit
    gaussian = lambda x, amplitude, avg, stddev, offset: \
        amplitude * np.exp(-((x - avg) ** 2) / (2 * stddev ** 2)) + offset

    # Fit feature with a gaussian
    start = [0.5, np.median(wave), 50., 0.]
    (amplitude, avg, stddev, offset), cov = curve_fit(
        gaussian, wave, -(flux - 1), p0=start, sigma=eflux)

    c = speed_of_light.to(unit).value
    return c * (
            ((((rest_frame - avg) / rest_frame) + 1) ** 2 - 1) /
            ((((rest_frame - avg) / rest_frame) + 1) ** 2 + 1)
    )


def calc_feature_properties(
        wave, flux, eflux, feat_name, feat_start, feat_end):
    """Calculate the properties of a single feature in a spectrum

    Args:
        wave     (ndarray): An array of wavelengths
        flux     (ndarray): An array of flux for each wavelength
        eflux    (ndarray): The error for each flux value
        feat_name    (str): The name of the feature
        feat_start (float): Starting wavelength of the feature
        feat_end   (float): Ending wavelength of the feature

    Returns:
        - The average velocity
        - The standard deviation in the velocity
        - The average equivalent width
        - The standard deviation in the equivalent width
        - The average feature area
        - The standard deviation in feature area
    """

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
            nef = eflux[sample_start_idx: sample_end_idx]

            # Determine feature properties
            area.append(feature_area(nw, nf))
            norm_flux, pew = feature_pew(nw, nf)
            pequiv_width.append(pew)
            velocity.append(feature_velocity(nw, norm_flux, nef, rest_frame))

    return (
        np.mean(velocity), np.std(velocity),
        np.mean(pequiv_width), np.std(pequiv_width),
        np.mean(area), np.std(area)
    )
