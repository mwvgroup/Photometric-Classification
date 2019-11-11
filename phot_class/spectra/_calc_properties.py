#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module calculates the properties of spectral features such as area,
velocity, and  equivalent width. All functions in this module are built to
support ``uarray`` objects from the ``uncertainties`` package as inputs.
"""

from pathlib import Path

import numpy as np
import sfdmap
import yaml
from astropy import units
from astropy.constants import c as speed_of_light
from astropy.table import Table
from scipy.optimize import curve_fit
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values, std_devs, uarray

# File paths for external data
_file_dir = Path(__file__).resolve().parent
dust_path = _file_dir.parent / 'schlegel98_dust_map'
line_locations_path = _file_dir / 'features.yml'

DUST_MAP = sfdmap.SFDMap(dust_path)
with open(line_locations_path) as infile:
    line_locations = yaml.load(infile, Loader=yaml.FullLoader)


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

    Error values can be provided using the ``eflux`` argument, or
    by specifying ``flux`` and an array of ``ufloat`` objects.

    Args:
        rest_frame     (float): The rest frame wavelength of the feature
        wave         (ndarray): A sorted array of wavelengths for the feature
        flux (ndarray, uarray): An array of flux values for each wavelength
        eflux        (ndarray): An array of error values for each flux value
        unit      (PrefixUnit): Astropy unit for returned velocity (default km/s)

    Returns:
        The velocity of the feature
    """

    if eflux is None and not any(std_devs(flux)):
        raise ValueError('No error values provided.')

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


def _get_peak_wavelength(
        wavelength, flux, lower_bound, upper_bound, behavior='min'):
    """Return wavelength of the maximum flux within given wavelength bounds

    The behavior argument can be used to select the 'min' or 'max' wavelength
    when there are multiple wavelengths having the same peak flux value. The
    default behavior is 'min'.

    Args:
        wavelength (ndarray): An array of wavelength values
        flux       (ndarray): An array of flux values
        lower_bound  (float): Lower wavelength boundary
        upper_bound  (float): Upper wavelength boundary
        behavior       (str): Return the 'min' or 'max' wavelength

    Returns:
        The wavelength for the maximum flux value
        The maximum flux value
    """

    # Make sure the given spectrum spans the given wavelength bounds
    if (min(wavelength) > lower_bound) or (upper_bound > max(wavelength)):
        raise ValueError('Feature not in spectral wavelength range.')

    # Select the portion of the spectrum within the given bounds
    feature_indices = (lower_bound <= wavelength) & (wavelength <= upper_bound)
    feature_flux = flux[feature_indices]
    feature_wavelength = wavelength[feature_indices]

    peak_indices = np.argwhere(feature_flux == np.max(feature_flux))
    behavior_func = getattr(np, behavior)
    return behavior_func(feature_wavelength[peak_indices])


def get_feature_bounds(wavelength, flux, feature):
    """Get the start and end wavelengths / flux for a given feature

    Args:
        wavelength (ndarray): An array of wavelength values
        flux       (ndarray): An array of flux values
        feature        (row): A dictionary defining feature parameters

    Returns:
        The starting wavelength of the feature
        The ending wavelength of the feature
    """

    feat_start = _get_peak_wavelength(
        wavelength, flux, feature['lower_blue'], feature['upper_blue'], 'min')

    feat_end = _get_peak_wavelength(
        wavelength, flux, feature['lower_red'], feature['upper_red'], 'max')

    return feat_start, feat_end


def calc_feature_properties(
        feat_name, wave, flux, feat_start, feat_end, eflux=None):
    """Calculate the properties of a single feature in a spectrum

    Velocity values are returned in km / s. Error values are determined
    both formally (summed in quadrature) and by re-sampling the feature
    boundaries five steps in either direction.

    Args:
        wave     (ndarray): An array of wavelengths
        flux     (ndarray): An array of flux for each wavelength
        feat_name    (str): The name of the feature
        feat_start (float): Starting wavelength of the feature
        feat_end   (float): Ending wavelength of the feature
        eflux    (ndarray): The optional error for each flux value

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


def _spectrum_properties(wave, flux, z, ra, dec, eflux=None):
    """Calculate the properties of multiple features in a spectrum

    Velocity, pseudo equivalent width, and area are returned for
    each feature in ``line_locations`` along with their respective errors.

    Args:
        wave  (ndarray): An array of wavelengths
        flux  (ndarray): An array of flux for each wavelength
        eflux (ndarray): The optional error for each flux value

    Returns:
        A list of measurements and errors for each feature
    """

    # rest-frame spectra
    wave /= 1 + z

    # correct for extinction
    mwebv = DUST_MAP.ebv(ra, dec)
    # Todo

    # Iterate over features
    out_data = []
    for feat_name, feat_properties in line_locations:
        feat_start, feat_end = get_feature_bounds(wave, flux, feat_properties)
        feat_properties = calc_feature_properties(
            feat_name, wave, flux, feat_start, feat_end, eflux)

        out_data += np.array(feat_properties).flatten().tolist()

    return out_data


def tabulate_spectral_properties(wave, flux, z, ra, dec, eflux=None):
    """Tabulate spectral properties for multiple spectra

    Args:
        wave  (iterable[ndarray, uarray]): Wavelengths for each spectrum
        flux          (iterable[ndarray]): Flux for each spectrum
        z               (iterable[float]): Redshift of each spectrum
        ra              (iterable[float]): Right Ascension for each spectrum
        dec             (iterable[float]): Declination for each spectrum
        eflux         (iterable[ndarray]): The optional error for each flux

    Returns:
        A Table with measurements for each spectrum and feature
    """

    # Calculate feature properties
    if eflux is None:
        data_iter = zip(wave, flux, z, ra, dec)

    else:
        data_iter = zip(wave, flux, z, ra, dec, eflux)

    rows = [_spectrum_properties(*args) for args in data_iter]

    # Format results as a table
    col_names = ['date']
    for feat_name in line_locations:
        col_names.append(feat_name)
        col_names.append(feat_name + '_err')
        col_names.append(feat_name + '_samperr')

    dtype = [float for _ in col_names]
    return Table(rows=rows, names=col_names, dtype=dtype)
