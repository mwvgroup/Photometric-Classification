#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module calculates the properties of spectral features such as area,
velocity, and  equivalent width. All functions in this module are built to
support ``uarray`` objects from the ``uncertainties`` package as inputs.
"""

from pathlib import Path

import extinction
import numpy as np
import sfdmap
import yaml
from astropy import units
from astropy.constants import c
from astropy.table import Table
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from uncertainties import nominal_value, std_dev, ufloat
from uncertainties.unumpy import nominal_values, std_devs

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


def feature_velocity(rest_frame, wave, flux, unit=None, plot=False):
    """Calculate the velocity of a feature

    Args:
        rest_frame     (float): The rest frame wavelength of the feature
        wave         (ndarray): A sorted array of wavelengths for the feature
        flux (ndarray, uarray): An array of flux values for each wavelength
        unit      (PrefixUnit): Astropy unit for returned velocity (default km/s)

    Returns:
        The velocity of the feature
    """

    eflux = std_devs(flux)
    flux = nominal_values(flux)
    unit = units.km / units.s if unit is None else unit

    # Fit feature with a gaussian
    def gaussian(x, depth, avg, stddev, offset):
        return -depth * np.exp(-((x - avg) ** 2) / (2 * stddev ** 2)) + offset

    (depth, avg, stddev, offset), cov = curve_fit(
        f=gaussian,
        xdata=wave,
        ydata=flux,
        p0=[0.5, np.median(wave), 50., 0],
        sigma=eflux if any(eflux) else None)

    if plot:
        plt.scatter(wave, flux, label='Measured', color='C0', s=5)
        plt.errorbar(wave, flux, yerr=eflux, color='C0', linestyle='')
        fit = gaussian(wave, depth, avg, stddev, offset)
        plt.plot(wave, fit, label='Fit', color='C1')
        plt.axvline(avg, color='k', linestyle='--', label=f'Average: {avg}')
        plt.legend()
        plt.show()

    if any(eflux):
        avg = ufloat(avg, np.sqrt(cov[1][1]))

    speed_of_light = c.to(unit).value
    return speed_of_light * (
            ((((rest_frame - avg) / rest_frame) + 1) ** 2 - 1) /
            ((((rest_frame - avg) / rest_frame) + 1) ** 2 + 1)
    )


def find_peak_wavelength(
        wave, flux, lower_bound, upper_bound, behavior='min'):
    """Return wavelength of the maximum flux within given wavelength bounds

    The behavior argument can be used to select the 'min' or 'max' wavelength
    when there are multiple wavelengths having the same peak flux value. The
    default behavior is 'min'.

    Args:
        wave       (ndarray): An array of wavelength values
        flux       (ndarray): An array of flux values
        lower_bound  (float): Lower wavelength boundary
        upper_bound  (float): Upper wavelength boundary
        behavior       (str): Return the 'min' or 'max' wavelength

    Returns:
        The wavelength for the maximum flux value
        The maximum flux value
    """

    # Make sure the given spectrum spans the given wavelength bounds
    if (min(wave) > lower_bound) or (upper_bound > max(wave)):
        raise ValueError('Feature not in spectral wavelength range.')

    # Select the portion of the spectrum within the given bounds
    feature_indices = (lower_bound <= wave) & (wave <= upper_bound)
    feature_flux = flux[feature_indices]
    feature_wavelength = wave[feature_indices]

    peak_indices = np.argwhere(feature_flux == np.max(feature_flux))
    behavior_func = getattr(np, behavior)
    return behavior_func(feature_wavelength[peak_indices])


def find_feature_bounds(wave, flux, feature):
    """Get the start and end wavelengths / flux for a given feature

    Args:
        wave (ndarray): An array of wavelength values
        flux (ndarray): An array of flux values
        feature (dict): A dictionary defining feature parameters

    Returns:
        The starting wavelength of the feature
        The ending wavelength of the feature
    """

    feat_start = find_peak_wavelength(
        wave, flux, feature['lower_blue'], feature['upper_blue'], 'min')

    feat_end = find_peak_wavelength(
        wave, flux, feature['lower_red'], feature['upper_red'], 'max')

    return feat_start, feat_end


def sample_feature_properties(
        feat_name, feat_start, feat_end, wave, flux, nstep=5, debug=False):
    """Calculate the properties of a single feature in a spectrum

    Velocity values are returned in km / s. Error values are determined
    both formally (summed in quadrature) and by re-sampling the feature
    boundaries ``nstep`` flux measurements in either direction.

    Args:
        feat_name        (str): The name of the feature
        feat_start     (float): Starting wavelength of the feature
        feat_end       (float): Ending wavelength of the feature
        wave         (ndarray): An array of wavelengths
        flux (ndarray, uarray): An array of flux for each wavelength
        nstep          (float): The number of steps to take in each direction
        debug           (bool): Return samples instead of the average values

    Returns:
        - (The line velocity, its formal error, and its sampling error)
        - (The equivalent width, its formal error, and its sampling error)
        - (The feature area, its formal error, and its sampling error)
    """

    # Get rest frame location of the specified feature
    rest_frame = line_locations[feat_name]['restframe']

    # Get indices for beginning and end of the feature
    idx_start = np.where(wave == feat_start)[0][0]
    idx_end = np.where(wave == feat_end)[0][0]
    if idx_end - idx_start <= 10:
        raise ValueError('Range too small. Please select a wider range')

    # We vary the beginning and end of the feature to estimate the error
    velocity, pequiv_width, area = [], [], []
    for i in np.arange(-nstep, nstep + 1):
        for j in np.arange(-nstep, nstep + 1):
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

    if debug:
        return velocity, pequiv_width, area

    avg_velocity = np.mean(velocity)
    avg_ew = np.mean(pequiv_width)
    avg_area = np.mean(area)

    return (
        (
            nominal_value(avg_velocity),
            std_dev(avg_velocity),
            np.std(nominal_values(avg_velocity))
        ),
        (
            nominal_value(avg_ew),
            std_dev(avg_ew),
            np.std(nominal_values(pequiv_width))
        ),
        (
            nominal_value(avg_area),
            std_dev(avg_area),
            np.std(nominal_values(area))
        )
    )


def _spectrum_properties(wave, flux, z, ra, dec, rv=3.1):
    """Calculate the properties of multiple features in a spectrum

    Velocity, pseudo equivalent width, and area are returned for
    each feature in ``line_locations`` along with their respective errors.
    Spectra are rest-framed and corrected for MW extinction.

    Args:
        wave  (ndarray): An array of wavelengths in angstroms
        flux  (ndarray): An array of flux for each wavelength
        eflux (ndarray): The optional error for each flux value

    Returns:
        A list of measurements and errors for each feature
    """

    # rest-frame spectra
    wave = wave / (1 + z)

    # correct for extinction
    mwebv = DUST_MAP.ebv(ra, dec, frame='fk5j2000', unit='degree')
    mag_ext = extinction.fitzpatrick99(wave, rv * mwebv, rv)
    flux = flux * 10 ** (0.4 * mag_ext)

    # Iterate over features
    out_data = []
    for feat_name, feat_properties in line_locations.items():
        feat_start, feat_end = find_feature_bounds(wave, flux, feat_properties)
        feat_properties = sample_feature_properties(
            feat_name, feat_start, feat_end, wave, flux)

        out_data += np.array(feat_properties).flatten().tolist()

    return out_data


def tabulate_spectral_properties(date, wave, flux, z, ra, dec):
    """Tabulate spectral properties for multiple spectra

    Spectra are rest-framed and corrected for MW extinction using the
    Schlegel et al. 98 dust map and the Fitzpatrick et al. 99 extinction law.

    Args:
        date           (iter[float]): The date of observation for each spectrum
        wave (iter[ndarray, uarray]): Wavelengths for each spectrum in angstroms
        flux         (iter[ndarray]): Flux for each spectrum in arbitrary units
        z              (iter[float]): Redshift of each spectrum
        ra             (iter[float]): The J2000 Ra for each spectrum in degrees
        dec            (iter[float]): The J2000 Dec for each spectrum in degrees

    Returns:
        A Table with measurements for each spectrum and feature
    """

    # Calculate feature properties
    data_iter = zip(date, wave, flux, z, ra, dec)
    rows = [[date] + _spectrum_properties(*args[1:]) for args in data_iter]
    if not rows:
        rows = None

    # Format results as a table
    col_names = ['date']
    for feat_name in line_locations:
        for value in ('_vel', '_pew', '_area'):
            col_names.append(feat_name + value)
            col_names.append(feat_name + value + '_err')
            col_names.append(feat_name + value + '_samperr')

    dtype = [float for _ in col_names]
    return Table(rows, names=col_names, dtype=dtype)