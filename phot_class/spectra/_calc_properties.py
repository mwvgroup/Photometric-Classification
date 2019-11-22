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
dust_dir = _file_dir.parent / 'schlegel98_dust_map'
line_locations_path = _file_dir / 'features.yml'

dust_map = sfdmap.SFDMap(dust_dir)
with open(line_locations_path) as infile:
    line_locations = yaml.load(infile, Loader=yaml.FullLoader)

plt.ion()


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
        - The value of the continuum for each wavelength
        - The normalized flux
        - The pseudo equivalent-width of the feature
    """

    # Fit a line to the end points
    x0, x1 = wave[0], wave[-1]
    y0, y1 = flux[0], flux[-1]
    m = (y0 - y1) / (x0 - x1)
    b = - m * x0 + y0

    continuum = m * wave + b
    norm_flux = flux / continuum
    pew = (x1 - x0) - np.trapz(y=norm_flux, x=wave)
    return continuum, norm_flux, pew


def feature_velocity(rest_frame, wave, flux, unit=None):
    """Calculate the velocity of a feature

    Args:
        rest_frame     (float): The rest frame wavelength of the feature
        wave         (ndarray): A sorted array of wavelengths for the feature
        flux (ndarray, uarray): An array of flux values for each wavelength
        unit      (PrefixUnit): Astropy unit for returned velocity (default km/s)

    Returns:
        - The velocity of the feature
        - The average of the Gaussian
        - The Gaussian evaluated for each wavelength
    """

    eflux = std_devs(flux)
    flux = nominal_values(flux)
    unit = units.km / units.s if unit is None else unit

    # Fit feature with a gaussian
    def gaussian(x, _depth, _avg, _std, _offset):
        return -_depth * np.exp(-((x - _avg) ** 2) / (2 * _std ** 2)) + _offset

    (depth, avg, stddev, offset), cov = curve_fit(
        f=gaussian,
        xdata=wave,
        ydata=flux,
        p0=[0.5, np.median(wave), 50., 0],
        sigma=eflux if any(eflux) else None)

    fit = gaussian(wave, depth, avg, stddev, offset)
    if any(eflux):
        avg = ufloat(avg, np.sqrt(cov[1][1]))

    speed_of_light = c.to(unit).value
    vel = speed_of_light * (
            ((((rest_frame - avg) / rest_frame) + 1) ** 2 - 1) /
            ((((rest_frame - avg) / rest_frame) + 1) ** 2 + 1)
    )

    return vel, avg, fit


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
        - The starting wavelength of the feature
        - The ending wavelength of the feature
    """

    feat_start = find_peak_wavelength(
        wave, flux, feature['lower_blue'], feature['upper_blue'], 'min')

    feat_end = find_peak_wavelength(
        wave, flux, feature['lower_red'], feature['upper_red'], 'max')

    return feat_start, feat_end


def sample_feature_properties(
        feat_name, feat_start, feat_end, wave, flux,
        nstep=5, plot=False, debug=False):
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
        nstep            (int): The number of steps to take in each direction
        plot            (bool): Plot live fit results
        debug           (bool): Return samples instead of the average values

    Returns:
        - (The line velocity, its formal error, and its sampling error)
        - (The equivalent width, its formal error, and its sampling error)
        - (The feature area, its formal error, and its sampling error)
    """

    if plot:
        plt.clf()

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
        for j in np.arange(nstep, -nstep - 1, -1):
            # Get sub-sampled wavelength/flux
            sample_start_idx = idx_start + i
            sample_end_idx = idx_end + j

            nw = wave[sample_start_idx: sample_end_idx]
            nf = flux[sample_start_idx: sample_end_idx]

            # Determine feature properties
            area.append(feature_area(nw, nf))
            continuum, norm_flux, pew = feature_pew(nw, nf)
            pequiv_width.append(pew)

            # Todo: Add back in the velocity calculation
            vel, avg, fit = 0, 0, 0  # feature_velocity(rest_frame, nw, norm_flux)
            velocity.append(vel)

            if plot and i == -nstep and j == nstep:
                plt.plot(nw, nf, color='k', zorder=1)
                plt.xlabel('Wavelength')
                plt.ylabel('Flux')

            if plot:
                feat_id = line_locations[feat_name]['feature_id']
                avg_pew = np.average(pequiv_width)
                std_pew = np.std(pequiv_width)
                plt.title(feat_id + f' (pEW = {avg_pew:.2f} +\- {std_pew:.2f})')
                plt.fill_between(nw, nf, continuum, color='grey', alpha=.2, zorder=0)
                plt.axvline(nw[0], color='grey', linestyle='--', alpha=.25, zorder=2)
                plt.axvline(nw[-1], color='grey', linestyle='--', alpha=.25, zorder=2)
                plt.plot(nw, continuum, color='C0', linestyle='--', alpha=.4, zorder=3)

                # Todo: Add back in the velocity calculation
                # plt.plot(nw, fit * continuum, label='Fit', color='C2', alpha=.25, zorder=4)
                # plt.axvline(avg, color='C1', linestyle=':', zorder=5)

                plt.draw()
                plt.pause(.001)

    if plot:
        plt.pause(.5)
        plt.clf()

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


def _spectrum_properties(wave, flux, nstep=5, plot=False):
    """Calculate the properties of multiple features in a spectrum

    Velocity, pseudo equivalent width, and area are returned for
    each feature in ``line_locations`` along with their respective errors.

    Args:
        wave  (ndarray): An array of wavelengths in angstroms
        flux  (ndarray): An array of flux for each wavelength
        nstep     (int): The number of sampling steps to take
        plot     (bool): Plot live fit results

    Returns:
        A list of measurements and errors for each feature
    """

    # Iterate over features
    out_data = []
    for feat_name, feat_definition in line_locations.items():
        try:
            feat_start, feat_end = find_feature_bounds(
                wave, flux, feat_definition)

            samp_results = sample_feature_properties(
                feat_name, feat_start, feat_end, wave, flux,
                nstep=nstep, plot=plot
            )

            feat_properties = \
                [feat_name] + np.array(samp_results).flatten().tolist() + ['']

        except KeyboardInterrupt:
            raise

        except Exception as msg:
            feat_properties = [feat_name] + np.full(9, np.nan).tolist() + [
                str(msg)]

        out_data.append(feat_properties)

    return out_data


def _correct_spectrum(wave, flux, ra, dec, z, rv=3.1):
    """Rest frame spectra and correct for MW extinction

    Spectra are rest-framed and corrected for MW extinction using the
    Schlegel et al. 98 dust map and the Fitzpatrick et al. 99 extinction law.
    if rv is not given, a value of 1.7 is used for E(B - V) > .3 and a value
    of 3.1 is used otherwise.

    Args:
        wave (ndarray): Array of wavelength values
        flux (ndarray): Array of flux values
        ra     (float): Ra coordinate of the object
        dec    (float): Dec coordinate of the object
        z      (float): Redshift of the object
        rv     (float): Rv value to use for extinction

    Returns:
        - The rest framed wavelengths
        - The flux corrected for extinction
    """

    mwebv = dust_map.ebv(ra, dec, frame='fk5j2000', unit='degree')
    mag_ext = extinction.fitzpatrick99(wave, rv * mwebv, rv)
    flux = flux * 10 ** (0.4 * mag_ext)
    rest_wave = wave / (1 + z)
    return rest_wave, flux


def tabulate_spectral_properties(data_iter, nstep=5, rv=3.1, plot=False):
    """Tabulate spectral properties for multiple spectra of the same object

    Spectra are rest-framed and corrected for MW extinction using the
    Schlegel et al. 98 dust map and the Fitzpatrick et al. 99 extinction law.

    Args:
        data_iter (iterable[Table]): Iterable of spectroscopic data tables
        nstep                 (int): The number of sampling steps to take
        rv                  (float): Rv value to use for extinction
        plot                 (bool): Plot live fit results

    Returns:
        A Table with measurements for each spectrum and feature
    """

    table_rows = []
    for spectrum in data_iter:
        obj_id = spectrum.meta['obj_id']
        wave = spectrum['wavelength']
        flux = spectrum['flux']
        z = spectrum.meta['z']
        ra = spectrum.meta['ra']
        dec = spectrum.meta['dec']
        sid = spectrum.meta.get('spec_id', '?')
        date = spectrum['date'][0]

        try:
            type = spectrum['type'][0]

        except IndexError:
            type = '?'

        rest_wave, corrected_flux = _correct_spectrum(
            wave, flux, ra, dec, z, rv=rv)

        # Tabulate properties and add object Id to each measurement
        spec_properties = _spectrum_properties(
            rest_wave, corrected_flux, nstep=nstep, plot=plot)

        table_rows += [[obj_id, sid, date, type] + r for r in spec_properties]

    if not table_rows:
        table_rows = None

    # Format results as a table
    col_names = ['obj_id', 'sid', 'date', 'type', 'feat_name']
    dtype = ['U100', 'U100', 'U100', 'U100', 'U20']
    for value in ('vel', 'pew', 'area'):
        col_names.append(value)
        col_names.append(value + '_err')
        col_names.append(value + '_samperr')
        dtype += [float, float, float]

    col_names.append('msg')
    dtype.append('U1000')

    return Table(rows=table_rows, names=col_names, dtype=dtype)
