# !/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module provides a graphical interface for calculating
the properties of spectral features.
"""

import numpy as np
from matplotlib import pyplot as plt
from uncertainties import nominal_value, std_dev
from uncertainties.unumpy import nominal_values

from .calc_properties import (
    bin_spectrum,
    correct_extinction,
    feature_area,
    feature_pew,
    guess_feature_bounds,
    line_locations)


def _draw_measurement(
        wave, flux, continuum, feat_name, eq_width, pause=.001):
    """Shade in the EW, continuum, and position of a spectral feature

    Args:
        wave      (ndarray): An array of wavelengths in angstroms
        flux      (ndarray): An array of flux for each wavelength
        continuum (ndarray): The continuum flux
        feat_name     (str): The name of the feature
        eq_width  (ndarray): Array of equivalent width measurements
        pause       (float): How long to pause after drawing
    """

    feat_id = line_locations[feat_name]['feature_id']
    avg_pew = np.average(eq_width)
    std_pew = np.std(eq_width)

    plt.title(feat_id + rf' (pEW = {avg_pew:.2f} $\pm$ {std_pew:.2f})')
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')

    plt.fill_between(wave, flux, continuum, color='grey', alpha=.2,
                     zorder=0)
    plt.axvline(wave[0], color='grey', linestyle='--', alpha=.25, zorder=2)
    plt.axvline(wave[-1], color='grey', linestyle='--', alpha=.25,
                zorder=2)
    plt.plot(wave, continuum, color='C0', linestyle='--', alpha=.4,
             zorder=3)

    # Todo: Add back in the velocity calculation
    # plt.plot(nw, fit * continuum, label='Fit', color='C2', alpha=.25, zorder=4)
    # plt.axvline(avg, color='C1', linestyle=':', zorder=5)

    plt.draw()
    plt.pause(pause)


class SpectrumInspector:

    def __init__(self, spectrum):

        # Meta data about the spectrum
        self.obj_id = spectrum.meta['obj_id']
        self.wave = spectrum['wavelength']
        self.flux = spectrum['flux']
        self.z = spectrum.meta['z']
        self.ra = spectrum.meta['ra']
        self.dec = spectrum.meta['dec']
        self.sid = spectrum.meta.get('spec_id', '?')
        self.date = spectrum['date'][0]
        try:
            self.spec_type = spectrum['type'][0]

        except IndexError:
            self.spec_type = '?'

        # Place holders for intermediate analysis results
        self.bin_wave, self.bin_flux = None, None
        self.corrected_flux, self.rest_wave = None, None

    def prepare_spectrum(self, bin_size, method, rv=None):
        """Bin, correct for extinction, and rest-frame the spectrum

        Args:
            bin_size (float): Bin size in units of Angstroms
            method     (str): Either 'avg' or 'sum' the values of each bin
            rv       (float): Rv value to use for extinction (Default: 3.1)
        """

        self.bin_wave, self.bin_flux = bin_spectrum(
            self.wave, self.flux, bin_size=bin_size, method=method)

        self.rest_wave, self.corrected_flux = correct_extinction(
            self.bin_wave, self.bin_flux, self.ra, self.dec, self.z, rv=rv)

    def ask_feature_bounds(self, wave, flux, feat_definition):

        gstart, gend = guess_feature_bounds(wave, flux, feat_definition)

        plt.clf()
        plt.plot(self.bin_wave, self.bin_flux, color='k')

        vline_style = dict(color='grey', linestyle='--', alpha=.25, zorder=2)
        low_line = plt.axvline(gstart, **vline_style)
        upper_line = plt.axvline(gend, **vline_style)
        plt.xlim(feat_definition['lower_blue'] - 20,
                 feat_definition['upper_red'] + 20)

        xy_low = plt.ginput(1)
        lower_bound = wave[(np.abs(wave - xy_low[0][0])).argmin()]
        low_line.remove()
        plt.axvline(lower_bound, **vline_style)

        xy_high = plt.ginput(1)
        upper_bound = wave[(np.abs(wave - xy_high[0][0])).argmin()]
        upper_line.remove()
        plt.axvline(upper_bound, **vline_style)

        return lower_bound, upper_bound

    @staticmethod
    def _sample_feature_properties(
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
            nstep            (int): Number of samples taken in each direction
            debug           (bool): Return samples instead of averaged values

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
                # vel, avg, fit = feature_velocity(rest_frame, nw, norm_flux)
                vel, avg, fit = 0, 0, 0
                velocity.append(vel)

                if i == -nstep and j == nstep:
                    plt.xlim(min(nw), max(nw))

                _draw_measurement(nw, nf, continuum, feat_name, pequiv_width)

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

    def _spectrum_properties(self, wave, flux, nstep=5):
        """Calculate the properties of multiple features in a spectrum

        Velocity, pseudo equivalent width, and area are returned for
        each feature in ``line_locations`` along with their respective errors.

        Args:
            wave  (ndarray): An array of wavelengths in angstroms
            flux  (ndarray): An array of flux for each wavelength
            nstep     (int): The number of sampling steps to take

        Returns:
            A list of measurements and errors for each feature
        """

        # Iterate over features
        out_data = []
        for feat_name, feat_definition in line_locations.items():
            try:

                feat_start, feat_end = self.ask_feature_bounds(
                    wave, flux, feat_definition)

                samp_results = self._sample_feature_properties(
                    feat_name, feat_start, feat_end, wave, flux, nstep=nstep
                )

                samp_results = np.array(samp_results).flatten().tolist()
                feat_properties = [feat_name] + samp_results + ['']

            except KeyboardInterrupt:
                raise

            except Exception as msg:
                masked_values = np.full(9, np.nan).tolist()
                feat_properties = [feat_name] + masked_values + [str(msg)]

            out_data.append(feat_properties)

        return out_data

    def run(self, bin_size=5, method='avg', nstep=5, rv=None):

        self.prepare_spectrum(bin_size, method, rv)

        # Tabulate spectral properties
        spec_properties = self._spectrum_properties(
            self.rest_wave, self.corrected_flux, nstep=nstep)

        plt.close()

        # Add object Id and other meta data to each measurement
        meta_data = [self.obj_id, self.sid, self.date, self.spec_type]
        return [meta_data + r for r in spec_properties]
