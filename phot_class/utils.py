#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""A collection of general utilities used across the parent package."""

import signal
from copy import deepcopy

import numpy as np
from astropy.table import Table


class timeout:
    """A timeout context manager"""

    def __init__(self, seconds=1, error_message='Timeout'):
        """A timeout context manager

        Args:
            seconds       (int): The number of seconds until timeout
            error_message (str): The TimeOutError message on timeout
        """

        self.seconds = seconds
        self.error_message = error_message

    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)

    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)

    def __exit__(self, type_, value, traceback):
        signal.alarm(0)


def calc_model_chisq(data, model):
    """Calculate the chi-squared for a given data table and model

    Args:
        data  (Table): An SNCosmo input table
        model (Model): An SNCosmo Model

    Returns:
        The un-normalized chi-squared
        The number of data points used in the calculation
    """

    data = deepcopy(data)
    data.sort('time')  # To keep sncosmo happy

    while True:
        try:
            # Model flux and keep only non-zero values
            model_flux = model.bandflux(data['band'], data['time'])

        except ValueError as err:
            if 'outside spectral range' not in str(err):
                raise  # Something went wrong that we aren't expecting

            # sncosmo doesn't like bands that are outside of the model's range
            data = data[data['band'] != err.args[0].split()[1][1:-1]]
            if len(data) == 0:
                raise RuntimeError('Ran out of data inside model range!')

        else:
            data = data[model_flux > 0]
            model_flux = model_flux[model_flux > 0]

            residuals = model_flux - data['flux']
            chisq = np.sum((residuals / data['fluxerr']) ** 2)
            return chisq, len(data)


def split_bands(bands, lambda_eff):
    """Split band-passes into collections of blue and red bands

    Blue bands have an effective wavelength < 5500 Ang. Red bands have an
    effective wavelength >= 5500 Ang.

    Args:
        bands        (array[str]): Name of band-passes
        lambda_eff (array[float]): Effective wavelength of band-passes

    Returns:
        An array of blue filter names
        An array of red filter names
    """

    is_blue = np.array(lambda_eff) < 5500
    band_array = np.array(bands)
    return band_array[is_blue], band_array[~is_blue]


def split_data(data_table, band_names, lambda_eff, z, cutoff=700):
    """Split a data table into blue and red data (by rest frame)

    Wavelengths are expected to be in angstroms. Split data by keeping filters
    that are red-ward or blue-ward of 5500 Ang. If the closest rest frame
    filter for an observation is more than ``cutoff`` angstroms away, drop the
    observation.

    Args:
        data_table (Table): An SNCosmo input table with column 'band'
        band_names  (iter): List of all bands available in the survey
        lambda_eff  (iter): The effective wavelength of each band in band_names
        z          (float): The redshift of the observed target

    Returns:
        A SNCosmo input table with only blue bands
        A SNCosmo input table with only red bands
    """

    # Check an effective wavelength was specified for each band in the
    # data table. This avoids a cryptic error message later on.
    observed_bands = np.unique(data_table['band'])
    band_has_lambda_eff = np.isin(observed_bands, band_names)
    if not band_has_lambda_eff.all():
        missing_bands = observed_bands[~band_has_lambda_eff]
        raise ValueError(f'Missing effective wavelength for: {missing_bands}')

    # Type cast to allow numpy indexing
    band_names = np.array(band_names)
    lambda_eff = np.array(lambda_eff)

    @np.vectorize
    def lambda_for_band(band):
        return lambda_eff[band_names == band]

    # Calculate rest frame effective wavelengths for each observation
    observed_lambda = lambda_for_band(data_table['band'])
    rest_frame_lambda = observed_lambda / (1 + z)

    # Get the name of the observer frame band with the smallest distance
    # to each rest frame lambda
    delta_lambda = np.array([
        np.abs(rest_frame_lambda - l_eff) for l_eff in lambda_eff])

    min_indx = np.argmin(delta_lambda, axis=0)
    rest_frame_filters = np.array(band_names)[min_indx]

    # Keep only the specified filters that are within 700 Angstroms of the
    # rest frame effective wavelength
    is_ok_diff = delta_lambda[min_indx, np.arange(delta_lambda.shape[1])] < cutoff

    # Split into blue and red band passes
    out_list = []
    for bands in split_bands(band_names, lambda_eff):
        is_in_bands = np.isin(rest_frame_filters, bands)
        indices = np.logical_and(is_in_bands, is_ok_diff)
        out_list.append(data_table[indices])

    return out_list


def classification_filter_factory(classifications):
    """Factory function that returns an sndata filter function

    Filter function limits data iterator to targets of a given fitting

    Args:
         classifications (list[str]): A list of fitting to allow

    Returns:
        A filter function for sndata
    """

    def filter_func(table):
        return 'fitting' not in table.meta or \
               table.meta['fitting'] in classifications

    return filter_func
