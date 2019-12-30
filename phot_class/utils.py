#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""A collection of general utilities used across the parent package."""

import signal
from copy import deepcopy

import numpy as np
import sncosmo
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


def parse_config_dict(obj_id, config_dict):
    """Return the priors and kwargs for a given object from a config file

    Args:
        obj_id       (str): The object id in the dictionary
        config_dict (dict): A dictionary with data from a config file

    Returns:
        - A dictionary with object priors for the hsiao_x1 model
        - A dictionary of fitting kwargs for the hsiao_x1 model
        - A dictionary with object priors for the sn91bg model
        - A dictionary of fitting kwargs for the sn91bg model
    """

    out_data = []
    for model in ('hsiao_x1', 'sn91bg'):
        for dtype in ('priors', 'kwargs'):
            object_data = config_dict[model].get(obj_id, {}).get(dtype, {})
            out_data.append(object_data)

    return tuple(out_data)


def calc_model_chisq(data, result, model):
    """Calculate the chi-squared for a given data table and model

    Chi-squareds are calculated using parameter values from ``model``. Degrees
    of freedom are calculated using the number of varied parameters specified
    is the ``result`` object.

    Args:
        data    (Table): An sncosmo input table
        model   (Model): An sncosmo Model
        result (Result): sncosmo fitting result

    Returns:
        The un-normalized chi-squared
        The number of data points used in the calculation
    """

    data = deepcopy(data)

    # Drop any data that is not withing the model's range
    min_band_wave = [sncosmo.get_bandpass(b).minwave() for b in data['band']]
    max_band_wave = [sncosmo.get_bandpass(b).maxwave() for b in data['band']]
    data = data[
        (data['time'] >= model.mintime()) &
        (data['time'] <= model.maxtime()) &
        (min_band_wave >= model.minwave()) &
        (max_band_wave <= model.maxwave())
        ]

    if len(data) == 0:
        raise ValueError('No data within model range')

    return sncosmo.chisq(data, model), len(data) - len(result.vparam_names)


def split_bands(bands, lambda_eff, redshift=0):
    """Split band-passes into collections of blue and red bands

    Blue bands have an rest frame effective wavelength < 5500 Ang. Red bands
    have a rest frame effective wavelength >= 5500 Ang.

    Args:
        bands        (array[str]): Name of band-passes
        lambda_eff (array[float]): Effective wavelength of band-passes
        redshift          (float): The redshift of the rest frame

    Returns:
        An array of blue filter names
        An array of red filter names
    """

    # Blueshift wavelengths to rest frame
    lambda_eff = np.array(lambda_eff) / (1 + redshift)

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
        cutoff     (float): The cutoff distance for dropping an observation

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
    within_dif_range = delta_lambda[
                           min_indx, np.arange(delta_lambda.shape[1])] < cutoff

    # Split into blue and red band passes
    out_list = []
    for bands in split_bands(band_names, lambda_eff):
        is_in_bands = np.isin(rest_frame_filters, bands)
        indices = np.logical_and(is_in_bands, within_dif_range)
        out_list.append(data_table[indices])

    return out_list


def classification_filter_factory(classifications, ftype='exclude'):
    """Returns function to determine whether data should be skipped/kept in an
     iterator based on its classification

    The function returned by this factory has signature
    ``returned_function(table: astropy.Table) -> boolean``. The boolean
    indicates whether the data should kept (i.e. not skipped). The class of
    each object is determined from ``table.meta['classification']``. If
    there is no '`classification'` key in the meta data, the return is True.

    Args:
         classifications (list[str]): A list of classifications to allow
         ftype                 (str): 'exclude' or 'include' the given classes

    Returns:
        A filter function for sndata
    """

    def filter_func(table):
        if 'classification' not in table.meta:
            return True

        if ftype == 'exclude':
            return table.meta['classification'] not in classifications

        elif ftype == 'include':
            return table.meta['classification'] in classifications

        else:
            raise ValueError(f'Unknown filter type: {ftype}')

    return filter_func
