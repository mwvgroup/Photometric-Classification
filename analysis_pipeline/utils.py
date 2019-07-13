#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""A collection of general utilities used across the analysis pipeline."""

import numpy as np
from astropy.table import Column, unique
from astropy.table import Table

from . import _paths


def _split_bands(bands, lambda_eff):
    """Split band-passes into collections of blue and red bands

    Blue bands have an effective wavelength < 5500 Ang. Red bands have an
    effective wavelength >= 5500 Ang.

    Args:
        bands        (array[str]): Name of band-passes
        lambda_eff (array[float]): Effective wavelength of band-passes
    """

    is_blue = np.array(lambda_eff) < 5500
    b_array = np.array(bands)
    return b_array[is_blue], b_array[~is_blue]


def split_data(data_table, band_names, lambda_eff):
    """Split a data table into blue and red data (by rest frame)

    Split data by keeping filters that are red-ward or blue-ward of 5500 Ang.

    Args:
        data_table (Table): An SNCosmo input table with column 'band'
        band_names  (iter): List of all bands available in the survey
        lambda_eff  (iter): The effective wavelength of each band in band_names

    Returns:
        A SNCosmo input table with only blue bands
        A SNCosmo input table with only red bands
    """

    # Type cast to allow numpy indexing
    band_names = np.array(band_names)
    lambda_eff = np.array(lambda_eff)

    @np.vectorize
    def lambda_for_band(band):
        return lambda_eff[band_names == band]

    # Rest frame effective wavelengths for each observation
    z = data_table.meta['redshift']
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
    is_ok_diff = delta_lambda[min_indx, np.arange(delta_lambda.shape[1])] < 700

    # Split into blue and red band passes
    out_list = []
    for bands in _split_bands(band_names, lambda_eff):
        is_in_bands = np.isin(rest_frame_filters, bands)
        indices = np.logical_and(is_in_bands, is_ok_diff)
        out_list.append(data_table[indices])

    return out_list


def get_fit_results(survey, model, num_params, bands):
    """Get light-curve fits for a given survey, model, and number of parameters

    Args:
        survey  (module): An SNData submodule for a particular data release
        model    (Model): The SNCosmo model used to fit light-curves
        num_params (int): The number of params that were fit

    Returns:
        A DataFrame of fits in all bands or None
        A DataFrame of fits in blue bands or None
        A DataFrame of fits in red bands or None
    """

    all_path, blue_path, red_path = \
        _paths.get_fit_result_paths(model, survey, num_params)

    path_index = ['all', 'blue', 'red'].index(bands)
    path = [all_path, blue_path, red_path][path_index]

    data = Table.read(path)
    data = data.to_pandas()
    data.set_index('obj_id', inplace=True)

    return data


def get_priors(survey, model):
    """Get light-curve priors for a given survey and model

    Args:
        survey  (module): An SNData submodule for a particular data release
        model    (Model): The SNCosmo model used to fit light-curves

    Returns:
        An astropy table
    """

    path = _paths.get_priors_path(model, survey)
    data = Table.read(path).to_pandas()
    data.set_index('obj_id', inplace=True)
    return data


def save_manual_priors(obj_id, survey, model, priors_dict, message='-'):
    """Save priors to the analysis pipeline's internal file structure

    Args:
        obj_id       (str): The ID of the object to save priors for
        survey    (module): An sndata data access module
        model      (Model): The SNCosmo model of the priors
        priors_dict (dict): Dictionary of prior values
        message (str): Message to include with new priors (Default: '-')
    """

    auto_priors_path = _paths.get_priors_path(model, survey, manual=True)
    manual_priors_path = _paths.get_priors_path(model, survey, manual=True)

    try:
        existing_data = Table.read(manual_priors_path)
        existing_data['message'] = Column(existing_data['message'],
                                          dtype='U100')

    except FileNotFoundError:
        auto_data = Table.read(auto_priors_path)
        names = auto_data.colnames
        dtype = [float for _ in names]
        dtype[0] = dtype[-1] = 'U1000'
        existing_data = Table(names=names, dtype=dtype)

    priors_dict['obj_id'] = obj_id
    priors_dict['message'] = message
    existing_data.add_row(priors_dict)
    new_data = unique(existing_data, keep='last', keys=['obj_id'])
    new_data.write(manual_priors_path, overwrite=True)
