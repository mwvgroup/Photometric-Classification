#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module provides utilities for selecting a subset of photometric
observations by rest frame wavelength.
"""

import numpy as np


@np.vectorize
def band_index_mapping(band_names, i):
    """Vectorized mapping between index and filter name sdss<"ugriz"[i]>

    Args:
        i (int or str): Index of bandpass or bandpass name

    Returns:
        If the argument is an integer, return the filter name sdss<"ugriz"[i]>
        If the argument is a filter name, return the index
    """

    if isinstance(i, int):  # If i is index return, band name
        return band_names[i]

    elif isinstance(i, str):  # If i is band name, return index
        return band_names.index(i)

    else:
        raise ValueError('Argument must be int or str')


def keep_restframe_bands(data_table, bands):
    """Return select rest-frame bandpasses from an SNCsomo input table

    Args:
        data_table (Table): An SNCosmo input table as given
                              by iter_sncosmo_input
        bands       (list): List of rest-frame bandpasses to keep

    Returns:
        A new input table for SNCosmo only containing select rest frame bands
    """

    # Effective wavelengths for SDSS filters ugriz in angstroms
    # https://www.sdss.org/instruments/camera/#Filters
    effective_lambda = np.array([3551, 4686, 6166, 7480, 8932])

    # Rest frame effective wavelengths for each observation
    z = data_table.meta['redshift']
    filter_indices = band_index_mapping(data_table['band'])
    rest_frame_lambda = effective_lambda[filter_indices] / (1 + z)

    # Get the name of the observer frame band with the smallest distance
    # to each rest frame lambda
    delta_lambda = np.array([np.abs(rest_frame_lambda - l_eff) for l_eff in effective_lambda])
    min_indx = np.argmin(delta_lambda, axis=0)
    rest_frame_filters = band_index_mapping(min_indx)

    # Keep only the specified filters that are within 1000 Angstroms of the
    # rest frame effective wavelength
    is_ok_diff = delta_lambda[min_indx, np.arange(delta_lambda.shape[1])] < 1000
    is_in_bands = np.isin(rest_frame_filters, bands)
    indices = np.logical_and(is_in_bands, is_ok_diff)

    return data_table[indices]
