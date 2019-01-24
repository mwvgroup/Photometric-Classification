#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module provides utilities for selecting a subset of photometric
observations by rest frame wavelength.
"""

import numpy as np


def keep_restframe_bands(data_table, bands, band_names, effective_lambda):
    """Return select rest-frame bandpasses from an SNCsomo input table

    Args:
        data_table (Table): An SNCosmo input table as given
                              by iter_sncosmo_input
        bands       (list): List of rest-frame bandpasses to keep

    Returns:
        A new input table for SNCosmo only containing select rest frame bands
    """

    band_names = np.array(band_names)

    # Rest frame effective wavelengths for each observation
    z = data_table.meta['redshift']
    lambda_per_band = [effective_lambda[band_names == band][0] for band in data_table['band']]
    rest_frame_lambda = lambda_per_band / (1 + z)

    # Get the name of the observer frame band with the smallest distance
    # to each rest frame lambda
    delta_lambda = np.array([np.abs(rest_frame_lambda - l_eff) for l_eff in effective_lambda])
    min_indx = np.argmin(delta_lambda, axis=0)
    rest_frame_filters = np.array(band_names)[min_indx]

    # Keep only the specified filters that are within 1000 Angstroms of the
    # rest frame effective wavelength
    is_ok_diff = delta_lambda[min_indx, np.arange(delta_lambda.shape[1])] < 1000
    is_in_bands = np.isin(rest_frame_filters, bands)
    indices = np.logical_and(is_in_bands, is_ok_diff)

    return data_table[indices]
