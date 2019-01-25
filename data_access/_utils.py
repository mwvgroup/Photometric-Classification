#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module provides utilities for selecting a subset of photometric
observations by rest frame wavelength.
"""

import numpy as np


def keep_restframe_bands(data_table, bands, band_names, effective_lambda):
    """Return select rest-frame bandpasses from a table of photometric data

    Args:
        data_table      (Table): An SNCosmo input table with column 'band'
        bands            (list): List of rest-frame bandpasses to keep
        band_names       (list): List of all bands available in the survey
        effective_lambda (list): The effective wavelength of each band in band_names

    Returns:
        A new input table for SNCosmo only containing select rest frame bands
    """

    # Type cast to allow numpy indexing
    band_names = np.array(band_names)
    effective_lambda = np.array(effective_lambda)

    @np.vectorize
    def lambda_for_band(band):
        return effective_lambda[band_names == band]

    # Rest frame effective wavelengths for each observation
    z = data_table.meta['redshift']
    observed_lambda = lambda_for_band(data_table['band'])
    rest_frame_lambda = observed_lambda / (1 + z)

    # Get the name of the observer frame band with the smallest distance
    # to each rest frame lambda
    delta_lambda = np.array([
        np.abs(rest_frame_lambda - l_eff) for l_eff in effective_lambda
    ])

    min_indx = np.argmin(delta_lambda, axis=0)
    rest_frame_filters = np.array(band_names)[min_indx]

    # Keep only the specified filters that are within 700 Angstroms of the
    # rest frame effective wavelength
    is_ok_diff = delta_lambda[min_indx, np.arange(delta_lambda.shape[1])] < 700
    is_in_bands = np.isin(rest_frame_filters, bands)
    indices = np.logical_and(is_in_bands, is_ok_diff)

    return data_table[indices]