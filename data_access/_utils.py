#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module provides utilities for selecting a subset of photometric
observations by rest frame wavelength.
"""

import numpy as np
from astropy.table import Table
import sncosmo


def keep_restframe_bands(data_table, bands, band_names, effective_lambda):
    """Return select rest-frame bandpasses from a table of photometric data

    Args:
        data_table      (Table): An SNCosmo input table with column 'band'
        bands            (iter): List of rest-frame band passes to keep
        band_names       (iter): List of all bands available in the survey
        effective_lambda (iter): The effective wavelength of each band in band_names

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


def parse_snoopy_data(path):
    """Return data from a snoopy file as an astropy table

    Args:
        path (str): The file path of a snoopy input file

    Returns:
        An astropy table with columns 'time', 'band', 'mag', and 'mag_err'
    """

    out_table = Table(names=['time', 'band', 'mag', 'mag_err'])
    with open(path) as ofile:
        # Get meta data from first line
        name, z, ra, dec = ofile.readline().split()
        out_table.meta['name'] = name
        out_table.meta['redshift'] = float(z)
        out_table.meta['ra'] = float(ra)
        out_table.meta['dec'] = float(dec)

        # Read photometric data from the rest of the file
        band = None
        for line in ofile.readlines():
            line_list = line.split()
            if line.startswith('filter'):
                band = line_list[1]
                continue

            time, mag, mag_err = line_list
            out_table.add_row([time, band, mag, mag_err])

    return out_table


@np.vectorize
def register_filter(file_path, filt_name):
    """Registers filter profiles with sncosmo if not already registered

    Args:
        file_path (str): Path of an ascii table with wavelength (Angstrom)
                          and transmission columns
        filt_name (str): The name of the registered filter.
    """

    available_bands = set(k[0] for k in sncosmo.bandpasses._BANDPASSES._loaders.keys())
    if filt_name not in available_bands:
        print('Registering band "{}" with sncosmo.'.format(filt_name))
        filt_data = np.genfromtxt(file_path).T
        band = sncosmo.Bandpass(filt_data[0], filt_data[1])
        sncosmo.register(band, filt_name, force=True)
