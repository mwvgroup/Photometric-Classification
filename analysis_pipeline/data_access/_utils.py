#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module provides utilities used by various submodules."""

import os
import tarfile

import numpy as np
import requests
import sncosmo
from astropy.table import Table


def keep_restframe_bands(data_table, bands, band_names, effective_lambda):
    """Return select rest-frame band-passes from a table of photometric data

    Args:
        data_table      (Table): An SNCosmo input table with column 'band'
        bands            (iter): List of rest-frame band-passes to keep
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

    out_table = Table(
        names=['time', 'band', 'mag', 'mag_err'],
        dtype=[float, object, float, float]
    )

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


def register_filter(file_path, filt_name):
    """Registers filter profiles with sncosmo if not already registered

    Args:
        file_path (str): Path of an ascii table with wavelength (Angstrom)
                          and transmission columns
        filt_name (str): The name of the registered filter.
    """

    # Get set of registered builtin and custom bandpasses
    available_bands = set(k[0] for k in sncosmo.bandpasses._BANDPASSES._loaders)
    available_bands.update(k[0] for k in sncosmo.bandpasses._BANDPASSES._instances)

    # Register the new bandpass
    if filt_name not in available_bands:
        filt_data = np.genfromtxt(file_path).T
        band = sncosmo.Bandpass(filt_data[0], filt_data[1])
        sncosmo.register(band, filt_name, force=False)


def _download_file(url, out_path):
    """Download a specified file to a given output path

    Any top level .tar.gz archives will be automatically unzipped.

    Args:
        url      (str): URL of the file to download
        out_path (str): The path where the downloaded file should be written
    """

    # Make temporary file path
    if os.path.isdir(out_path):
        temp_path = os.path.join(out_path, '.temp')
        out_dir = out_path

    else:
        temp_path = out_path + '.temp'
        out_dir = os.path.dirname(out_path)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Download data
    response = requests.get(url)
    response.raise_for_status()
    with open(temp_path, 'wb') as ofile:
        ofile.write(response.content)

    # Unzip file if its an archive
    if url.endswith(".tar.gz") or url.endswith(".tgz"):
        with tarfile.open(temp_path, "r:gz") as data:
            data.extractall(out_dir)

        os.remove(temp_path)
        for (dirpath, dirnames, filenames) in os.walk(out_dir):
            for file in filenames:
                if file.endswith(".tar.gz") or file.endswith(".tgz"):
                    path = os.path.join(dirpath, file)
                    with tarfile.open(path, "r:gz") as data:
                        data.extractall(dirpath)

    else:
        os.rename(temp_path, out_path)


def download_data(base_url, out_dir, remote_name, check_local_name=None):
    """Downloads data files from a given url and unzip if it is a .tar.gz

    If check_local_names is provided, check if <out_dir>/<check_local_name[i]>
    exists first and don't download the file if it does.

    Args:
        base_url               (str): Url to download files from
        out_dir                (str): Directory to save files into
        remote_name      (list[str]): Name of files to download
        check_local_name (list[str]): Names of file to check for
    """

    for i, f_name in enumerate(remote_name):
        out_path = os.path.join(out_dir, f_name)
        if check_local_name is not None:
            check_path = os.path.join(out_dir, check_local_name[i])
            if os.path.exists(check_path):
                continue

        print('downloading', f_name)
        url = requests.compat.urljoin(base_url, f_name)
        _download_file(url, out_path)
