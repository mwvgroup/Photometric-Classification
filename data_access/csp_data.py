#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module provides access to data from the third data release of the
Carnegie Supernova Project (CSP).

For more information on CSP see: https://csp.obs.carnegiescience.edu
"""

from os import path as _path

import numpy as np
from astropy.table import Table
from tqdm import tqdm

from ._download_data import download_data
from ._utils import keep_restframe_bands, parse_snoopy_data

CSP_URL = 'https://csp.obs.carnegiescience.edu/data/'

# Define local paths of CSP data
DATA_DIR = _path.join(_path.dirname(_path.realpath(__file__)), 'csp_data')
PHOT_DIR = _path.join(DATA_DIR, 'DR3')                 # DR3 Light Curves
FILT_DIR = _path.join(DATA_DIR, 'CSP_filter_package')  # DR3 Light Curves
MASTER_PTH = _path.join(PHOT_DIR, 'tab1.dat')          # Master table

# Download data if it does not exist
download_data(
    base_url=CSP_URL,
    out_dir=DATA_DIR,
    remote_name=['CSP_Photometry_DR3.tgz', 'CSP_filter_package.tgz'],
    check_local_name=[PHOT_DIR, FILT_DIR])

master_table = Table.read(MASTER_PTH, format='ascii')


def get_data_for_id(cid):
    """Returns photometric data for a supernova candidate in a given filter

    Args:
        cid (int): The Candidate ID of the desired object

    Returns:
        An astropy table of photometric data for the given candidate ID
    """

    file_path = _path.join(PHOT_DIR, '{}_snpy.txt'.format(cid))
    return parse_snoopy_data(file_path)


def get_input_for_id(cid, bands):
    """Returns an SNCosmo input table a given SDSS object ID

    Args:
        cid         (int): The ID of the desired object
        bands (list[str]): Optionally only return select bands (eg. 'desg')

    Returns:
        An astropy table of photometric data formatted for use with SNCosmo
    """

    # Todo: Add and register CSP bands
    csp_bands = []
    lambda_effective = []

    sn_data = get_data_for_id(cid)
    sn_data['flux'] = sn_data['mag']
    sn_data['fluxerr'] = sn_data['mag_err']
    sn_data['zp'] = np.full(len(sn_data), 27.5)
    sn_data['zpsys'] = np.full(len(sn_data), 'ab')
    sn_data.remove_columns(['mag', 'mag_err'])

    if bands is not None:
        sn_data = keep_restframe_bands(
            sn_data, bands, csp_bands, lambda_effective)

    return sn_data


def iter_sncosmo_input(bands=None, verbose=False):
    """Iterate through SDSS supernova and yield the SNCosmo input tables

    To return a select collection of band passes, specify the band argument.

    Args:
        bands      (list): Optional list of band passes to return
        verbose    (bool): Whether to display a progress bar while iterating

    Yields:
        An astropy table formatted for use with SNCosmo
    """

    # Yield an SNCosmo input table for each target
    iter_data = tqdm(master_table['SN']) if verbose else master_table['SN']
    for target_name in iter_data:
        sncosmo_table = get_input_for_id(target_name, bands)
        if sncosmo_table:
            yield sncosmo_table
