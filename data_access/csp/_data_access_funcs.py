#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module defines functions for accessing locally available data files."""

from os import path as _path

import numpy as np
from astropy.table import Table
from tqdm import tqdm

from . import _module_paths as paths
from .._utils import keep_restframe_bands, parse_snoopy_data

master_table = Table.read(paths.master_path, format='ascii')


def get_data_for_id(cid):
    """Returns photometric data for a supernova candidate in a given filter

    Args:
        cid (int): The Candidate ID of the desired object

    Returns:
        An astropy table of photometric data for the given candidate ID
    """

    file_path = _path.join(paths.photometry_dir, f'{cid}_snpy.txt')
    return parse_snoopy_data(file_path)


def get_input_for_id(cid, bands=None):
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

    To return a select collection of band-passes, specify the band argument.

    Args:
        bands (iter[str]): Optional list of band-passes to return
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
