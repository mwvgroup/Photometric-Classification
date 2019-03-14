#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module defines functions for accessing locally available data files."""

from os import path as _path

import numpy as np
from astropy.table import Table
from tqdm import tqdm

from . import _module_meta_data as meta_data
from .._utils import keep_restframe_bands, parse_snoopy_data

master_table = Table.read(
    meta_data.master_path,
    names=['SN', 'N', 'epoch_range', 'SNID', 'Wang', 'Branch',
           'zHelio', 'zCMB', 'T0', 'Bmag', 'e_Bmag', 'Dm15',
           'eDm15', 'E(B-V)', 'B-V', 'e_B-V', 'low_reddening',
           'DBMag', 'e_DBMag'],
    format='ascii')


def get_data_for_id(cid):
    """Returns photometric data for a supernova candidate in a given filter

    No data cuts are applied to the returned data.

    Args:
        cid (int): The Candidate ID of the desired object

    Returns:
        An astropy table of photometric data for the given candidate ID
    """

    file_path = _path.join(meta_data.photometry_dir, f'SN{cid}_snpy.txt')
    data_table = parse_snoopy_data(file_path)
    data_table['band'] = '91bg_proj_csp_' + data_table['band']
    data_table.meta['cid'] = cid

    return data_table


def _get_zp_for_band(band):
    sorter = np.argsort(meta_data.band_names)
    indices = sorter[np.searchsorted(meta_data.band_names, band, sorter=sorter)]
    return np.array(meta_data.zero_point)[indices]


def get_input_for_id(cid, bands=None):
    """Returns an SNCosmo input table a given CSP object ID

    No data cuts are applied to the returned data.

    Args:
        cid         (int): The ID of the desired object
        bands (list[str]): Optionally only return select bands
                             (eg. '91bg_proj_csp_V0')

    Returns:
        An astropy table of photometric data formatted for use with SNCosmo
    """

    sn_data = get_data_for_id(cid)
    sn_data['zp'] = _get_zp_for_band(sn_data['band'])
    sn_data['zpsys'] = np.full(len(sn_data), 'ab')
    sn_data['flux'] = 10 ** ((sn_data['mag'] - sn_data['zp']) / -2.5)
    sn_data['fluxerr'] = np.log(10) * sn_data['flux'] * sn_data['mag_err'] / 2.5
    sn_data.remove_columns(['mag', 'mag_err'])

    if bands is not None:
        sn_data = keep_restframe_bands(
            sn_data, bands, meta_data.band_names, meta_data.lambda_effective)

    return sn_data


def iter_sncosmo_input(bands=None, verbose=False):
    """Iterate through CSP supernova and yield the SNCosmo input tables

    To return a select collection of band-passes, specify the band argument.
    No data cuts are applied to the returned data.

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
