#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module provides access to data from the SDSS-II SN Catalog Data
Release.
"""

from os import path as _path

import numpy as np
from astropy.table import Table
from tqdm import tqdm

from ._download_data import download_data
from ._utils import keep_restframe_bands

# Define local paths of SDSS data
DATA_DIR = _path.join(_path.dirname(_path.realpath(__file__)), 'sdss_data')
SNANA_FILES = _path.join(DATA_DIR, 'SDSS_dataRelease-snana/')  # SNANA files
MASTER_PTH = _path.join(DATA_DIR, 'master_data.txt')           # Master table
SMP_DIR = _path.join(DATA_DIR, 'SMP_Data/')                    # SMP data files

# Download data if it does not exist
SDSS_URL = 'https://data.sdss.org/sas/dr10/boss/papers/supernova/'
download_data(
    base_url=SDSS_URL,
    out_dir=DATA_DIR,
    remote_name=['master_data.txt', 'SMP_Data.tar.gz', 'SDSS_dataRelease-snana.tar.gz'],
    check_local_name=[MASTER_PTH, SMP_DIR, SNANA_FILES])

master_table = Table.read(MASTER_PTH, format='ascii')


def get_data_for_id(cid):
    """Returns photometric data for a supernova candidate in a given filter

    Args:
        cid (int): The Candidate ID of the desired object

    Returns:
        An astropy table of photometric data for the given candidate ID
    """

    # Read in ascci data table for specified object
    file_path = _path.join(SMP_DIR, 'SMP_{:06d}.dat'.format(cid))
    all_data = Table.read(file_path, format='ascii')

    # Rename columns using header data from file
    col_names = all_data.meta['comments'][-1].split()
    for i, name in enumerate(col_names):
        all_data['col{}'.format(i + 1)].name = name

    meta_data = master_table[master_table['CID'] == cid]
    all_data.meta['redshift'] = meta_data['zspecHelio'][0]
    all_data.meta['ra'] = meta_data['RA'][0]
    all_data.meta['dec'] = meta_data['DEC'][0]
    all_data.meta['classification'] = meta_data['Classification'][0]
    all_data.meta['name'] = meta_data['IAUName'][0]

    return all_data


def get_input_for_id(cid, bands):
    """Returns an SNCosmo input table a given SDSS object ID

    Args:
        cid         (int): The ID of the desired object
        bands (list[str]): Optionally only return select bands (eg. 'desg')

    Returns:
        An astropy table of photometric data formatted for use with SNCosmo
    """

    # Effective wavelengths for SDSS filters ugriz in angstroms
    # https://www.sdss.org/instruments/camera/#Filters
    sdss_bands = ('sdssu', 'sdssg', 'sdssr', 'sdssi', 'sdssz')
    lambda_effective = np.array([3551, 4686, 6166, 7480, 8932])

    # Format table
    all_sn_data = get_data_for_id(cid)
    sncosmo_table = Table()
    sncosmo_table['time'] = all_sn_data['MJD']
    sncosmo_table['band'] = [sdss_bands[i] for i in all_sn_data['FILT']]
    sncosmo_table['flux'] = all_sn_data['FLUX']
    sncosmo_table['fluxerr'] = all_sn_data['FLUXERR']
    sncosmo_table['zp'] = np.full(len(all_sn_data), 25)
    sncosmo_table['zpsys'] = np.full(len(all_sn_data), 'ab')
    sncosmo_table.meta = all_sn_data.meta
    sncosmo_table.meta['cid'] = cid

    # Keep only specified band-passes
    if bands is not None:
        sncosmo_table = keep_restframe_bands(
            sncosmo_table, bands, sdss_bands, lambda_effective)

    return sncosmo_table


def iter_sncosmo_input(bands=None, skip_types=(), verbose=False):
    """Iterate through SDSS supernova and yield the SNCosmo input tables

    To return a select collection of band passes, specify the band argument.

    Args:
        bands      (list): Optional list of bandpasses to return
        skip_types (list): List of case sensitive classifications to skip
        verbose    (bool): Whether to display a progress bar while iterating

    Yields:
        An astropy table formatted for use with SNCosmo
    """

    # Create iterable without unwanted data
    skip_data_indx = np.isin(master_table['Classification'], skip_types)
    cut_data = master_table[np.logical_not(skip_data_indx)]

    # Yield an SNCosmo input table for each target
    iter_data = tqdm(cut_data['CID']) if verbose else cut_data['CID']
    for cid in iter_data:
        sncosmo_table = get_input_for_id(cid, bands)
        if sncosmo_table:
            yield sncosmo_table
