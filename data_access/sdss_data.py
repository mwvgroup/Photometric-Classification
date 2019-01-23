#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module provides access to data from the SDSS-II SN Catalog Data
Release.
"""

from os import path as _path

from astropy.table import Table

from _download_data import download_data

SDSS_URL = 'https://data.sdss.org/sas/dr10/boss/papers/supernova/'
DATA_DIR = _path.join(_path.dirname(_path.realpath(__file__)), 'sdss_data')
MASTER_PTH = _path.join(DATA_DIR, 'master_table.txt')  # Master table path
SMP_DIR = _path.join(DATA_DIR, 'SMP_Data/')  # SMP data files

# Download data if it does not exist
download_data(SDSS_URL, DATA_DIR, ['master_data.txt',
                                   'SMP_Data.tar.gz',
                                   'SDSS_dataRelease-snana.tar.gz'])

master_table = Table.read(MASTER_PTH, format='ascii')


def get_cid_data(cid, filter_name=None):
    """Returns photometric data for a supernova candidate in a given filter

    Args:
        cid         (int): The Candidate ID of the desired object
        filter_name (str): Optionally return data only for a given filter ugriz

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

    if filter_name is not None:
        filter_index = 'ugriz'.index(filter_name)
        all_data = all_data[all_data['FILT'] == filter_index]
        all_data.remove_column('FILT')

    meta_data = master_table[master_table['CID'] == cid]
    all_data.meta['redshift'] = meta_data['zspecHelio'][0]
    all_data.meta['ra'] = meta_data['RA'][0]
    all_data.meta['dec'] = meta_data['DEC'][0]
    all_data.meta['classification'] = meta_data['Classification'][0]
    all_data.meta['name'] = meta_data['IAUName'][0]

    return all_data
