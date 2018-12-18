#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This document defines functions for parsing data from the SDSS-II SN Catalog
Data Release.
"""

import os

from astropy.table import Table

from download_data import download_sdss_data

FILE_DIR = os.path.dirname(os.path.realpath(__file__))
MASTER_PTH = os.path.join(FILE_DIR, 'data/master_table.txt')  # Master table path
SMP_DIR = os.path.join(FILE_DIR, 'data/SMP_Data/')            # SMP data files

# Download data if not available
for path in (MASTER_PTH, SMP_DIR):
    if not os.path.exists(path):
        print('SDSS data not found. Downloading it to ./data/')
        data_dir = os.path.join(FILE_DIR, 'data/')
        download_sdss_data(data_dir)

MASTER_TABLE = Table.read(MASTER_PTH, format='ascii')


def get_cid_data(cid, filter_name=None):
    """Returns photometric data for a supernova candidate in a given filter

    Args:
        cid         (int): The Candidate ID of the desired object
        filter_name (str): Optionally return data only for a given filter ugriz

    Returns:
        An astropy table of photometric data for the given candidate ID
    """

    # Read in ascci data table for specified object
    file_path = os.path.join(SMP_DIR, 'SMP_{:06d}.dat'.format(cid))
    all_data = Table.read(file_path, format='ascii')

    # Rename columns using header data from file
    col_names = all_data.meta['comments'][-1].split()
    for i, name in enumerate(col_names):
        all_data['col{}'.format(i + 1)].name = name

    if filter_name is not None:
        filter_index = 'ugriz'.index(filter_name)
        all_data = all_data[all_data['FILT'] == filter_index]
        all_data.remove_column('FILT')

    meta_data = MASTER_TABLE[MASTER_TABLE['CID'] == cid]
    all_data.meta['redshift'] = meta_data['zspecHelio'][0]
    all_data.meta['ra'] = meta_data['RA'][0]
    all_data.meta['dec'] = meta_data['DEC'][0]

    return all_data
