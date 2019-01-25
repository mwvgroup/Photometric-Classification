#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module provides access to data from the year three DES supernova
cosmology paper.
"""

from os import path as _path

import numpy as np
from astropy.table import Table
from tqdm import tqdm

from ._download_data import download_data
from ._utils import keep_restframe_bands

DES_URL = 'http://desdr-server.ncsa.illinois.edu/despublic/sn_files/y3/tar_files/'
DATA_DIR = _path.join(_path.dirname(_path.realpath(__file__)), 'des_data')
FILT_DIR = _path.join(DATA_DIR, '01-FILTERS')
PHOT_DIR = _path.join(DATA_DIR, '02-DATA_PHOTOMETRY/DES-SN3YR_DES')
FITS_DIR = _path.join(DATA_DIR, '04-BBCFITS')

# Download data if it does not exist
download_data(DES_URL, DATA_DIR,
              ['01-FILTERS.tar.gz', '02-DATA_PHOTOMETRY.tar.gz', '04-BBCFITS.tar.gz'],
              [FILT_DIR, PHOT_DIR, FITS_DIR])


def get_data_for_id(obj_id):
    """Returns photometric data for a supernova candidate in a given filter

    Args:
        obj_id (int): The Candidate ID of the desired object

    Returns:
        An astropy table of photometric data for the given candidate ID
    """

    # Read in ascci data table for specified object
    file_path = _path.join(PHOT_DIR, 'des_{:08d}.dat'.format(obj_id))
    all_data = Table.read(
        file_path, format='ascii',
        data_start=34, data_end=-1,
        names=['VARLIST:', 'MJD', 'BAND', 'FIELD', 'FLUXCAL', 'FLUXCALERR',
               'ZPFLUX', 'PSF', 'SKYSIG', 'GAIN', 'PHOTFLAG', 'PHOTPROB'])

    # Add meta data to table
    with open(file_path) as ofile:
        meta_data = ofile.readlines()
        all_data.meta['ra'] = float(meta_data[7].split()[1])
        all_data.meta['dec'] = float(meta_data[8].split()[1])
        all_data.meta['PEAKMJD'] = float(meta_data[12].split()[1])
        all_data.meta['redshift'] = float(meta_data[13].split()[1])
        all_data.meta['redshift_err'] = float(meta_data[13].split()[3])
        del all_data.meta['comments']

    return all_data


def iter_sncosmo_input(bands=None, verbose=False):
    """Iterate through SDSS supernova and yield the SNCosmo input tables

    To return a select collection of band passes, specify the band argument.

    Args:
        bands   (list): Optional list of bandpasses to return
        verbose (bool): Whether to display a progress bar while iterating

    Yields:
        An astropy table formatted for use with SNCosmo
    """

    # Effective wavelengths taken from http://www.mso.anu.edu.au/~brad/filters.html
    des_bands = ('desg', 'desr', 'desi', 'desz', 'desy')
    lambda_effective = np.array([5270, 6590, 7890, 9760, 10030])

    # Load list of all target ids
    file_path = _path.join(PHOT_DIR, 'DES-SN3YR_DES.LIST')
    file_list = np.genfromtxt(file_path, dtype=str)

    # Yield an SNCosmo input table for each target
    iter_data = tqdm(file_list) if verbose else file_list
    for obj_id in iter_data:
        obj_id_int = int(obj_id.lstrip('des_').rstrip('.dat'))
        all_sn_data = get_data_for_id(obj_id_int)

        sncosmo_table = Table()
        sncosmo_table['time'] = all_sn_data['MJD']
        sncosmo_table['band'] = ['des' + s for s in all_sn_data['BAND']]
        sncosmo_table['flux'] = all_sn_data['FLUXCAL']
        sncosmo_table['fluxerr'] = all_sn_data['FLUXCALERR']
        sncosmo_table['zp'] = all_sn_data['ZPFLUX']
        sncosmo_table['zpsys'] = np.full(len(all_sn_data), 'ab')
        sncosmo_table.meta = all_sn_data.meta
        sncosmo_table.meta['obj_id'] = obj_id

        if bands is not None:
            sncosmo_table = keep_restframe_bands(
                sncosmo_table, bands, des_bands, lambda_effective)

        yield sncosmo_table
