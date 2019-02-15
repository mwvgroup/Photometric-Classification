#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module defines functions for accessing locally available data files."""

from os import path as _path

import numpy as np
from astropy.table import Table
from tqdm import tqdm

import _module_paths as paths
from data_access._utils import keep_restframe_bands

# Effective wavelengths taken from
# http://www.mso.anu.edu.au/~brad/filters.html
band_names = ('desg', 'desr', 'desi', 'desz', 'desy')
lambda_effective = (5270, 6590, 7890, 9760, 10030)

master_table = Table.read(
    paths.master_table_path,
    format='ascii',
    data_start=4,
    comment='#',
    exclude_names=['dummy_col'],
    names=['dummy_col', 'CID', 'CIDint', 'IDSURVEY', 'TYPE', 'FIELD',
           'CUTFLAG_SNANA', 'zHEL', 'zHELERR', 'zCMB', 'zCMBERR',
           'zHD', 'zHDERR', 'VPEC', 'VPECERR', 'HOST_LOGMASS',
           'HOST_LOGMASS_ERR', 'SNRMAX1', 'SNRMAX2', 'SNRMAX3', 'PKMJD',
           'PKMJDERR', 'x1', 'x1ERR', 'c', 'cERR', 'mB', 'mBERR', 'x0',
           'x0ERR', 'COV_x1_c', 'COV_x1_x0', 'COV_c_x0', 'NDOF',
           'FITCHI2', 'FITPROB', 'RA', 'DECL', 'TGAPMAX', 'TrestMIN',
           'TrestMAX', 'MWEBV', 'm0obs_i', 'm0obs_r', 'em0obs_i', 'em0obs_r',
           'MU', 'MUMODEL', 'MUERR', 'MUERR_RAW', 'MURES', 'MUPULL', 'M0DIF',
           'ERRCODE', 'biasCor_mu', 'biasCorErr_mu', 'biasCor_mB',
           'biasCor_x1', 'biasCor_c', 'biasScale_muCOV', 'IDSAMPLE'])


def get_data_for_id(cid):
    """Returns DES photometric data for a given object ID

    Args:
        cid (int): The ID of the desired object

    Returns:
        An astropy table of photometric data for the given candidate ID
    """

    # Read in ascci data table for specified object
    file_path = _path.join(paths.photometry_dir, 'des_{:08d}.dat'.format(cid))
    all_data = Table.read(
        file_path, format='ascii',
        data_start=27, data_end=-1,
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


def get_input_for_id(cid, bands=None):
    """Returns an SNCosmo input table a given DES object ID

    Args:
        cid         (int): The ID of the desired object
        bands (list[str]): Optionally only return select bands (eg. 'desg')

    Returns:
        An astropy table of photometric data formatted for use with SNCosmo
    """

    all_sn_data = get_data_for_id(cid)
    sncosmo_table = Table()
    sncosmo_table['time'] = all_sn_data['MJD']
    sncosmo_table['band'] = ['des' + s for s in all_sn_data['BAND']]
    sncosmo_table['flux'] = all_sn_data['FLUXCAL']
    sncosmo_table['fluxerr'] = all_sn_data['FLUXCALERR']
    sncosmo_table['zp'] = np.full(len(all_sn_data), 27.5)
    sncosmo_table['zpsys'] = np.full(len(all_sn_data), 'ab')
    sncosmo_table.meta = all_sn_data.meta
    sncosmo_table.meta['cid'] = cid

    if bands is not None:
        sncosmo_table = keep_restframe_bands(
            sncosmo_table, bands, band_names, lambda_effective)

    return sncosmo_table


def iter_sncosmo_input(bands=None, verbose=False):
    """Iterate through SDSS supernova and yield the SNCosmo input tables

    To return a select collection of band-passes, specify the band argument.

    Args:
        bands (iter[str]): Optional list of band-passes to return
        verbose    (bool): Whether to display a progress bar while iterating

    Yields:
        An astropy table formatted for use with SNCosmo
    """

    # Load list of all target ids
    target_list_path = _path.join(paths.photometry_dir, 'DES-SN3YR_DES.LIST')
    file_list = np.genfromtxt(target_list_path, dtype=str)

    # Yield an SNCosmo input table for each target
    iter_data = tqdm(file_list) if verbose else file_list
    for file_name in iter_data:
        cid_int = int(file_name.lstrip('des_').rstrip('.dat'))
        sncosmo_table = get_input_for_id(cid_int, bands)
        if sncosmo_table:
            yield sncosmo_table
