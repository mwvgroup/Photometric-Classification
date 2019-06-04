#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module defines functions for accessing locally available data files."""

import os

import numpy as np
from astropy.table import Column, Table
from tqdm import tqdm

from . import _module_meta_data as meta_data
from .._utils import keep_restframe_bands

master_table = Table.read(meta_data.master_table_path, format='ascii')
master_table['obj_id'] = Column(master_table['CID'], dtype=str)


def _get_outliers():
    """Return a dictionary of data points marked by SDSS II as outliers

    Returns:
        A dictionary {<obj_id>: [<MJD of bad data point>, ...], ...}
    """

    out_dict = dict()
    with open(meta_data.outlier_path) as ofile:
        for line in ofile.readlines():
            if line.startswith('IGNORE:'):
                line_list = line.split()
                obj_id, mjd, band = line_list[1], line_list[2], line_list[3]
                if obj_id not in out_dict:
                    out_dict[str(obj_id)] = []

                out_dict[str(obj_id)].append(mjd)

    return out_dict


outlier_mjd = _get_outliers()


@np.vectorize
def construct_band_name(filter_id, ccd_id):
    """Return the sncosmo band name given filter and CCD id

    Args:
        filter_id (int): Filter index 1 through 5 for 'ugriz'
        ccd_id    (int): Column number 1 through 6

    Args:
        The name of the filter registered with sncosmo
    """

    return f'91bg_proj_sdss_{"ugriz"[filter_id]}{ccd_id}'


def get_data_for_id(obj_id):
    """Returns published photometric data for a SDSS observed object

    No data cuts are applied to the returned data.

    Args:
        obj_id (str): The Candidate ID of the desired object

    Returns:
        An astropy table of photometric data for the given candidate ID
    """

    # Read in ascii data table for specified object
    file_path = os.path.join(meta_data.smp_dir, f'SMP_{int(obj_id):06d}.dat')
    all_data = Table.read(file_path, format='ascii')

    # Rename columns using header data from file
    col_names = all_data.meta['comments'][-1].split()
    for i, name in enumerate(col_names):
        all_data[f'col{i + 1}'].name = name

    table_meta_data = master_table[master_table['obj_id'] == obj_id]
    all_data.meta['redshift'] = table_meta_data['zCMB'][0]
    all_data.meta['redshift_err'] = table_meta_data['zerrCMB'][0]
    all_data.meta['ra'] = table_meta_data['RA'][0]
    all_data.meta['dec'] = table_meta_data['DEC'][0]
    all_data.meta['classification'] = table_meta_data['Classification'][0]
    all_data.meta['name'] = table_meta_data['IAUName'][0]

    return all_data


def get_input_for_id(obj_id, bands=None):
    """Returns an SNCosmo input table a given SDSS object ID

    Only data points with a published photometric quality flag < 1024 are
    included in the returned table. Data points flagged in the SDSS II release
    as outliers are also removed.

    Args:
        obj_id         (str): The ID of the desired object
        bands (iter[str]): Optionally only return select rest frame
                             bands (eg. '91bg_proj_sdss_u1')

    Returns:
        An astropy table of photometric data formatted for use with SNCosmo
    """

    # Format table
    obj_id = str(obj_id)
    phot_data = get_data_for_id(obj_id)
    phot_data = phot_data[phot_data['FLAG'] < 1024]

    outlier_list = outlier_mjd.get(obj_id, [])
    if outlier_list:
        keep_indices = ~np.isin(phot_data['MJD'], outlier_list)
        phot_data = phot_data[keep_indices]

    if not phot_data:
        return Table(names=['time', 'band', 'zp', 'flux', 'fluxerr', 'zpsys'])

    sncosmo_table = Table()
    sncosmo_table.meta = phot_data.meta
    sncosmo_table['time'] = phot_data['MJD']
    sncosmo_table['band'] = construct_band_name(phot_data['FILT'], phot_data['IDCCD'])
    sncosmo_table['zp'] = np.full(len(phot_data), 2.5 * np.log10(3631))
    sncosmo_table['flux'] = phot_data['FLUX'] * 1E-6
    sncosmo_table['fluxerr'] = phot_data['FLUXERR'] * 1E-6
    sncosmo_table['zpsys'] = np.full(len(phot_data), 'ab')
    sncosmo_table.meta['obj_id'] = obj_id

    # Keep only specified band-passes
    if bands is not None:
        sncosmo_table = keep_restframe_bands(
            sncosmo_table,
            bands,
            meta_data.band_names,
            meta_data.lambda_effective)

    return sncosmo_table


def get_target_ids(keep_types=(), skip_types=()):
    """Return a list of object ID values

    Args:
        keep_types (iter[str]): Optional case sensitive classifications to keep
        skip_types (iter[str]): Optional case sensitive classifications to skip

    Returns:
        A list of object ID values
    """

    data = master_table
    if keep_types:
        keep_data_indx = np.isin(master_table['Classification'], keep_types)
        data = data[keep_data_indx]

    if skip_types:
        skip_data_indx = ~np.isin(master_table['Classification'], skip_types)
        data = data[skip_data_indx]

    return list(data['obj_id'])


def iter_sncosmo_input(bands=None, keep_types=(), skip_types=(), verbose=False):
    """Iterate through SDSS supernova and yield the SNCosmo input tables

    To return a select collection of band-passes, specify the band argument.
    Only data points with a published photometric quality flag < 1024 are
    included in the returned tables.

    Args:
        bands      (iter[str]): Optional list of band-passes to return
        keep_types (iter[str]): Optional case sensitive classifications to keep
        skip_types (iter[str]): Optional case sensitive classifications to skip
        verbose         (bool): Optionally display progress bar while iterating

    Yields:
        An astropy table formatted for use with SNCosmo
    """

    # Get list of IDS without unwanted types
    ids = get_target_ids(keep_types, skip_types)

    # Customize iterable of data
    if isinstance(verbose, dict):
        iter_data = tqdm(ids, **verbose)

    elif verbose:
        iter_data = tqdm(ids)

    else:
        iter_data = ids

    # Yield an SNCosmo input table for each target
    for obj_id in iter_data:
        sncosmo_table = get_input_for_id(obj_id, bands)
        if sncosmo_table:
            yield sncosmo_table
