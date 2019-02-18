#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module defines functions for accessing locally available data files."""

import os


import numpy as np
from astropy.table import Table
from tqdm import tqdm

from . import _module_meta_data as meta_data
from .._utils import keep_restframe_bands

master_table = Table.read(meta_data.master_table_path, format='ascii')


@np.vectorize
def construct_band_name(filter_id, ccd_id):
    """Return the sncosmo band name given filter and CCD id

    Args:
        filter_id (int): Filter index 1 through 5 for 'ugriz'
        ccd_id    (int): Column number 1 through 6

    Args:
        The name of the filter registered with sncosmo
    """

    return f'doi_2010_{"ugriz"[filter_id]}{ccd_id}'


def get_data_for_id(cid):
    """Returns published photometric data for a SDSS observed object

    No data cuts are applied to the returned data.

    Args:
        cid (int): The Candidate ID of the desired object

    Returns:
        An astropy table of photometric data for the given candidate ID
    """

    # Read in ascii data table for specified object
    file_path = os.path.join(meta_data.paths.smp_dir, f'SMP_{cid:06d}.dat')
    all_data = Table.read(file_path, format='ascii')

    # Rename columns using header data from file
    col_names = all_data.meta['comments'][-1].split()
    for i, name in enumerate(col_names):
        all_data[f'col{i + 1}'].name = name

    meta_data = master_table[master_table['CID'] == cid]
    all_data.meta['redshift'] = meta_data['zCMB'][0]
    all_data.meta['ra'] = meta_data['RA'][0]
    all_data.meta['dec'] = meta_data['DEC'][0]
    all_data.meta['classification'] = meta_data['Classification'][0]
    all_data.meta['name'] = meta_data['IAUName'][0]

    return all_data


def get_input_for_id(cid, bands=None):
    """Returns an SNCosmo input table a given SDSS object ID

    Only data points with a published photometric quality flag < 1024 are
    included in the returned table. Data is dropped for epochs not between
    -15 and 45 days.

    Args:
        cid         (int): The ID of the desired object
        bands (iter[str]): Optionally only return select bands (eg. 'desg')

    Returns:
        An astropy table of photometric data formatted for use with SNCosmo
    """

    # Format table
    phot_data = get_data_for_id(cid)
    peak_mjd = master_table[master_table['CID'] == cid]['MJDatPeakrmag'][0]
    phot_data = phot_data[phot_data['FLAG'] < 1024]
    phot_data = phot_data[phot_data['MJD'] < peak_mjd + 45]
    phot_data = phot_data[phot_data['MJD'] > peak_mjd - 15]
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
    sncosmo_table.meta['cid'] = cid

    # Keep only specified band-passes
    if bands is not None:
        sncosmo_table = keep_restframe_bands(
            sncosmo_table,
            bands,
            meta_data.band_names,
            meta_data.lambda_effective)

    return sncosmo_table


def iter_sncosmo_input(bands=None, keep_types=(), verbose=False):
    """Iterate through SDSS supernova and yield the SNCosmo input tables

    To return a select collection of band-passes, specify the band argument.
    Only data points with a published photometric quality flag < 1024 are
    included in the returned tables. Data is dropped for epochs not between
    -15 and 45 days.

    Args:
        bands      (iter[str]): Optional list of band-passes to return
        keep_types (iter[str]): List of case sensitive classifications to skip
        verbose         (bool): Whether to a display progress bar while iterating

    Yields:
        An astropy table formatted for use with SNCosmo
    """

    # Create iterable without unwanted data
    skip_data_indx = np.isin(master_table['Classification'], keep_types)
    cut_data = master_table[skip_data_indx]

    # Yield an SNCosmo input table for each target
    iter_data = tqdm(cut_data['CID']) if verbose else cut_data['CID']
    for cid in iter_data:
        sncosmo_table = get_input_for_id(cid, bands)
        if sncosmo_table:
            yield sncosmo_table
