#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module defines functions for accessing locally available data files."""

import os

import numpy as np
from astropy.table import Table
from tqdm import tqdm

import _module_paths as paths
from data_access._utils import keep_restframe_bands

master_table = Table.read(paths.master_table_path, format='ascii')


def get_data_for_id(cid):
    """Returns published photometric data for a SDSS observed object

    No data cuts are applied to the returned data.

    Args:
        cid (int): The Candidate ID of the desired object

    Returns:
        An astropy table of photometric data for the given candidate ID
    """

    # Read in ascci data table for specified object
    file_path = os.path.join(paths.smp_dir, 'SMP_{:06d}.dat'.format(cid))
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


@np.vectorize
def sdss_mag_to_ab_flux(mag, band):
    """For a given sdss magnitude return the AB flux

    Args:
        mag (float): An SDSS asinh magnitude
        band  (str): The band of the magnitude doi_2010_<ugriz><123456>

    Return:
        The equivalent AB magnitude
    """

    if band[-2] == 'u':
        offset = -0.679

    elif band[-2] == 'g':
        offset = 0.0203

    elif band[-2] == 'r':
        offset = 0.0049

    elif band[-2] == 'i':
        offset = 0.0178

    elif band[-2] == 'z':
        offset = 0.0102

    else:
        ValueError('Unknown band {}'.format(band))

    return 3631 * 10 ** ((mag + offset) / -2.5)


def calc_err(sigma_sdss_mag, flux_ab):
    """Calculate the error of the AB magnitude equivilent for an SDSS asinh mag

    Args:
        sigma_sdss_mag (float): Error in the SDSS magnitude
        flux_ab        (float): AB flux of the measurement

    Returns:
        The error in the equivalent AB flux
    """

    return sigma_sdss_mag * flux_ab * np.log(10) / 2.5


@np.vectorize
def band_name(filt, idccd):
    """Return the sncosmo band name given filter and CCD id

    Args:
        filt  (str): Filter name <ugriz>
        idccd (int): Column number 1 through 6

    Args:
        The name of the filter registered with sncosmo
    """

    return 'doi_2010_{}{}'.format('ugriz'[filt], idccd)


def get_input_for_id(cid, bands=None):
    """Returns an SNCosmo input table a given SDSS object ID

    Only data points with a published photometric quality flag < 1024 are
    included in the returned table.

    Args:
        cid         (int): The ID of the desired object
        bands (iter[str]): Optionally only return select bands (eg. 'desg')

    Returns:
        An astropy table of photometric data formatted for use with SNCosmo
    """

    # Effective wavelengths for SDSS filters ugriz in angstroms
    # https://www.sdss.org/instruments/camera/#Filters
    sdss_bands = ('sdssu', 'sdssg', 'sdssr', 'sdssi', 'sdssz')
    lambda_effective = np.array([3551, 4686, 6166, 7480, 8932])

    # Format table
    phot_data = get_data_for_id(cid)
    phot_data = phot_data[phot_data['FLAG'] < 1024]

    sncosmo_table = Table()
    sncosmo_table.meta = phot_data.meta
    sncosmo_table['time'] = phot_data['MJD']
    sncosmo_table['band'] = band_name(phot_data['FILT'], phot_data['IDCCD'])
    sncosmo_table['zp'] = np.full(len(phot_data), 2.5 * np.log10(3631))
    sncosmo_table['flux'] = sdss_mag_to_ab_flux(phot_data['MAG'], sncosmo_table['band'])
    sncosmo_table['fluxerr'] = calc_err(phot_data['MERR'], sncosmo_table['flux'])
    sncosmo_table['zpsys'] = np.full(len(phot_data), 'ab')
    sncosmo_table.meta['cid'] = cid

    # Keep only specified band-passes
    if bands is not None:
        sncosmo_table = keep_restframe_bands(
            sncosmo_table, bands, sdss_bands, lambda_effective)

    return sncosmo_table


def iter_sncosmo_input(bands=None, skip_types=(), verbose=False):
    """Iterate through SDSS supernova and yield the SNCosmo input tables

    To return a select collection of band-passes, specify the band argument.
    Only data points with a published photometric quality flag < 1024 are
    included in the returned tables.

    Args:
        bands      (iter[str]): Optional list of band-passes to return
        skip_types (iter[str]): List of case sensitive classifications to skip
        verbose         (bool): Whether to a display progress bar while iterating

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
