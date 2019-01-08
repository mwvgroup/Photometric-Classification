#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module fits SDSS light curves using sncosmo"""

import os

import numpy as np
import sncosmo
from astropy.table import Table
from tqdm import tqdm

import sys; sys.path.append('../')
from parse_sn_data import get_cid_data,master_table


@np.vectorize
def band_index_mapping(i):
    """Vectorized mapping between index and filter name sdss<"ugriz"[i]>

    Args:
        i (int or str): Index of bandpass or bandpass name

    Returns:
        If the argument is an integer, return the filter name sdss<"ugriz"[i]>
        If the argument is a filtername, return the index
    """

    if isinstance(i, int):
        return 'sdss' + 'ugriz'[i]

    elif isinstance(i, str):
        return ['sdssu', 'sdssg', 'sdssr', 'sdssi', 'sdssz'].index(i)

    else:
        raise ValueError('Argument must be int or str')


def keep_restframe_bands(data_table, bands):
    """Return select rest-frame bandpasses from an SNCsomo input table

    Args:
        data_table (Table): An SNCosmo input table as given
                              by iter_sncosmo_input
        bands       (list): List of rest-frame bandpasses to keep

    Returns:
        A new input table for SNCosmo only containing select rest frame bands
    """

    # Effective wavelengths for SDSS filters ugriz in angstroms
    effective_lambda = np.array([3551, 4686, 6166, 7480, 8932])

    # Rest frame effective wavelengths for each observation
    z = data_table.meta['redshift']
    filter_indices = band_index_mapping(data_table['band'])
    rest_frame_lambda = effective_lambda[filter_indices] / (1 + z)

    # Find rest frame filter names
    band_index_array = np.abs(effective_lambda - rest_frame_lambda).argmin()
    rest_frame_filters = band_index_mapping(band_index_array)

    # Keep only the specified filters
    indices = np.isin(rest_frame_filters, bands)
    return data_table[indices]


def iter_sncosmo_input(bands=None):
    """Iterate through SDSS supernova and yield the SNCosmo input tables

    To return a select collection of band passes, specify the band argument.

    Args:
        bands (list): Optional list of bandpasses to return

    Yields:
        An astropy table formatted for use with SNCosmo
    """

    has_redshift = master_table['zspecHelio'] != -9
    iter_data = master_table[has_redshift]
    for cid in tqdm(iter_data['CID']):
        all_sn_data = get_cid_data(cid)

        sncosmo_table = Table()
        sncosmo_table['time'] = all_sn_data['MJD']
        sncosmo_table['band'] = band_index_mapping(all_sn_data['FILT'])
        sncosmo_table['flux'] = all_sn_data['FLUX']
        sncosmo_table['fluxerr'] = all_sn_data['FLUXERR']
        sncosmo_table['zp'] = np.full(len(all_sn_data), 25)
        sncosmo_table['zpsys'] = np.full(len(all_sn_data), 'ab')
        sncosmo_table.meta = all_sn_data.meta
        sncosmo_table.meta['cid'] = cid

        if bands is not None:
            sncosmo_table = keep_restframe_bands(sncosmo_table, bands)

        yield sncosmo_table


def run_fit_for_object(input_table, model_name, params_to_fit):
    """Fit a single light curve with SNCosmo

    Args:
        input_table  (Table): An SNCosmo input table
        model_name     (str): The name of the SNCosmo model to use
        params_to_fit (list): List of parameters to fit

    Returns:
        The output dictionary from SNCosmo with fitting results
    """

    # Configure model for fitting
    model = sncosmo.Model(source=model_name)
    z = input_table.meta['redshift']
    params_to_fit = list(params_to_fit)

    # Tell SNCosmo to fit for the redshift if it is not given
    if z == -9.:
        params_to_fit.insert(0, 'z')
        bounds = {'z': (0.002, 1)}

    else:
        model.set(z=z)
        bounds = None

    # Run fit
    result, fitted_model = sncosmo.fit_lc(
        input_table, model, params_to_fit, bounds=bounds)

    return result


def fit_sdss_data(out_path, model_name='salt2', bands=None,
                  params_to_fit=('t0', 'x0', 'x1', 'c')):
    """Fit SDSS light curves with SNCosmo

    Files are named as <out_dir>/<target cid>.txt

    Args:
        out_path       (str): Where to write fit results
        model_name     (str): Model to use for fitting. Default = salt2
        params_to_fit (list): List of parameters to fit
        bands         (list): Optional list of bandpasses to fit
    """

    out_dir = os.path.dirname(out_path)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # To store fit results
    names = ['cid', 'num_bands', 'num_points', 'class', 'fit_z', 'z']
    names.extend(params_to_fit)
    names.append('z_err')
    names.extend((v + '_err' for v in params_to_fit))
    names.extend(('chi', 'dof', 'message'))
    out_table = Table(names=names, dtype=[object for _ in names])

    for input_table in iter_sncosmo_input(bands=bands):
        z_was_fit = int(input_table.meta['redshift'] == -9)
        new_row = [input_table.meta['cid'],
                   len(set(input_table['band'])),
                   len(input_table),
                   input_table.meta['classification'],
                   z_was_fit]

        try:
            result = run_fit_for_object(input_table, model_name, params_to_fit)

        except (sncosmo.fitting.DataQualityError, RuntimeError, ValueError) as e:
            mask = np.full(len(out_table.colnames) - len(new_row) - 1, np.NAN).tolist()
            new_row.extend(mask)
            new_row.append(str(e))

        else:
            new_row.extend(result.parameters)
            if 'z' not in result.errors:
                new_row.append(np.NAN)

            new_row.extend(result.errors.values())
            new_row.append(result.chisq)
            new_row.append(result.ndof)
            new_row.append(result.message)

        out_table.add_row(new_row)
        out_table.write(out_path, overwrite=True)


if __name__ == '__main__':
    print('Fitting type Ia model in all bands')
    fit_sdss_data('./sncosmo_results/snia_ugriz.csv')

    print('\n\nFitting type Ia model in ug')
    fit_sdss_data('./sncosmo_results/snia_ug.csv',
                  bands=['sdssu', 'sdssg'])

    print('\n\nFitting type Ia model in riz')
    fit_sdss_data('./sncosmo_results/snia_riz.csv',
                  bands=['sdssr', 'sdssi', 'sdssz'])

    # 91bg model
    print('\n\nFitting 91bg model in all bands')
    fit_sdss_data('./sncosmo_results/91bg_ugriz.csv',
                  model_name='nugent-sn91bg',
                  params_to_fit=['t0', 'amplitude'])

    print('\n\nFitting 91bg model in ug')
    fit_sdss_data('./sncosmo_results/91bg_ug.csv',
                  model_name='nugent-sn91bg',
                  bands=['sdssu', 'sdssg'],
                  params_to_fit=['t0', 'amplitude'])

    print('\n\nFitting 91bg model in riz')
    fit_sdss_data('./sncosmo_results/91bg_riz.csv',
                  model_name='nugent-sn91bg',
                  bands=['sdssr', 'sdssi', 'sdssz'],
                  params_to_fit=['t0', 'amplitude'])


