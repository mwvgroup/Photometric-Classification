#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module fits SDSS light curves using sncosmo"""

import os

import numpy as np
import sncosmo
from astropy.table import Table
from tqdm import tqdm

from parse_sn_data import get_cid_data, MASTER_TABLE


@np.vectorize
def map_index_to_band(i):
    return 'sdss' + 'ugriz'[i]


def iter_sncosmo_input(bands=None):
    """Iterate through SDSS supernova and yield the SNCosmo input tables

    Yields:
        An astropy table formatted for use with SNCosmo
    """

    for cid in tqdm(MASTER_TABLE['CID']):
        sncosmo_table = Table()
        all_sn_data = get_cid_data(cid)

        sncosmo_table['time'] = all_sn_data['MJD']
        sncosmo_table['band'] = map_index_to_band(all_sn_data['FILT'])
        sncosmo_table['flux'] = all_sn_data['FLUX']
        sncosmo_table['fluxerr'] = all_sn_data['FLUXERR']
        sncosmo_table['zp'] = np.full(len(all_sn_data), 25)
        sncosmo_table['zpsys'] = np.full(len(all_sn_data), 'ab')
        sncosmo_table.meta = all_sn_data.meta
        sncosmo_table.meta['cid'] = cid

        if bands is not None:
            indices = np.isin(sncosmo_table['band'], bands)
            sncosmo_table = sncosmo_table[indices]

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


def fit_sdss_data(out_dir, model_name='salt2', bands=None,
                  params_to_fit=('t0', 'x0', 'x1', 'c')):
    """Fit SDSS light curves with SNCosmo

    Files are named as <out_dir>/<target cid>.txt

    Args:
        out_dir        (str): Where to write fit results
        model_name     (str): Model to use for fitting. Default = salt2
        params_to_fit (list): List of parameters to fit
    """

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # To store fit results
    names = ['cid', 'class', 'fit_z', 'z']
    names.extend(params_to_fit)
    names.append('z_err')
    names.extend((v + '_err' for v in params_to_fit))
    names.extend(('chi', 'dof', 'message'))
    out_table = Table(names=names, dtype=[object for _ in names])
    out_path = os.path.join(out_dir, 'summary.csv')

    for input_table in iter_sncosmo_input(bands=bands):
        z_was_fit = int(input_table.meta['redshift'] == -9)
        new_row = [input_table.meta['cid'], input_table.meta['classification'], z_was_fit]

        try:
            result = run_fit_for_object(input_table, model_name, params_to_fit)

        except (sncosmo.fitting.DataQualityError, RuntimeError, ValueError) as e:
            mask = np.full(len(out_table.colnames) - 4, np.NAN).tolist()
            new_row.extend(mask)
            new_row.append(e)

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
    print('Fitting type Ia model in ugriz')
    fit_sdss_data('./sncosmo_results/snia_all_bands')

    print('Fitting type Ia model in ug')
    fit_sdss_data('./sncosmo_results/snia_ug', bands=['sdssu', 'sdssg'])

    print('Fitting type Ia model in riz')
    fit_sdss_data('./sncosmo_results/snia_riz', bands=['sdssr', 'sdssi', 'sdssz'])

    # 91bg model
    print('Fitting 91bg model in ugriz')
    fit_sdss_data('./sncosmo_results/91bg_all_bands',
                  model_name='nugent-sn91bg',
                  params_to_fit=['t0', 'amplitude'])

    print('Fitting 91bg model in ug')
    fit_sdss_data('./sncosmo_results/91bg_ug',
                  model_name='nugent-sn91bg',
                  bands=['sdssu', 'sdssg'],
                  params_to_fit=['t0', 'amplitude'])

    print('Fitting 91bg model in riz')
    fit_sdss_data('./sncosmo_results/91bg_riz',
                  model_name='nugent-sn91bg',
                  bands=['sdssr', 'sdssi', 'sdssz'],
                  params_to_fit=['t0', 'amplitude'])


