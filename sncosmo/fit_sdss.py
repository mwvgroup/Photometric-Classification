#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module fits SDSS light curves using sncosmo"""

import os

import numpy as np
import sncosmo
from astropy.table import Table
from sncosmo.fitting import DataQualityError

import sys; sys.path.insert(0, '../')
from data_access.sdss_data import iter_sncosmo_input

SDSS_BANDS = ('sdssu', 'sdssg', 'sdssr', 'sdssi', 'sdssz')


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


def create_empty_summary_table(bands, params_to_fit):
    """Returns a table with columns:

         cid, class, num_points_ + *bands, fit_z, z,
         *params_to_fit, z_err, *params_to_fit + _err,
         chi, dof, message

    Args:
        bands         (list): List of SDSS bandpasses sdss<ugriz>
        params_to_fit (list): List of fit parameters
    """

    names = ['cid', 'class']
    names.extend(['num_points_' + band for band in bands])
    names.append('fit_z')
    names.append('z')
    names.extend(params_to_fit)
    names.append('z_err')
    names.extend((v + '_err' for v in params_to_fit))
    names.extend(('chi', 'dof', 'message'))
    out_table = Table(names=names, dtype=[object for _ in names])

    return out_table


def fit_sdss_data(out_path,
                  model_name='salt2',
                  bands=SDSS_BANDS,
                  params_to_fit=('t0', 'x0', 'x1', 'c'),
                  skip_types=()):
    """Fit SDSS light curves with SNCosmo

    Files are named as <out_dir>/<target cid>.txt

    Args:
        out_path       (str): Where to write fit results
        model_name     (str): Model to use for fitting. Default = salt2
        params_to_fit (list): List of parameters to fit
        skip_types    (list): List of case sensitive classifications to skip
        bands         (list): Optional list of bandpasses to fit
    """

    out_dir = os.path.dirname(out_path)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Run fit for each target
    out_table = create_empty_summary_table(bands, params_to_fit)
    for input_table in iter_sncosmo_input(
            bands=bands, skip_types=skip_types, verbose=True):

        # Determine if redshift is fit or given
        z = input_table.meta['redshift']
        z_was_fit = int(z == -9)

        # Determine number of data points per band
        band_names, band_counts = np.unique(input_table['band'], return_counts=True)
        count_dict = dict(zip(band_names, band_counts))
        num_data_points = [count_dict.get(band, 0) for band in bands]

        # Create a new, incomplete row for the table
        new_row = [input_table.meta['obj_id'], input_table.meta['classification']]
        new_row.extend(num_data_points)
        new_row.append(z_was_fit)

        try:
            result = run_fit_for_object(input_table, model_name, params_to_fit)

        except (DataQualityError, RuntimeError, ValueError) as e:
            mask_length = len(out_table.colnames) - len(new_row) - 2
            mask = np.full(mask_length, np.NAN).tolist()
            new_row.append(z)
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
    print('\n\nFitting type Ia model in ug')
    fit_sdss_data('./sncosmo_results/snia_ug.csv',
                  skip_types=['Variable', 'AGN'],
                  bands=['sdssu', 'sdssg'])

    print('\n\nFitting type Ia model in riz')
    fit_sdss_data('./sncosmo_results/snia_riz.csv',
                  skip_types=['Variable', 'AGN'],
                  bands=['sdssr', 'sdssi', 'sdssz'])

    print('\n\nFitting 91bg model in ug')
    fit_sdss_data('./sncosmo_results/91bg_ug.csv',
                  skip_types=['Variable', 'AGN'],
                  model_name='nugent-sn91bg',
                  bands=['sdssu', 'sdssg'],
                  params_to_fit=['t0', 'amplitude'])

    print('\n\nFitting 91bg model in riz')
    fit_sdss_data('./sncosmo_results/91bg_riz.csv',
                  skip_types=['Variable', 'AGN'],
                  model_name='nugent-sn91bg',
                  bands=['sdssr', 'sdssi', 'sdssz'],
                  params_to_fit=['t0', 'amplitude'])

    print('Fitting type Ia model in all bands')
    fit_sdss_data('./sncosmo_results/snia_ugriz.csv',
                  skip_types=['Variable', 'AGN'])
