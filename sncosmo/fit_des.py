#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This script fits DES light curves using sncosmo"""

import os

import numpy as np
import sncosmo
from astropy.table import Table
from matplotlib import pyplot as plt
from sncosmo.fitting import DataQualityError

import sys; sys.path.insert(0, '../')
from data_access.des_data import iter_sncosmo_input

DES_BANDS = ('desg', 'desr', 'desi', 'desz', 'desy')


def create_empty_summary_table(bands, params_to_fit):
    """Returns a table with columns:

         cid, class, num_points_ + *bands, fit_z, z,
         *params_to_fit, z_err, *params_to_fit + _err,
         chi, dof, message

    Args:
        bands         (list): List of SDSS bandpasses sdss<ugriz>
        params_to_fit (list): List of fit parameters
    """

    names = ['cid', 'z', 'z_err']
    names.extend(['num_points_' + band for band in bands])
    names.extend(params_to_fit)
    names.extend((p + '_err' for p in params_to_fit))
    names.extend(('chi', 'dof', 'message'))
    out_table = Table(names=names, dtype=[object for _ in names])

    return out_table


def fit_des_data(out_path,
                 model_name='salt2',
                 bands=DES_BANDS,
                 params_to_fit=('t0', 'x0', 'x1', 'c')):
    """Fit DES light curves with SNCosmo

    Files are named as <out_dir>/<target cid>.txt

    Args:
        out_path       (str): Where to write fit results
        model_name     (str): Model to use for fitting. Default = salt2
        params_to_fit (list): List of parameters to fit
        bands         (list): Optional list of bandpasses to fit
    """

    out_dir = os.path.dirname(out_path)
    out_fname = os.path.splitext(os.path.basename(out_path))[0]
    fig_dir = os.path.join(out_dir, out_fname + '_figs')
    if not os.path.exists(fig_dir):
        os.makedirs(fig_dir)

    # Run fit for each target
    out_table = create_empty_summary_table(bands, params_to_fit)
    for input_table in iter_sncosmo_input(bands=bands, verbose=True):

        # Determine number of data points per band
        band_names, band_counts = np.unique(input_table['band'],
                                            return_counts=True)
        count_dict = dict(zip(band_names, band_counts))
        num_data_points = [count_dict.get(band, 0) for band in bands]

        # Create a new, incomplete row for the table
        new_row = [input_table.meta['cid'],
                   input_table.meta['redshift'],
                   input_table.meta['redshift_err']]

        new_row.extend(num_data_points)

        try:
            model = sncosmo.Model(source=model_name)
            model.set(z=input_table.meta['redshift'])
            result, fitted_model = sncosmo.fit_lc(
                input_table, model, list(params_to_fit), bounds=None)

            sncosmo.plot_lc(input_table, model=fitted_model,
                            errors=result.errors)

            f_name = '{}.pdf'.format(input_table.meta['cid'])
            fig_path = os.path.join(fig_dir, f_name)
            plt.savefig(fig_path)

        except (DataQualityError, RuntimeError, ValueError) as e:
            mask_length = len(out_table.colnames) - len(new_row) - 1
            mask = np.full(mask_length, np.NAN).tolist()
            new_row.extend(mask)
            new_row.append(str(e))

        else:
            new_row.extend(result.parameters[1:])  # Slice to remove redshift
            new_row.extend(result.errors.values())
            new_row.append(result.chisq)
            new_row.append(result.ndof)
            new_row.append(result.message)

        out_table.add_row(new_row)
        out_table.write(out_path, overwrite=True)


if __name__ == '__main__':
    print('\n\nFitting type Ia model in all bands')
    fit_des_data('./des_results/snia_ugriz.csv')

    print('Fitting type Ia model in ug')
    fit_des_data('./des_results/snia_ug.csv',
                 bands=['desu', 'desg'])

    print('\n\nFitting type Ia model in riz')
    fit_des_data('./des_results/snia_riz.csv',
                 bands=['desr', 'desi', 'desz'])

    print('\n\nFitting 91bg model in ug')
    fit_des_data('./des_results/91bg_ug.csv',
                 model_name='nugent-sn91bg',
                 bands=['desu', 'desg'],
                 params_to_fit=['t0', 'amplitude'])

    print('\n\nFitting 91bg model in riz')
    fit_des_data('./des_results/91bg_riz.csv',
                 model_name='nugent-sn91bg',
                 bands=['desr', 'desi', 'desz'],
                 params_to_fit=['t0', 'amplitude'])
