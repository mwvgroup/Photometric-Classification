#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This script fits DES light curves using sncosmo"""

import os
import sys

import numpy as np
import sncosmo
from astropy.table import Table
from sncosmo.fitting import DataQualityError

sys.path.insert(0, '../')
from data_access import des


def create_empty_summary_table():
    """Returns a table with columns:

         cid, z, t0, x0, x1, z_err, t0_err, x0_err, x1_err, c_err,
         chi, dof, message.
    """

    names = ['cid']
    names.extend(['num_points_' + band for band in des.band_names])

    param_names = ['z', 't0', 'x0', 'x1', 'c']
    names.extend(param_names)
    names.extend((p + '_err' for p in param_names))
    names.extend(('chi', 'dof', 'message'))
    out_table = Table(names=names, dtype=[object for _ in names])

    return out_table


def count_data_points_per_band(band_list):
    # Determine number of data points per band
    band_names, band_counts = np.unique(band_list, return_counts=True)
    count_dict = dict(zip(band_names, band_counts))
    return [count_dict.get(band, 0) for band in des.band_names]


def fit_des_data(out_path, model, bands=None, **kwargs):
    """Fit DES light curves with SNCosmo

    Args:
        out_path      (str): Where to write fit results
        model       (model): Model to use for fitting
        bands        (list): Optional list of rest frame band-passes to fit

        Additionally any arguments for sncosmo.fit_lc
    """

    out_dir = os.path.dirname(out_path)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Run fit for each target
    out_table = create_empty_summary_table()
    for input_table in des.iter_sncosmo_input(bands=bands, verbose=True):

        # Create a new, incomplete row for the table
        new_row = [input_table.meta['cid']]
        new_row.extend(count_data_points_per_band(input_table['band']))

        try:
            model.set(z=input_table.meta['redshift'])
            result, fitted_model = sncosmo.fit_lc(
                data=input_table,
                model=model,
                vparam_names=['t0', 'x0', 'x1', 'c'],
                **kwargs)

        except (DataQualityError, RuntimeError, ValueError) as e:
            new_row.append(input_table.meta['redshift_err'])
            new_row.extend(np.full(4, np.NAN).tolist())
            new_row.append(input_table.meta['redshift_err'])
            new_row.extend(np.full(6, np.NAN).tolist())
            new_row.append(str(e))

        else:
            new_row.extend(result.parameters)
            new_row.append(input_table.meta['redshift_err'])
            new_row.extend(result.errors.values())
            new_row.append(result.chisq)
            new_row.append(result.ndof)
            new_row.append(result.message)

        out_table.add_row(new_row)
        out_table.write(out_path, overwrite=True)


if __name__ == '__main__':
    salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))

    print('Fitting type Ia model in all bands')
    fit_des_data('./des_results/snia_ugriz.csv',
                 model=salt_2_4,
                 bounds=None,
                 modelcov=True,
                 phase_range=[-15, 45],
                 minsnr=5,
                 warn=False)

    print('\n\nFitting type Ia model in ug')
    fit_des_data('./des_results/snia_ug.csv',
                 model=salt_2_4,
                 bands=['desu', 'desg'],
                 bounds=None,
                 modelcov=True,
                 phase_range=[-15, 45],
                 minsnr=5,
                 warn=False)

    print('\n\nFitting type Ia model in riz')
    fit_des_data('./des_results/snia_riz.csv',
                 model=salt_2_4,
                 bands=['desu', 'desg'],
                 bounds=None,
                 modelcov=True,
                 phase_range=[-15, 45],
                 minsnr=5,
                 warn=False)
