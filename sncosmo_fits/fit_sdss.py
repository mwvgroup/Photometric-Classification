#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This script fits SDSS light curves using sncosmo"""

import os
import sys
from itertools import product

import numpy as np
import sncosmo
from astropy.table import Table
from sncosmo.fitting import DataQualityError

sys.path.insert(0, '../')
from data_access import sdss


def create_empty_summary_table(params_to_fit):
    """Returns a table with columns:

         cid, class, z, *params_to_fit, z_err, *params_to_fit + _err,
         chi, dof, message.

    Args:
        params_to_fit (list): List of fit parameters
    """

    names = ['cid', 'class', 'z']
    names.extend(params_to_fit)
    names.append('z_err')
    names.extend((v + '_err' for v in params_to_fit))
    names.extend(('chi', 'dof', 'message'))
    out_table = Table(names=names, dtype=[object for _ in names])

    return out_table


def fit_sdss_data(out_path, model, vparam_names, bands=None, fit_types=(),
                  **kwargs):
    """Fit SDSS light curves with SNCosmo

    Args:
        out_path      (str): Where to write fit results
        model       (model): Model to use for fitting
        vparam_names (list): List of parameters to fit
        bands        (list): Optional list of rest frame band-passes to fit
        fit_types    (list): Optional include only certain SDSS II target
                               classifications (case sensitive)

        Additionally any arguments for sncosmo.fit_lc
    """

    out_dir = os.path.dirname(out_path)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Run fit for each target
    out_table = create_empty_summary_table(vparam_names)
    for input_table in sdss.iter_sncosmo_input(
            bands=bands, keep_types=fit_types, verbose=True):

        # Only fit target with published redshift
        z = input_table.meta['redshift']
        if z < 0:
            continue

        # Create a new, incomplete row for the table
        new_row = [input_table.meta['cid'], input_table.meta['classification']]

        try:
            model.set(z=z)
            result, fitted_model = sncosmo.fit_lc(
                input_table, model, ['t0', 'x0', 'x1', 'c'], **kwargs)

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
    salt_2_0 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.0'))
    salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))
    blue_bands = [f'91bg_proj_sdss_{b}{c}' for b, c in product('ug', '123456')]
    red_bands = [f'91bg_proj_sdss_{b}{c}' for b, c in product('riz', '123456')]

    print('Fitting type Ia targets in all bands (Salt 2.0)')
    fit_sdss_data('./sdss_results/snia_ugriz.csv',
                  model=salt_2_0,
                  vparam_names=['t0', 'x0', 'x1', 'c'],
                  fit_types=['zSNIa', 'pSNIa', 'SNIa', 'SNIa?'],
                  bounds=None,
                  modelcov=True,
                  phase_range=[-15, 45],
                  minsnr=5,
                  warn=False)

    print('\nFitting all targets in all bands (Salt 2.4)')
    fit_sdss_data('./sdss_results/all_ugriz.csv',
                  model=salt_2_4,
                  vparam_names=['t0', 'x0', 'x1', 'c'],
                  bounds=None,
                  modelcov=True,
                  phase_range=[-15, 45],
                  minsnr=5,
                  warn=False)

    print('\nFitting all targets in ug (Salt 2.4)')
    fit_sdss_data('./sdss_results/all_ug.csv',
                  model=salt_2_4,
                  vparam_names=['t0', 'x0', 'x1', 'c'],
                  bands=blue_bands,
                  bounds=None,
                  modelcov=True,
                  phase_range=[-15, 45],
                  minsnr=5,
                  warn=False)

    print('\nFitting all targets in riz (Salt 2.4)')
    fit_sdss_data('./sdss_results/all_riz.csv',
                  model=salt_2_4,
                  vparam_names=['t0', 'x0', 'x1', 'c'],
                  bands=red_bands,
                  bounds=None,
                  modelcov=True,
                  phase_range=[-15, 45],
                  minsnr=5,
                  warn=False)
