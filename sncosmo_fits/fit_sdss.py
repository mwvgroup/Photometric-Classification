#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This script fits SDSS light curves using sncosmo"""

import os

import numpy as np
import sncosmo
from astropy.table import Table
from sncosmo.fitting import DataQualityError

import sys; sys.path.insert(0, '../')
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


def fit_sdss_data(out_path,
                  model,
                  bands=None,
                  params_to_fit=('t0', 'x0', 'x1', 'c'),
                  fit_types=()):
    """Fit SDSS light curves with SNCosmo

    Files are named as <out_dir>/<target cid>.txt

    Args:
        out_path       (str): Where to write fit results
        model        (model): Model to use for fitting. Default = salt2
        params_to_fit (list): List of parameters to fit
        fit_types    (list): List of case sensitive classifications to skip
        bands         (list): Optional list of band-passes to fit
    """

    out_dir = os.path.dirname(out_path)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Run fit for each target
    out_table = create_empty_summary_table(params_to_fit)
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
                input_table, model, params_to_fit, bounds=None)

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
    classifications_to_fit = ['zSNIa', 'pSNIa', 'SNIa', 'SNIa?']

    source = sncosmo.get_source('salt2', version='2.0')
    fitting_model = sncosmo.Model(source=source)
    print('Fitting type Ia model in all bands')
    fit_sdss_data('./sdss_results/snia_ugriz.csv',
                  model=fitting_model,
                  fit_types=classifications_to_fit)
