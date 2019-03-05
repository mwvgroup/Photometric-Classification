#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This script fits DES light curves using sncosmo"""

import os
import sys

import numpy as np
import sncosmo
from sncosmo.fitting import DataQualityError

sys.path.insert(0, '../')
from data_access import des
from _utils import create_empty_summary_table, count_points_per_band


def fit_des_data(out_path, model, rest_bands=None, **kwargs):
    """Fit DES light curves with SNCosmo

    Args:
        out_path      (str): Where to write fit results
        model       (model): Model to use for fitting
        rest_bands   (list): Optional list of rest frame band-passes to fit

        Additionally any arguments for sncosmo.fit_lc
    """

    out_dir = os.path.dirname(out_path)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Run fit for each target
    out_table = create_empty_summary_table(des.band_names)
    for input_table in des.iter_sncosmo_input(bands=rest_bands, verbose=True):

        # Create a new, incomplete row for the table
        new_row = [input_table.meta['cid']]
        band_count = count_points_per_band(input_table['band'], des.band_names)
        new_row.extend(band_count)

        try:
            model.set(z=input_table.meta['redshift'])
            result, fitted_model = sncosmo.fit_lc(
                data=input_table,
                model=model,
                vparam_names=['t0', 'x0', 'x1', 'c'],
                **kwargs)

        except (DataQualityError, RuntimeError, ValueError) as e:
            new_row.append(input_table.meta['redshift'])
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
    sncosmo_args = dict(bounds=None,
                        modelcov=True,
                        minsnr=5,
                        warn=False)

    print('Fitting type Ia model in all bands', flush=True)
    fit_des_data('./des_results/snia_ugriz.csv',
                 model=salt_2_4,
                 **sncosmo_args)

    print('\n\nFitting type Ia model in ug', flush=True)
    fit_des_data('./des_results/snia_ug.csv',
                 model=salt_2_4,
                 rest_bands=['desu', 'desg'],
                 **sncosmo_args)

    print('\n\nFitting type Ia model in riz', flush=True)
    fit_des_data('./des_results/snia_riz.csv',
                 model=salt_2_4,
                 rest_bands=['desu', 'desg'],
                 **sncosmo_args)
