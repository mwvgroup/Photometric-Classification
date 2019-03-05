#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This script fits SDSS light curves using sncosmo"""

import os
import sys
from itertools import product

import numpy as np
import sncosmo
from sncosmo.fitting import DataQualityError

sys.path.insert(0, '../')
from data_access import sdss
from _utils import create_empty_summary_table, count_points_per_band


def fit_sdss_data(out_path, model, rest_bands=None, fit_types=(), **kwargs):
    """Fit SDSS light curves with SNCosmo

    Args:
        out_path      (str): Where to write fit results
        model       (model): Model to use for fitting
        rest_bands   (list): Optional list of rest frame band-passes to fit
        fit_types    (list): Optional include only certain SDSS II target
                               classifications (case sensitive)

        Additionally any arguments for sncosmo.fit_lc
    """

    out_dir = os.path.dirname(out_path)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Run fit for each target
    out_table = create_empty_summary_table(sdss.band_names, ['class'])
    for input_table in sdss.iter_sncosmo_input(
            bands=rest_bands, keep_types=fit_types, verbose=True):

        # Only fit target with published redshift
        z = input_table.meta['redshift']
        if z < 0:
            continue

        # Create a new, incomplete row for the table
        new_row = [input_table.meta['classification'], input_table.meta['cid']]
        band_count = count_points_per_band(input_table['band'], sdss.band_names)
        new_row.extend(band_count)

        try:
            model.set(z=z)
            result, fitted_model = sncosmo.fit_lc(
                data=input_table,
                model=model,
                vparam_names=['t0', 'x0', 'x1', 'c'],
                **kwargs)

        except (DataQualityError, RuntimeError, ValueError) as e:
            new_row.append(z)
            new_row.extend(np.full(4, np.NAN).tolist())
            new_row.append(input_table.meta['redshift'])
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
    salt_2_0 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.0'))
    salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))

    classifications_to_fit = ['zSNIa', 'pSNIa', 'SNIa', 'SNIa?']
    all_bands = [f'91bg_proj_sdss_{b}{c}' for b, c in product('ugriz', '123456')]
    blue_bands = all_bands[:12]
    red_bands = all_bands[12:]

    sncosmo_args = dict(bounds=None,
                        modelcov=True,
                        phase_range=[-15, 45],
                        minsnr=5,
                        warn=False)

    print('Fitting type Ia targets in all bands (Salt 2.0)', flush=True)
    fit_sdss_data('./sdss_results/snia_ugriz_2_0.csv',
                  model=salt_2_0,
                  fit_types=classifications_to_fit,
                  **sncosmo_args)

    print('\n\nFitting all targets in all bands (Salt 2.4)', flush=True)
    fit_sdss_data('./sdss_results/snia_ugriz.csv',
                  model=salt_2_4,
                  fit_types=classifications_to_fit,
                  **sncosmo_args)

    print('\n\nFitting all targets in ug (Salt 2.4)', flush=True)
    fit_sdss_data('./sdss_results/snia_ug.csv',
                  model=salt_2_4,
                  rest_bands=blue_bands,
                  **sncosmo_args)

    print('\n\nFitting all targets in riz (Salt 2.4)', flush=True)
    fit_sdss_data('./sdss_results/snia_riz.csv',
                  model=salt_2_4,
                  rest_bands=red_bands,
                  fit_types=classifications_to_fit,
                  **sncosmo_args)
