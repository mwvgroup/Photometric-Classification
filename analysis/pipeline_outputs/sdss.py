#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This script fits SDSS light curves using sncosmo"""

import sys
from itertools import product

import sncosmo

sys.path.insert(0, '../')
from analysis_pipeline.lc_fitting.sdss import fit_data

if __name__ == '__main__':
    salt_2_0 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.0'))
    salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))
    nugent_91bg = sncosmo.Model(source=sncosmo.get_source('nugent-sn91bg'))

    classes_to_skip = ['AGN', 'SLSN', 'SNII', 'Variable', 'pSNII', 'zSNII']
    all_bands = [f'91bg_proj_sdss_{b}{c}' for b, c in
                 product('ugriz', '123456')]
    blue_bands = all_bands[:12]
    red_bands = all_bands[12:]

    sncosmo_args = dict(bounds=None,
                        modelcov=True,
                        phase_range=[-15, 45],
                        minsnr=5,
                        warn=False)

    print('Fitting Salt 2.0 in all bands', flush=True)
    fit_data('./sdss_salt20_all.csv',
             model=salt_2_0,
             keep_types=['zSNIa', 'pSNIa', 'SNIa', 'SNIa?'],
             **sncosmo_args)

    print('\n\nFitting Salt 2.4 in all bands', flush=True)
    fit_data('./sdss_salt24_all.csv',
             model=salt_2_4,
             skip_types=classes_to_skip,
             **sncosmo_args)

    print(f'\n\nFitting Salt 2.4 in {blue_bands}', flush=True)
    fit_data('./sdss_salt24_blue.csv',
             model=salt_2_4,
             rest_bands=blue_bands,
             skip_types=classes_to_skip,
             **sncosmo_args)

    print(f'\n\nFitting Salt 2.4 in {red_bands}', flush=True)
    fit_data('./sdss_salt24_red.csv',
             model=salt_2_4,
             rest_bands=red_bands,
             skip_types=classes_to_skip,
             **sncosmo_args)
