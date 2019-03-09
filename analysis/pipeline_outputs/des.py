#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This script fits DES light curves using sncosmo"""

import sys

import sncosmo

sys.path.insert(0, '../')
from analysis_pipeline.lc_fitting.des import fit_data

if __name__ == '__main__':
    salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))

    blue_bands = ['desg', 'desr']
    red_bands = ['desi', 'desz', 'desy']

    sncosmo_args = dict(bounds=None,
                        modelcov=True,
                        minsnr=5,
                        warn=False)

    print('Fitting Salt 2.4 in all bands', flush=True)
    fit_data('./des_salt24_all.csv',
             model=salt_2_4,
             **sncosmo_args)

    print(f'\n\nFitting Salt 2.4 in {blue_bands}', flush=True)
    fit_data('./des_salt24_blue.csv',
             model=salt_2_4,
             rest_bands=blue_bands,
             **sncosmo_args)

    print(f'\n\nFitting Salt 2.4 in {red_bands}', flush=True)
    fit_data('./des_salt24_red.csv',
             model=salt_2_4,
             rest_bands=red_bands,
             **sncosmo_args)
