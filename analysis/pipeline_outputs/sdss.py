#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This script fits SDSS light curves using sncosmo"""

import sys
from itertools import product

import sncosmo

sys.path.insert(0, '../')
from analysis_pipeline.lc_fitting import fit_5_param, fit_4_param
from analysis_pipeline.data_access import sdss

salt_2_0 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.0'))
salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))

all_bands = [f'91bg_proj_sdss_{b}{c}' for b, c in product('ugriz', '123456')]
blue_bands = all_bands[:12]
red_bands = all_bands[12:]

vlaidation_classes = ['zSNIa', 'pSNIa', 'SNIa', 'SNIa?']
classes_to_skip = ['AGN', 'SLSN', 'SNII', 'Variable', 'pSNII', 'zSNII']

validation_data = sdss.iter_sncosmo_input(keep_types=vlaidation_classes)
all_band_data = sdss.iter_sncosmo_input(skip_types=classes_to_skip)
blue_band_data = sdss.iter_sncosmo_input(skip_types=classes_to_skip, bands=blue_bands)
red_band_data = sdss.iter_sncosmo_input(skip_types=classes_to_skip, bands=red_bands)

sncosmo_args = dict(bounds={'z': (0.1, 0.8)},
                    modelcov=True,
                    phase_range=[-15, 45],
                    minsnr=5,
                    warn=False)

print('Fitting Salt 2.0 - 4 param in all bands', flush=True)
fit_4_param('./sdss_salt20_4_all.csv',
            validation_data,
            sdss.band_names,
            model=salt_2_0,
            **sncosmo_args)

print('\n\nFitting Salt 2.4 - 4 param in all bands', flush=True)
fit_4_param('./sdss_salt24_4_all.csv',
            all_band_data,
            sdss.band_names,
            model=salt_2_4,
            **sncosmo_args)

print(f'\n\nFitting Salt 2.4 - 4 param in {blue_bands}', flush=True)
fit_4_param('./sdss_salt24_4_blue.csv',
            blue_band_data,
            sdss.band_names,
            model=salt_2_4,
            **sncosmo_args)

print(f'\n\nFitting Salt 2.4 - 4 param in {red_bands}', flush=True)
fit_4_param('./sdss_salt24_4_red.csv',
            red_band_data,
            sdss.band_names,
            model=salt_2_4,
            **sncosmo_args)

print('\n\nFitting Salt 2.4 - 5 param in all bands', flush=True)
fit_5_param('./sdss_salt24_4_all.csv',
            all_band_data,
            sdss.band_names,
            model=salt_2_4,
            **sncosmo_args)

print(f'\n\nFitting Salt 2.4 - 5 param in {blue_bands}', flush=True)
fit_5_param('./sdss_salt24_4_blue.csv',
            blue_band_data,
            sdss.band_names,
            model=salt_2_4,
            **sncosmo_args)

print(f'\n\nFitting Salt 2.4 - 5 param in {red_bands}', flush=True)
fit_5_param('./sdss_salt24_4_red.csv',
            red_band_data,
            sdss.band_names,
            model=salt_2_4,
            **sncosmo_args)
