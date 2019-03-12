#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This script fits CSP light curves using sncosmo"""

import sys

import sncosmo

sys.path.insert(0, '../')
from analysis_pipeline.lc_fitting import fit_5_param, fit_4_param
from analysis_pipeline.data_access import csp

# Get models for fitting
salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))

# Define red and blue bandpasses
blue_bands = ['u', 'g', 'B', 'V0', 'V', 'Y']
blue_bands = [f'91bg_proj_csp_{f}' for f in blue_bands]
red_bands = ['r', 'i', 'H', 'J', 'Jrc2', 'Ydw', 'Jdw', 'Hdw']
red_bands = [f'91bg_proj_csp_{f}' for f in red_bands]

# Define arguments for SNCosmo
sncosmo_args = dict(bounds={'z': (0.1, 0.8)},
                    modelcov=True,
                    minsnr=5,
                    warn=False)

# Run Fitting
print('Fitting Salt 2.4 - 4 param in all bands', flush=True)
fit_4_param('./csp_salt24_4_all.csv',
            csp.iter_sncosmo_input(),
            csp.band_names,
            model=salt_2_4,
            **sncosmo_args)

print(f'\n\nFitting Salt 2.4 - 4 param in {blue_bands}', flush=True)
fit_4_param('./csp_salt24_4_blue.csv',
            csp.iter_sncosmo_input(bands=blue_bands),
            csp.band_names,
            model=salt_2_4,
            **sncosmo_args)

print(f'\n\nFitting Salt 2.4 - 4 param in {red_bands}', flush=True)
fit_4_param('./csp_salt24_4_red.csv',
            csp.iter_sncosmo_input(bands=red_bands),
            csp.band_names,
            model=salt_2_4,
            **sncosmo_args)

print('Fitting Salt 2.4 - 5 param in all bands', flush=True)
fit_5_param('./csp_salt24_5_all.csv',
            csp.iter_sncosmo_input(),
            csp.band_names,
            model=salt_2_4,
            **sncosmo_args)

print(f'\n\nFitting Salt 2.4 - 5 param in {blue_bands}', flush=True)
fit_5_param('./csp_salt24_5_blue.csv',
            csp.iter_sncosmo_input(bands=blue_bands),
            csp.band_names,
            model=salt_2_4,
            **sncosmo_args)

print(f'\n\nFitting Salt 2.4 - 5 param in {red_bands}', flush=True)
fit_5_param('./csp_salt24_5_red.csv',
            csp.iter_sncosmo_input(bands=red_bands),
            csp.band_names,
            model=salt_2_4,
            **sncosmo_args)
