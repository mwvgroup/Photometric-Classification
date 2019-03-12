#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This script fits DES light curves using sncosmo"""

import sys

import sncosmo

sys.path.insert(0, '../')
from analysis_pipeline.lc_fitting import fit_5_param, fit_4_param
from analysis_pipeline.data_access import des

# Get models for fitting
salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))

# Define red and blue bandpasses
blue_bands = ['91bg_proj_des_g', '91bg_proj_des_r']
red_bands = ['91bg_proj_des_i', '91bg_proj_des_z', '91bg_proj_des_y']

# Define arguments for SNCosmo
sncosmo_args = dict(bounds=None,
                    modelcov=True,
                    minsnr=5,
                    warn=False)

# Run Fitting
print('Fitting Salt 2.4 - 4 param in all bands', flush=True)
fit_4_param('./des_salt24_4_all.csv',
            des.iter_sncosmo_input(),
            des.band_names,
            model=salt_2_4,
            **sncosmo_args)

print(f'\n\nFitting Salt 2.4 - 4 param in {blue_bands}', flush=True)
fit_4_param('./des_salt24_4_blue.csv',
            des.iter_sncosmo_input(bands=blue_bands),
            des.band_names,
            model=salt_2_4,
            **sncosmo_args)

print(f'\n\nFitting Salt 2.4 - 4 param in {red_bands}', flush=True)
fit_4_param('./des_salt24_4_red.csv',
            des.iter_sncosmo_input(bands=red_bands),
            des.band_names,
            model=salt_2_4,
            **sncosmo_args)

print('Fitting Salt 2.4 - 5 param in all bands', flush=True)
fit_5_param('./des_salt24_5_all.csv',
            des.iter_sncosmo_input(),
            des.band_names,
            model=salt_2_4,
            **sncosmo_args)

print(f'\n\nFitting Salt 2.4 - 5 param in {blue_bands}', flush=True)
fit_5_param('./des_salt24_5_blue.csv',
            des.iter_sncosmo_input(bands=blue_bands),
            des.band_names,
            model=salt_2_4,
            **sncosmo_args)

print(f'\n\nFitting Salt 2.4 - 5 param in {red_bands}', flush=True)
fit_5_param('./des_salt24_5_red.csv',
            des.iter_sncosmo_input(bands=red_bands),
            des.band_names,
            model=salt_2_4,
            **sncosmo_args)
