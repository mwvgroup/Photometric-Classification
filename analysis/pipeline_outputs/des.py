#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This script fits DES light curves using sncosmo"""

import sys

import sncosmo
from tqdm import tqdm

sys.path.insert(0, '../')
from analysis_pipeline.lc_fitting import fit_5_param, fit_4_param
from analysis_pipeline.data_access import des

# Get models for fitting
salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))

# Define red and blue bandpasses
blue_bands = ['desg', 'desr']
red_bands = ['desi', 'desz', 'desy']

# Define arguments for SNCosmo
sncosmo_args = dict(bounds=None,
                    modelcov=True,
                    minsnr=5,
                    warn=False)

# Run Fitting
tqdm.write('Fitting Salt 2.4 - 4 param in all bands')
fit_4_param('./des_salt24_4param_all.csv',
            des.iter_sncosmo_input(),
            des.band_names,
            model=salt_2_4,
            **sncosmo_args)

tqdm.write(f'\n\nFitting Salt 2.4 - 4 param in {blue_bands}')
fit_4_param('./des_salt24_4param_blue.csv',
            des.iter_sncosmo_input(bands=blue_bands),
            des.band_names,
            model=salt_2_4,
            **sncosmo_args)

tqdm.write(f'\n\nFitting Salt 2.4 - 4 param in {red_bands}')
fit_4_param('./des_salt24_4param_red.csv',
            des.iter_sncosmo_input(bands=red_bands),
            des.band_names,
            model=salt_2_4,
            **sncosmo_args)

# Update sncosmo params for 5 parameter fit
sncosmo_args['bounds'] = {'z': (0.01, 0.9)}

tqdm.write('Fitting Salt 2.4 - 5 param in all bands')
fit_5_param('./des_salt24_5param_all.csv',
            des.iter_sncosmo_input(),
            des.band_names,
            model=salt_2_4,
            **sncosmo_args)

tqdm.write(f'\n\nFitting Salt 2.4 - 5 param in {blue_bands}')
fit_5_param('./des_salt24_5param_blue.csv',
            des.iter_sncosmo_input(bands=blue_bands),
            des.band_names,
            model=salt_2_4,
            **sncosmo_args)

tqdm.write(f'\n\nFitting Salt 2.4 - 5 param in {red_bands}')
fit_5_param('./des_salt24_5param_red.csv',
            des.iter_sncosmo_input(bands=red_bands),
            des.band_names,
            model=salt_2_4,
            **sncosmo_args)
