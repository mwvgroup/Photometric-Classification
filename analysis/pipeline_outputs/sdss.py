#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This script fits SDSS light curves using sncosmo"""

import sys
from itertools import product

import sncosmo
from tqdm import tqdm

sys.path.insert(0, '../')
from analysis_pipeline.lc_fitting import fit_5_param, fit_4_param
from analysis_pipeline.data_access import sdss

salt_2_0 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.0'))
salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))

all_bands = [f'91bg_proj_sdss_{b}{c}' for b, c in product('ugriz', '123456')]
blue_bands = all_bands[:12]
red_bands = all_bands[12:]

validation_classes = ['zSNIa', 'pSNIa', 'SNIa', 'SNIa?']
classes_to_skip = ['AGN', 'SLSN', 'SNII', 'Variable', 'pSNII', 'zSNII']

validation_data = sdss.iter_sncosmo_input(keep_types=validation_classes)
all_band_data = sdss.iter_sncosmo_input(skip_types=classes_to_skip)
blue_band_data = sdss.iter_sncosmo_input(skip_types=classes_to_skip, bands=blue_bands)
red_band_data = sdss.iter_sncosmo_input(skip_types=classes_to_skip, bands=red_bands)

sncosmo_args = dict(bounds={'t0': (53600, 54500),
                            'x0': (0, 0.015),
                            'x1': (-5, 5),
                            'c': (-.5, 1)},

                    modelcov=True,
                    minsnr=5,
                    warn=False)

tqdm.write('Fitting Salt 2.0 - 4 param in all bands')
fit_4_param('./sdss_salt20_4param_all.csv',
            validation_data,
            sdss.band_names,
            model=salt_2_0,
            **sncosmo_args)

tqdm.write('\n\nFitting Salt 2.4 - 4 param in all bands')
fit_4_param('./sdss_salt24_4param_all.csv',
            all_band_data,
            sdss.band_names,
            model=salt_2_4,
            **sncosmo_args)

tqdm.write(f'\n\nFitting Salt 2.4 - 4 param in {blue_bands}')
fit_4_param('./sdss_salt24_4param_blue.csv',
            blue_band_data,
            sdss.band_names,
            model=salt_2_4,
            **sncosmo_args)

tqdm.write(f'\n\nFitting Salt 2.4 - 4 param in {red_bands}')
fit_4_param('./sdss_salt24_4param_red.csv',
            red_band_data,
            sdss.band_names,
            model=salt_2_4,
            **sncosmo_args)

# Update sncosmo params for 5 parameter fit
sncosmo_args['bounds']['z'] = (0.00001, 6.5)

tqdm.write('\n\nFitting Salt 2.4 - 5 param in all bands')
fit_5_param('./sdss_salt24_5param_all.csv',
            all_band_data,
            sdss.band_names,
            model=salt_2_4,
            **sncosmo_args)

tqdm.write(f'\n\nFitting Salt 2.4 - 5 param in {blue_bands}')
fit_5_param('./sdss_salt24_5param_blue.csv',
            blue_band_data,
            sdss.band_names,
            model=salt_2_4,
            **sncosmo_args)

tqdm.write(f'\n\nFitting Salt 2.4 - 5 param in {red_bands}')
fit_5_param('./sdss_salt24_5param_red.csv',
            red_band_data,
            sdss.band_names,
            model=salt_2_4,
            **sncosmo_args)
