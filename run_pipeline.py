#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Run the analysis pipeline for all data (CSP, DES, and SDSS) and all models
(SALT 2.4 and 91bg).
"""

import os

import tqdm

from analysis_pipeline.lc_fitting import LCFitting

file_dir = os.path.dirname(os.path.abspath(__file__))
lc_fitting = LCFitting('./fitting_params.yml')

csp_dir = os.path.join(file_dir, 'pipeline_outputs/csp')
des_dir = os.path.join(file_dir, 'pipeline_outputs/des')
sdss_dir = os.path.join(file_dir, 'pipeline_outputs/sdss')

for path in (csp_dir, des_dir, sdss_dir):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)

tqdm.tqdm.write('Fitting SDSS for comparison with published values.')
classes_to_skip = ['AGN', 'SLSN', 'SNII', 'Variable']
lc_fitting.fit_sdss(sdss_dir, models=['salt_2_0'],
                    num_params=[4],
                    bands=['all'],
                    skip_types=classes_to_skip)

tqdm.tqdm.write('Fitting CSP')
lc_fitting.fit_csp(csp_dir,
                   models=['salt_2_4', 'sn_91bg'],
                   num_params=[4, 5],
                   bands=['all', 'blue', 'red'])

tqdm.tqdm.write('Fitting DES')
lc_fitting.fit_des(des_dir,
                   models=['salt_2_4', 'sn_91bg'],
                   num_params=[4, 5],
                   nest=False,
                   bands=['all', 'blue', 'red'])

tqdm.tqdm.write('Fitting SDSS')
lc_fitting.fit_sdss(sdss_dir,
                    models=['salt_2_4', 'sn_91bg'],
                    num_params=[4, 5],
                    nest=False,
                    bands=['all', 'blue', 'red'],
                    skip_types=classes_to_skip)
