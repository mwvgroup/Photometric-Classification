#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Run the analysis pipeline for all data (CSP, DES, and SDSS) and all models
(SALT 2.4 and 91bg).
"""

import os
import tqdm

from analysis_pipeline.lc_fitting import fit_csp, fit_des, fit_sdss

file_dir = os.path.dirname(os.path.abspath(__file__))
csp_dir = os.path.join(file_dir, 'pipeline_outputs/csp')
des_dir = os.path.join(file_dir, 'pipeline_outputs/des')
sdss_dir = os.path.join(file_dir, 'pipeline_outputs/sdss')

for path in (csp_dir, des_dir, sdss_dir):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)

tqdm.tqdm.write('Fitting CSP')
fit_csp(csp_dir)

tqdm.tqdm.write('Fitting DES')
fit_des(des_dir)

tqdm.tqdm.write('Fitting SDSS')
fit_sdss(sdss_dir)
