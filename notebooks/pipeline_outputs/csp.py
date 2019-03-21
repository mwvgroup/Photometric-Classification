#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This script fits CSP light curves using sncosmo"""

import os
import sys

import sncosmo
from tqdm import tqdm

sys.path.insert(0, '../')
from analysis_pipeline.lc_fitting import fit_n_params
from analysis_pipeline.data_access import csp
from sncosmo_91bgmodel import SN91bgSource

# Get models for fitting
salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))
sn_91bg = sncosmo.Model(source=SN91bgSource())


def fit_csp(out_dir):
    # Define red and blue bandpasses
    blue_bands = ['u', 'g', 'B', 'V0', 'V', 'Y']
    blue_bands = [f'91bg_proj_csp_{f}' for f in blue_bands]
    red_bands = ['r', 'i', 'H', 'J', 'Jrc2', 'Ydw', 'Jdw', 'Hdw']
    red_bands = [f'91bg_proj_csp_{f}' for f in red_bands]

    # Define model specific arguments for fitting
    salt24_4param = dict(bounds=None, modelcov=True)
    salt24_5param = dict(bounds={'z': (0.002, 0.085)}, modelcov=True)
    sn91bg_4param = dict(bounds={'x1': (0.65, 1.25), 'c': (0, 1)})
    sn91bg_5param = sn91bg_4param.copy()
    sn91bg_5param['z'] = (0.002, 0.085)

    modeling_data = zip(
        (salt24_4param, salt24_5param, sn91bg_4param, sn91bg_5param),
        (salt_2_4, salt_2_4, sn_91bg, sn_91bg)
    )

    path_pattern = '{}_{}_{}param_{}.csv'
    for model_args, model in modeling_data:
        for num_param in (4, 5):
            for bands in (None, blue_bands, red_bands):
                bands_str = bands if bands else 'all'
                model_name = model.source.name
                tqdm.write(f'Fitting {model_name} - {num_param} params in {bands_str}')

                fname = path_pattern.format('csp', model_name, num_param, bands_str)
                out_path = os.path.join(out_dir, fname)
                fit_n_params(out_path,
                             num_params=num_param,
                             inputs=csp.iter_sncosmo_input(bands),
                             bands=csp.band_names,
                             model=model,
                             **model_args)


if __name__ == '__main__':
    fit_csp('./')
