#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module provides functions for fitting CSP, DES, and SDSS data with a
SALT2 and 91bg model in all bands, the red bands, and the blue bands.
"""

import os
import time
from itertools import product

import sncosmo
from tqdm import tqdm

from . import fit_n_params
from ..data_access import csp, sdss, des
from ..sn91bg_model import SN91bgSource

# Get models for fitting
salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))
sn_91bg = sncosmo.Model(source=SN91bgSource())


def _run_fit_set(iter_inputs_func, out_dir, blue_bands, red_bands, modeling_data):
    """Run a four and five parameter fit using the Salt2 and 91bg model in
        all bands, the rest frame blue, and the rest frame red.

        Args:
            iter_inputs_func                  (func): Function returning an iterable of SNCosmo input tables
            out_dir                            (str): Directory where results are saved
            blue_bands                   (list[str]): List of blue band-passes
            red_bands                    (list[str]): List of blue band-passes
            modeling_data (iter[tuple[dict, Model]]): Iterable of sncosmo models and their arguments
        """

    path_pattern = '{}_{}param_{}.csv'
    band_data = zip(('all', 'blue', 'red'), (None, blue_bands, red_bands))
    for model_args, model in modeling_data:
        for num_param in (4, 5):
            for bands_str, bands in band_data:
                model_name = model.source.name
                tqdm.write(f'Fitting {num_param} params in {bands_str}')
                time.sleep(0.5)  # Give output time to flush

                fname = path_pattern.format(model_name, num_param, bands_str)
                out_path = os.path.join(out_dir, fname)
                fit_n_params(out_path,
                             num_params=num_param,
                             inputs=iter_inputs_func(bands),
                             bands=csp.band_names,
                             model=model,
                             **model_args)

                tqdm.write('\n')


def fit_csp(out_dir):
    """Fit CSP data and save result to file

    Args:
        out_dir (str): Directory where results are saved
    """

    # Define red and blue bandpasses
    blue_bands = ['u', 'g', 'B', 'V0', 'V', 'Y']
    blue_bands = [f'91bg_proj_csp_{f}' for f in blue_bands]
    red_bands = ['r', 'i', 'H', 'J', 'Jrc2', 'Ydw', 'Jdw', 'Hdw']
    red_bands = [f'91bg_proj_csp_{f}' for f in red_bands]

    # Define model specific arguments for fitting
    salt24_4param = dict(bounds=None, modelcov=True)
    sn91bg_4param = dict(bounds={'x1': (0.65, 1.25), 'c': (0, 1)})

    salt24_5param = salt24_4param.copy()
    salt24_5param['bounds'] = {'z': (0.002, 0.085)}

    sn91bg_5param = sn91bg_4param.copy()
    sn91bg_5param['bounds']['z'] = (0.002, 0.085)

    modeling_data = zip(
        (salt24_4param, salt24_5param, sn91bg_4param, sn91bg_5param),
        (salt_2_4, salt_2_4, sn_91bg, sn_91bg)
    )

    _run_fit_set(csp.iter_sncosmo_input,
                 out_dir,
                 blue_bands,
                 red_bands,
                 modeling_data)


def fit_des(out_dir):
    """Fit DES data and save result to file

    Args:
        out_dir (str): Directory where results are saved
    """

    # Define red and blue bandpasses
    blue_bands = ['desg', 'desr']
    red_bands = ['desi', 'desz', 'desy']

    # Define model specific arguments for fitting
    salt24_4param = dict(bounds={'t0': (51900, 57420),
                                 'x0': (0, 0.05),
                                 'x1': (-5, 5),
                                 'c': (-.5, 1)},
                         modelcov=True)

    sn91bg_4param = dict(bounds={'x1': (0.65, 1.25), 'c': (0, 1)})

    salt24_5param = salt24_4param.copy()
    salt24_5param['bounds']['z'] = (0.01, 0.9)

    sn91bg_5param = sn91bg_4param.copy()
    sn91bg_5param['bounds']['z'] = (0.01, 0.9)

    modeling_data = zip(
        (salt24_4param, salt24_5param, sn91bg_4param, sn91bg_5param),
        (salt_2_4, salt_2_4, sn_91bg, sn_91bg)
    )

    _run_fit_set(des.iter_sncosmo_input,
                 out_dir,
                 blue_bands,
                 red_bands,
                 modeling_data)


def fit_sdss(out_dir):
    """Fit SDSS data and save result to file

    Args:
        out_dir (str): Directory where results are saved
    """

    # Define red and blue bandpasses
    blue_bands = [f'91bg_proj_sdss_{b}{c}' for b, c in product('ug', '123456')]
    red_bands = [f'91bg_proj_sdss_{b}{c}' for b, c in product('riz', '123456')]

    # Define model specific arguments for fitting
    salt24_4param = dict(bounds={'t0': (53600, 54500),
                                 'x0': (0, 0.015),
                                 'x1': (-5, 5),
                                 'c': (-.5, 1)},
                         modelcov=True)

    sn91bg_4param = dict(bounds={'x1': (0.65, 1.25), 'c': (0, 1)})

    salt24_5param = salt24_4param.copy()
    salt24_5param['bounds']['z'] = (0.00001, 6.5)

    sn91bg_5param = sn91bg_4param.copy()
    sn91bg_5param['bounds']['z'] = (0.00001, 6.5)

    modeling_data = zip(
        (salt24_4param, salt24_5param, sn91bg_4param, sn91bg_5param),
        (salt_2_4, salt_2_4, sn_91bg, sn_91bg)
    )

    _run_fit_set(sdss.iter_sncosmo_input,
                 out_dir,
                 blue_bands,
                 red_bands,
                 modeling_data)
