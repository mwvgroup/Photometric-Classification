#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module provides functions for fitting CSP, DES, and SDSS data with a
SALT2 and 91bg model in all bands, the red bands, and the blue bands.
"""

import os
import time
from itertools import product

import sncosmo
import yaml
from tqdm import tqdm

from . import fit_n_params
from ..data_access import csp, des, sdss
from ..sn91bg_model import SN91bgSource

# Get models for fitting
salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))
sn_91bg = sncosmo.Model(source=SN91bgSource())


def _run_fit_set(iter_inputs_func, out_dir, blue_bands, red_bands, params):
    """Run a four and five parameter fit using the Salt2 and 91bg model in
        all bands, the rest frame blue, and the rest frame red.

        Args:
            iter_inputs_func (func): Function returning iterable of SNCosmo input tables
            out_dir           (str): Directory where results are saved
            blue_bands  (list[str]): List of blue band-passes
            red_bands   (list[str]): List of blue band-passes
            params           (dict): SNCosmo fitting arguments per fit
        """

    modeling_data = zip(
        (4, 5, 4, 5),  # Number of fitting params
        (salt_2_4, salt_2_4, sn_91bg, sn_91bg),  # Models
        (params['salt24_4param'], params['salt24_5param'],  # Fitting arguments
         params['sn91bg_4param'], params['sn91bg_5param'])
    )

    path_pattern = '{}_{}param_{}.csv'
    for num_param, model, model_args in modeling_data:
        band_data = zip(('all', 'blue', 'red'), (None, blue_bands, red_bands))
        for bands_str, bands in band_data:
            model_name = model.source.name
            tqdm.write(f'{num_param} param {model_name} in {bands_str} bands:')
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


class LCFitting():

    def __init__(self, params_file=None):
        self.fitting_params = None
        if params_file is not None:
            self.configure_params(params_file)

    def configure_params(self, params_file):
        """Load a yaml configuration file with fit parameters

        Args:
            params_file (Stf): Path of the file to load
        """

        with open(params_file) as ofile:
            self.fitting_params = yaml.load(ofile)

    def fit_csp(self, out_dir):
        """Fit CSP data and save result to file

        Args:
            out_dir (str): Directory where results are saved
        """

        # Define red and blue bandpasses
        blue_bands = ['u', 'g', 'B', 'V0', 'V', 'Y']
        blue_bands = [f'91bg_proj_csp_{f}' for f in blue_bands]
        red_bands = ['r', 'i', 'H', 'J', 'Jrc2', 'Ydw', 'Jdw', 'Hdw']
        red_bands = [f'91bg_proj_csp_{f}' for f in red_bands]

        params = self.fitting_params['csp']
        _run_fit_set(csp.iter_sncosmo_input,
                     out_dir,
                     blue_bands,
                     red_bands,
                     params)

    def fit_des(self, out_dir):
        """Fit DES data and save result to file

        Args:
            out_dir (str): Directory where results are saved
        """

        # Define red and blue bandpasses
        blue_bands = ['desg', 'desr']
        red_bands = ['desi', 'desz', 'desy']

        params = self.fitting_params['des']
        _run_fit_set(des.iter_sncosmo_input,
                     out_dir,
                     blue_bands,
                     red_bands,
                     params)

    def fit_sdss(self, out_dir):
        """Fit SDSS data and save result to file

        Args:
            out_dir (str): Directory where results are saved
        """

        # Define red and blue bandpasses
        blue_bands = [f'91bg_proj_sdss_{b}{c}' for b, c in product('ug', '123456')]
        red_bands = [f'91bg_proj_sdss_{b}{c}' for b, c in product('riz', '123456')]

        params = self.fitting_params['sdss']
        _run_fit_set(sdss.iter_sncosmo_input,
                     out_dir,
                     blue_bands,
                     red_bands,
                     params)
