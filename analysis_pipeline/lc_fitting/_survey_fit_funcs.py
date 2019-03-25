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
from .._sn91bg_model._model import SN91bgSource

# Get models for fitting
salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))
salt_2_0 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.0'))
sn_91bg = sncosmo.Model(source=SN91bgSource())
models_dict = dict(salt_2_4=salt_2_4, salt_2_0=salt_2_0, sn_91bg=sn_91bg)


def _run_fit_set(data, out_dir, models, num_params, bands_to_fit, params, **kwargs):
    """
    Iterate over a set of light-curve fits using different models, number of
    parameters, and rest frame bands

    Args:
        data       (module): A submodule of the data_access module (eg. csp)
        out_dir       (str): Directory where results are saved
        num_params   (List): List of SNCosmo models to fit
        bands_to_fit (dict): A dictionary of bands for fit (eg. {'blue': 'ugriz'})
        params       (dict): Dictionary of SNCosmo fitting arguments
    """

    modeling_data = product(models, num_params, bands_to_fit.items())

    path_pattern = '{}_{}param_{}.ecsv'
    for model, num_param, (band_name, band_lists) in modeling_data:
        model_name = model.source.name + '_' + model.source.version
        tqdm.write(f'{num_param} param {model_name} in {band_name} bands')
        time.sleep(0.5)  # Give output time to flush

        model_args = params[model_name][num_param]
        fname = path_pattern.format(model_name, num_param, band_name)
        out_path = os.path.join(out_dir, fname)
        fit_n_params(out_path,
                     num_params=num_param,
                     inputs=data.iter_sncosmo_input(band_lists, **kwargs),
                     bands=data.band_names,
                     model=model,
                     **model_args)

        tqdm.write('\n')


class LCFitting:

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

    def fit_csp(self, out_dir, models, num_params, bands, **kwargs):
        """Fit CSP data and save result to file

        Args:
            out_dir (str): Directory where results are saved
        """

        # Define red and blue bandpasses
        blue_bands = ['u', 'g', 'B', 'V0', 'V', 'Y']
        blue_bands = [f'91bg_proj_csp_{f}' for f in blue_bands]
        red_bands = ['r', 'i', 'H', 'J', 'Jrc2', 'Ydw', 'Jdw', 'Hdw']
        red_bands = [f'91bg_proj_csp_{f}' for f in red_bands]
        bands_dict = dict(blue=blue_bands,
                          red=red_bands,
                          all=blue_bands + red_bands)

        # Create dictionary of only user specified bands
        bands_to_fit = {b_name: bands_dict[b_name] for b_name in bands}
        models = [models_dict[model] for model in models]

        # Run fit
        params = self.fitting_params['csp']
        _run_fit_set(csp,
                     out_dir,
                     models,
                     num_params,
                     bands_to_fit,
                     params,
                     **kwargs)

    def fit_des(self, out_dir, models, num_params, bands, **kwargs):
        """Fit DES data and save result to file

        Args:
            out_dir (str): Directory where results are saved
        """

        # Define red and blue bandpasses
        blue_bands = ['desg', 'desr']
        red_bands = ['desi', 'desz', 'desy']
        bands_dict = dict(blue=blue_bands, red=red_bands,
                          all=blue_bands + red_bands)

        # Create dictionary of only user specified bands
        bands_to_fit = {b_name: bands_dict[b_name] for b_name in bands}
        models = [models_dict[model] for model in models]

        # Run fit
        params = self.fitting_params['des']
        _run_fit_set(des,
                     out_dir,
                     models,
                     num_params,
                     bands_to_fit,
                     params,
                     **kwargs)

    def fit_sdss(self, out_dir, models, num_params, bands, **kwargs):
        """Fit SDSS data and save result to file

        Args:
            out_dir (str): Directory where results are saved
        """

        # Define red and blue bandpasses
        blue_bands = [f'91bg_proj_sdss_{b}{c}' for b, c in product('ug', '123456')]
        red_bands = [f'91bg_proj_sdss_{b}{c}' for b, c in product('riz', '123456')]
        bands_dict = dict(blue=blue_bands,
                          red=red_bands,
                          all=blue_bands + red_bands)

        # Create dictionary of only user specified bands
        bands_to_fit = {b_name: bands_dict[b_name] for b_name in bands}
        models = [models_dict[model] for model in models]

        # Run fit
        params = self.fitting_params['sdss']
        _run_fit_set(sdss,
                     out_dir,
                     models,
                     num_params,
                     bands_to_fit,
                     params,
                     **kwargs)
