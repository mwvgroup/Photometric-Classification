#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module provides functions for fitting CSP, DES, and SDSS data with a
SALT2 and 91bg model in all bands, the red bands, and the blue bands.
"""

import os
import time
from itertools import product

import numpy as np
import sncosmo
import yaml
from tqdm import tqdm

from . import fit_n_params
from .._sn91bg_model._model import SN91bgSource
from ..data_access import csp, des, sdss

# Get models for fitting
salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))
salt_2_0 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.0'))
sn_91bg = sncosmo.Model(source=SN91bgSource())
models_dict = dict(salt_2_4=salt_2_4, salt_2_0=salt_2_0, sn_91bg=sn_91bg)


class LCFitting:

    def __init__(self, params_file=None):
        self.fitting_params = None
        if params_file is not None:
            self.load_config_file(params_file)

    def load_config_file(self, params_file):
        """Load a yaml configuration file with fit parameters

        Args:
            params_file (Stf): Path of the file to load
        """

        with open(params_file) as ofile:
            self.fitting_params = yaml.load(ofile)

    @staticmethod
    def split_bands(bands, lambda_eff):
        """Split band-passes into blue (< 5500 A) and red (>= 5500 A) bands

        Args:
            bands        (array[str]): Name of band-passes
            lambda_eff (array[float]): Effective wavelength of band-passes
        """

        is_blue = np.array(lambda_eff) < 5500
        b_array = np.array(bands)
        return b_array[is_blue], b_array[~is_blue]

    def _generic_fit(self, module, out_dir, models, num_params, bands,
                     verbose=True, **kwargs):
        """Iterate over a set of light-curve fits using different models,
        number of parameters, and rest frame bands

        Args:
            module        (module): A submodule of the data_access module
            out_dir          (str): Directory where results are saved
            models     (list[str]): Names of models to fit
            num_params (list[int]): Number of params to fit
            bands      (list[str]): List specifying bandpass collections
            verbose         (bool): Whether to show progress bar
            Any other arguments for module.iter_sncosmo_input
        """

        # Define red and blue bandpasses
        blue_bands, red_bands = self.split_bands(
            module.band_names, module.lambda_effective)

        # Create dictionary of only user specified bands
        bands_dict = dict(blue=blue_bands, red=red_bands, all=None)
        bands_to_fit = {b_name: bands_dict[b_name] for b_name in bands}
        models = [models_dict[model] for model in models]

        # Define data for running fits
        params = self.fitting_params[module.__name__.split('.')[-1]]
        path_pattern = '{}_{}param_{}.ecsv'
        modeling_data = product(models, num_params, bands_to_fit.items())

        # Run fit
        for model, num_param, (band_name, band_lists) in modeling_data:
            model_name = model.source.name + '_' + model.source.version
            tqdm.write(f'{num_param} param {model_name} in {band_name} bands')
            time.sleep(0.5)  # Give output time to flush

            model_args = params[model_name][num_param]
            fname = path_pattern.format(model_name, num_param, band_name)
            out_path = os.path.join(out_dir, fname)
            fit_n_params(
                out_path,
                num_params=num_param,
                inputs=module.iter_sncosmo_input(
                    band_lists, verbose=verbose, **kwargs),
                model=model,
                **model_args)

            tqdm.write('\n')

    def fit_csp(self, out_dir, models, num_params, bands, **kwargs):
        """Fit CSP data and save result to file

        Acceptable models to fit include 'salt_2_4', 'salt_2_0', and 'sn_91bg'.
        Acceptable number of fit params are 4 and 5.
        Acceptable bands are 'all', 'blue', and 'red'.

        Args:
            out_dir          (str): Directory where results are saved
            models     (list[str]): Names of models to fit
            num_params (list[int]): Number of params to fit
            bands      (list[str]): List specifying bandpass collections
            Any other arguments for csp.iter_sncosmo_input
        """

        self._generic_fit(csp, out_dir, models, num_params, bands, **kwargs)

    def fit_des(self, out_dir, models, num_params, bands, **kwargs):
        """Fit DES data and save result to file

        Acceptable models to fit include 'salt_2_4', 'salt_2_0', and 'sn_91bg'.
        Acceptable number of fit params are 4 and 5.
        Acceptable bands are 'all', 'blue', and 'red'.

        Args:
            out_dir          (str): Directory where results are saved
            models     (list[str]): Names of models to fit
            num_params (list[int]): Number of params to fit
            bands      (list[str]): List specifying bandpass collections
            Any other arguments for des.iter_sncosmo_input
        """

        self._generic_fit(des, out_dir, models, num_params, bands, **kwargs)

    def fit_sdss(self, out_dir, models, num_params, bands, **kwargs):
        """Fit SDSS data and save result to file

        Acceptable models to fit include 'salt_2_4', 'salt_2_0', and 'sn_91bg'.
        Acceptable number of fit params are 4 and 5.
        Acceptable bands are 'all', 'blue', and 'red'.

        Args:
            out_dir          (str): Directory where results are saved
            models     (list[str]): Names of models to fit
            num_params (list[int]): Number of params to fit
            bands      (list[str]): List specifying bandpass collections
            Any other arguments for sdss.iter_sncosmo_input
        """

        self._generic_fit(sdss, out_dir, models, num_params, bands, **kwargs)
