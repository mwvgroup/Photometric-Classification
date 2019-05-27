#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Provides a class based interface for iteratively fitting light-curves for
multiple surveys with multiple models.

The LCFitting class provides survey specific fitting functions and is capable
of fitting light-curves exclusively in the rest-frame blue or red.
"""

import numpy as np
import yaml

from ._fit_n_params import fit_n_params
from ..data_access import csp, des, sdss


def iter_all_fits(out_dir, module, models, num_params, kwargs):
    """Iteratively fit data for a given survey

    Args:
        out_dir          (str): Directory to write fit results to
        module        (module): A data access module
        models   (list[Model]): List of models to fit
        num_params (list[int]): Number of params to fit
        kwargs          (dict): A dictionary of kwargs for nest_lc AND fit_lc
    """

    for model in models:
        for num_params in num_params:
            # Get kwargs
            model_key = f'{model.source.name}_{model.source.version}'
            kwargs_this = kwargs[model_key][num_params]

            kwargs_this['warn'] = kwargs_this.get('warn', False)
            fit_n_params(out_dir, num_params, module, model, kwargs_this)


class LCFitting:

    def __init__(self, params_file=None):
        """Provides survey specific light-curve fitting functions

        Args:
            params_file (str): The path of a yaml file specifying SNCosmo args.

        Attributes:
            fit_params: Dictionary of settings used when fitting light-curves.

        Methods:
            fit_csp: Iteratively fit CSP data
            fit_des: Iteratively fit DES data
            fit_sdss: Iteratively fit SDSS data
            split_bands: Split band-passes into blue (< 5500 A) and
                         red (>= 5500 A) bands
        """

        self._fitting_params = None
        if params_file is not None:
            with open(params_file) as ofile:
                self._fitting_params = yaml.load(ofile, Loader=yaml.FullLoader)

    @property
    def fit_params(self):
        return self._fitting_params

    @staticmethod
    def split_bands(bands, lambda_eff):
        """Split band-passes into collections of blue and red bands

        Blue bands have an effective wavelength < 5500 Ang. Red bands have an
        effective wavelength >= 5500 Ang.

        Args:
            bands        (array[str]): Name of band-passes
            lambda_eff (array[float]): Effective wavelength of band-passes
        """

        is_blue = np.array(lambda_eff) < 5500
        b_array = np.array(bands)
        return b_array[is_blue], b_array[~is_blue]

    def fit_csp(self, out_dir, models, num_params):
        """Fit CSP data and save result to file

        Acceptable models to fit include 'salt_2_4', 'salt_2_0', and 'sn_91bg'.
        Acceptable number of fit params are 4 and 5.
        Acceptable bands are 'all', 'blue', and 'red'.

        Args:
            out_dir          (str): Directory where results are saved
            models   (list[Model]): List of models to fit
            num_params (list[int]): Number of params to fit
        """

        kwargs = self._fitting_params[csp.survey_name.lower()]
        iter_all_fits(out_dir, csp, models, num_params, kwargs)

    def fit_des(self, out_dir, models, num_params):
        """Fit DES data and save result to file

        Acceptable models to fit include 'salt_2_4', 'salt_2_0', and 'sn_91bg'.
        Acceptable number of fit params are 4 and 5.
        Acceptable bands are 'all', 'blue', and 'red'.

        Args:
            out_dir          (str): Directory where results are saved
            models   (list[Model]): List of models to fit
            num_params (list[int]): Number of params to fit
        """

        kwargs = self._fitting_params[des.survey_name.lower()]
        iter_all_fits(out_dir, des, models, num_params, kwargs)

    def fit_sdss(self, out_dir, models, num_params):
        """Fit SDSS data and save result to file

        Acceptable models to fit include 'salt_2_4', 'salt_2_0', and 'sn_91bg'.
        Acceptable number of fit params are 4 and 5.
        Acceptable bands are 'all', 'blue', and 'red'.

        Args:
            out_dir          (str): Directory where results are saved
            models   (list[Model]): List of models to fit
            num_params (list[int]): Number of params to fit
        """

        kwargs = self._fitting_params[sdss.survey_name.lower()]
        iter_all_fits(out_dir, sdss, models, num_params, kwargs)
