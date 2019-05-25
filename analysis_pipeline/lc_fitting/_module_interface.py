#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Provides a class based interface for iteratively fitting light-curves for
multiple surveys with multiple models.

The LCFitting class provides survey specific fitting functions and is capable
of fitting light-curves exclusively in the rest-frame blue or red.
"""

import numpy as np
import sncosmo
import yaml

from ._fit_n_params import fit_n_params
from .._sn91bg_model._model import SN91bgSource
from ..data_access import csp, des, sdss

# Get models for fitting
salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))
salt_2_0 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.0'))
sn_91bg = sncosmo.Model(source=SN91bgSource())
models_dict = dict(salt_2_4=salt_2_4, salt_2_0=salt_2_0, sn_91bg=sn_91bg)


def iter_all_fits(out_dir, module, kwargs):
    """Iteratively fit data for a given survey

    Args:
        out_dir   (str): Directory to write fit results to
        module (module): A data access module
        nkwargs  (dict): A dictionary of kwargs for sncosmo.nest_lc
        fkwargs  (dict): A dictionary of kwargs for sncosmo.fit_lc
    """

    salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))
    salt_2_0 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.0'))
    sn_91bg = sncosmo.Model(source=SN91bgSource())
    models_dict = {
        'salt2_2.4': salt_2_4,
        'salt2_2.0': salt_2_0,
        'sn91bg_color_interpolation': sn_91bg}

    for model_name, model in models_dict.items():
        for num_params in (4, 5):
            # Get kwargs
            kwargs_this = kwargs[module.survey_name.lower()][model_name][
                num_params]

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
                self._fitting_params = yaml.load(ofile)

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

    def fit_csp(self, out_dir):
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
            nest (bool): Use nested sampling to determine initial guess values
        """

        iter_all_fits(out_dir, csp, self._fitting_params)

    def fit_des(self, out_dir):
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
            nest (bool): Use nested sampling to determine initial guess values
        """

        iter_all_fits(out_dir, des, self._fitting_params)

    def fit_sdss(self, out_dir):
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
            nest (bool): Use nested sampling to determine initial guess values
        """

        iter_all_fits(out_dir, sdss, self._fitting_params)
