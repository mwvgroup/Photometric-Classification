#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module provides functions for fitting light curves using sncosmo.

The fit_lc function acts as a wrapper for sncosmo.fit_lc that returns fit
results as a list.

The fit_n_params function iteratively fits a collection of light-curves using a
4 or 5 parameter fit.

The LCFitting class provides survey specific fitting functions and is capable
of fitting light-curves exclusively in the rest-frame blue or red.
"""

import os
from copy import deepcopy
from itertools import product

import numpy as np
import sncosmo
import yaml
from astropy.table import Table
from sncosmo.fitting import DataQualityError
from tqdm import tqdm

from ._sn91bg_model._model import SN91bgSource
from .data_access import csp, des, sdss

# Get models for fitting
salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))
salt_2_0 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.0'))
sn_91bg = sncosmo.Model(source=SN91bgSource())
models_dict = dict(salt_2_4=salt_2_4, salt_2_0=salt_2_0, sn_91bg=sn_91bg)


def _create_empty_summary_table():
    """Returns a table with columns:

         cid, num_points, z, t0, x0, x1, z_err, t0_err, x0_err,
         x1_err, c_err, chi, dof, tmin, tmax, pre_max, post_max, message
    """

    names = ['cid', 'num_points']
    dtype = ['U20', int]

    param_names = ('z', 't0', 'x0', 'x1', 'c')
    names.extend(param_names)
    names.extend((p + '_err' for p in param_names))
    dtype.extend((float for _ in range(len(names) - 2)))

    names.extend(
        ('chi', 'dof', 'b_max',
         'delta_15', 'tmin', 'tmax',
         'pre_max', 'post_max', 'message'))

    dtype.extend([float, float, float, float, float, float, int, int, 'U1000'])
    return Table(names=names, dtype=dtype)


def _count_pre_and_post_max(obs_times, t_max):
    """Count the number of light-curve data points pre and post maximum

    Args:
        obs_times (list[float]): Times of observations
        t_max (float): Time of maximum

    Returns:
        The minimum time value
        The maximum time value
        The number of data points before maximum
        The number of data points after maximum
    """

    times_arr = np.array(obs_times)
    pre_max = sum(times_arr < t_max)
    post_max = sum(times_arr > t_max)
    return min(times_arr), max(times_arr), pre_max, post_max


def fit_lc(data, model, vparam_names, nest=False, **kwargs):
    """A wrapper for sncosmo.fit_lc that returns results as a list

    Exceptions raised by sncosmo.fit_lc are caught and stored as the exit
    message. By default, initial parameters are determine using nested
    sampling. Mutable arguments are not guaranteed to go unchanged.

    Args:
        Any arguments for sncosmo.fit_lc
        nest (bool): Whether to use nested sampling to determine initial values

    Returns:
        A list of values for 'z', 't0', 'x0', 'x1', 'c', their respective
        errors, the fit chi-squared, number of DOF, and SNCosmo exit message.
    """

    try:
        if nest:
            nest_result, _ = sncosmo.nest_lc(
                data, model, vparam_names,
                bounds=kwargs['bounds'],
                verbose=True,
                maxiter=2000
            )

            # Set initial parameters in model
            model_args = dict()
            for vp, p in zip(nest_result.vparam_names, nest_result.parameters):
                model_args[vp] = p

            model.set(**model_args)
            kwargs['guess_amplitude'] = False
            kwargs['guess_t0'] = False
            kwargs['guess_z'] = False

        result, fitted_model = sncosmo.fit_lc(
            data=data, model=model, vparam_names=vparam_names, **kwargs)

    # If the fit fails fill out_data with place holder values (NANs and zeros)
    except (DataQualityError, RuntimeError, ValueError) as e:
        out_data = np.full(16, np.NAN).tolist()
        out_data.extend((0, 0, str(e).replace('\n', ' ')))
        if 'z' not in vparam_names:
            out_data[0] = data.meta['redshift']
            out_data[5] = data.meta.get('redshift_err', np.NAN)

    else:
        params = {p: v for p, v in zip(result.param_names, result.parameters)}

        # Create list of parameter values
        out_data = [params[p] for p in ('z', 't0', 'x0', 'x1', 'c')]
        for param in ('z', 't0', 'x0', 'x1', 'c'):
            out_data.append(result.errors.get(param, 0))

        # Determine peak magnitude and decline rate
        b_max = fitted_model.source_peakabsmag('bessellb', 'ab')
        peak_phase = fitted_model.source.peakphase('bessellb')
        b_0 = fitted_model.source.bandmag('bessellb', 'ab', peak_phase)
        b_15 = fitted_model.source.bandmag('bessellb', 'ab', peak_phase + 15)
        delta_15 = b_15 - b_0

        out_data.append(result.chisq)
        out_data.append(result.ndof)
        out_data.append(b_max)
        out_data.append(delta_15)
        out_data.extend(_count_pre_and_post_max(data['time'], params['t0']))
        out_data.append(result.message.replace('\n', ' '))

    return out_data


def fit_n_params(
        out_path, num_params, inputs, model, warn=False, nest=False, **kwargs):
    """Fit light curves with a 4 parameter Salt2-like model using SNCosmo

    Redshift values are taken from the meta data of input tables using the
    'redshift' key. Inputs with negative, false, or missing redshift values
    are skipped.

    Args:
        out_path       (str): Where to write fit results
        num_params     (int): Number of parameters to fit. Either 4 or 5.
        inputs (iter[Table]): Iterable of SNCosmo input tables
        model        (model): SNCosmo model to use for fitting
        warn          (bool): Show sncosmo warnings (default = False)
        nest (bool): Use nested sampling to determine initial guess values

        Additionally any arguments for sncosmo.fit_lc not mentioned above
    """

    if num_params not in (4, 5):
        raise ValueError("Parameter 'num_params' must be either 4 or 5")

    out_dir = os.path.dirname(out_path)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    kwargs['warn'] = warn
    out_table = _create_empty_summary_table()
    out_table.meta = kwargs

    # Define list of parameters to fit
    vparam_names = ['z', 't0', 'x0', 'x1', 'c']
    if num_params == 4:
        del vparam_names[0]

    # Run fit for each target
    for input_table in inputs:
        model_this = deepcopy(model)
        kwargs_this = deepcopy(kwargs)

        # Set redshift in model for 4 param fits
        if num_params == 4:
            z = input_table.meta.get('redshift', False)
            if z < 0 or not z:
                continue

            model_this.set(z=z)

        fit_results = fit_lc(
            data=input_table,
            model=model_this,
            vparam_names=vparam_names,
            nest=nest,
            **kwargs_this)

        new_row = [input_table.meta['cid'], len(input_table)]
        new_row.extend(fit_results)
        out_table.add_row(new_row)
        out_table.write(out_path, overwrite=True)


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

    def _generic_fit(self, module, out_dir, models, num_params, bands,
                     nest=False, **kwargs):
        """Iterate over a set of light-curve fits using different models,
        number of parameters, and rest frame bands

        Args:
            module        (module): A submodule of the data_access module
            out_dir          (str): Directory where results are saved
            models     (list[str]): Names of models to fit
            num_params (list[int]): Number of params to fit
            bands      (list[str]): List specifying bandpass collections
            Any other arguments for module.iter_sncosmo_input
        """

        # Define red and blue band-passes
        blue_bands, red_bands = self.split_bands(
            module.band_names, module.lambda_effective)

        # Create dictionary of only user specified bands
        bands_dict = dict(blue=blue_bands, red=red_bands, all=None)
        bands_to_fit = {b_name: bands_dict[b_name] for b_name in bands}

        # Get model instances to use for fitting
        models = (models_dict[model] for model in models)

        # Define data for running fits
        params = self._fitting_params[module.__name__.split('.')[-1]]
        path_pattern = '{}_{}param_{}.ecsv'
        modeling_data = product(models, num_params, bands_to_fit.items())

        # Run fits
        for model, num_param, (band_color, band_list) in modeling_data:
            model_name = model.source.name + '_' + model.source.version
            fname = path_pattern.format(model_name, num_param, band_color)

            pbar_txt = f'{num_param} param {model_name} in {band_color} bands'
            pbar_pos = 1 if nest else 0
            inputs = module.iter_sncosmo_input(
                band_list,
                verbose={'desc': pbar_txt, 'position': pbar_pos},
                **kwargs)

            model_args = params[model_name][num_param]
            fit_n_params(
                out_path=os.path.join(out_dir, fname),
                num_params=num_param,
                inputs=inputs,
                model=model,
                nest=nest,
                **model_args)

            tqdm.write('\n')

    def fit_csp(
            self, out_dir, models, num_params, bands, nest=False, **kwargs):
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

        self._generic_fit(
            csp, out_dir, models, num_params, bands, nest=nest, **kwargs)

    def fit_des(
            self, out_dir, models, num_params, bands, nest=False, **kwargs):
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

        self._generic_fit(
            des, out_dir, models, num_params, bands, nest=nest, **kwargs)

    def fit_sdss(
            self, out_dir, models, num_params, bands, nest=False, **kwargs):
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

        self._generic_fit(
            sdss, out_dir, models, num_params, bands, nest=nest, **kwargs)
