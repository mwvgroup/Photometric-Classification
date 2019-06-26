#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module handles the iterative fitting of light-curves for multiple
models, surveys, and band pass collections.
"""

from copy import deepcopy
from pathlib import Path

import numpy as np
from astropy.table import Table

from ._fit_funcs import create_empty_summary_table, fit_lc, get_sampled_model


def _split_bands(bands, lambda_eff):
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


def split_data(data_table, band_names, lambda_eff):
    """Split a data table into blue and red data (by rest frame)

    Split data by keeping filters that are red-ward or blue-ward of 5500 Ang.

    Args:
        data_table (Table): An SNCosmo input table with column 'band'
        band_names  (iter): List of all bands available in the survey
        lambda_eff  (iter): The effective wavelength of each band in band_names

    Returns:
        A SNCosmo input table with only blue bands
        A SNCosmo input table with only red bands
    """

    # Type cast to allow numpy indexing
    band_names = np.array(band_names)
    lambda_eff = np.array(lambda_eff)

    @np.vectorize
    def lambda_for_band(band):
        return lambda_eff[band_names == band]

    # Rest frame effective wavelengths for each observation
    z = data_table.meta['redshift']
    observed_lambda = lambda_for_band(data_table['band'])
    rest_frame_lambda = observed_lambda / (1 + z)

    # Get the name of the observer frame band with the smallest distance
    # to each rest frame lambda
    delta_lambda = np.array([
        np.abs(rest_frame_lambda - l_eff) for l_eff in lambda_eff])

    min_indx = np.argmin(delta_lambda, axis=0)
    rest_frame_filters = np.array(band_names)[min_indx]

    # Keep only the specified filters that are within 700 Angstroms of the
    # rest frame effective wavelength
    is_ok_diff = delta_lambda[min_indx, np.arange(delta_lambda.shape[1])] < 700

    # Split into blue and red band passes
    out_list = []
    for bands in _split_bands(band_names, lambda_eff):
        is_in_bands = np.isin(rest_frame_filters, bands)
        indices = np.logical_and(is_in_bands, is_ok_diff)
        out_list.append(data_table[indices])

    return out_list


def _create_table_paths(out_dir, model, survey, num_params):
    """Create a list of file paths for all, blue, and red band fit results

    Args:
        out_dir (str): Out directory of file paths
        model (Model): SNCosmo model
        survey  (str): Name of a survey

    Returns:
        ["<out_dir>/<model>_<model version><survey>_<bands>", ...]
    """

    def create_fname(band):
        model_name = f'{model.source.name}_{model.source.version}'
        return f'{survey}_{num_params}_{model_name}_{band}.ecsv'

    out_dir = Path(out_dir)
    return [out_dir / create_fname(band) for band in ('all', 'blue', 'red')]


def _iter_fit_bands(out_dir, module, model, params_to_fit, kwargs, verbose,
                    time_out, skip_types=()):
    """Fit light curves from a given data module in all bandpass collections

    Args:
        out_dir             (str): Directory of output files
        module           (module): data access module for a particular survey
        model             (Model): SNCosmo model
        params_to_fit (list[str]): Parameters to fit for
        kwargs             (dict): Any arguments for nest_lc and fit_lc
        verbose            (dict): Optional arguments for tqdm progress bar
        time_out            (int): Seconds before nested sampling times out
        skip_types    (list[str]): Classifications to skip if provided by survey
    """

    Path(out_dir).mkdir(exist_ok=True)

    # Create separate tables for each band's fit results
    out_tables = [create_empty_summary_table() for _ in range(3)]
    out_paths = _create_table_paths(
        out_dir,
        model,
        module.survey_abbrev.lower(),
        len(params_to_fit))

    filter_func = lambda x: x.meta.get('classification', '') not in skip_types
    iterable = module.iter_data(
        verbose=verbose, format_sncosmo=True, filter_func=filter_func)

    for data in iterable:
        model_this = deepcopy(model)
        kwargs_this = deepcopy(kwargs)

        # If not fitting for redshift and it is available, set z
        if 'z' not in params_to_fit:
            z = data.meta.get('redshift', -9)
            if z <= 0:
                continue

            model_this.set(z=z)

        # Fit light-curves
        sampled_model = get_sampled_model(
            survey_name=module.survey_abbrev,
            data=data,
            model=model_this,
            vparam_names=params_to_fit,
            time_out=time_out,
            **kwargs_this)

        kwargs_this['guess_amplitude'] = False
        kwargs_this['guess_t0'] = False
        kwargs_this['guess_z'] = False

        blue, red = split_data(data, module.band_names,
                               module.lambda_effective)
        iter_data = zip(out_tables, out_paths, [data, blue, red])
        for table, path, input_table in iter_data:
            fit_results = fit_lc(
                data=input_table,
                model=deepcopy(sampled_model),
                vparam_names=params_to_fit,
                **deepcopy(kwargs_this))

            table.add_row(fit_results)
            table.write(path, overwrite=True)


def _fit_n_params(out_dir, num_params, module, model, kwargs,
                  time_out, skip_types=()):
    """Fit light curves from in all bandpass collections with 4 or 5 parameters

    Args:
        out_dir          (str): Directory of output files
        num_params       (int): Number of parameters to fit (4 or 5)
        module        (module): data access module for a particular survey
        model          (Model): SNCosmo model
        kwargs          (dict): Any arguments for nest_lc and fit_lc
        time_out         (int): Seconds before nested sampling times out
        skip_types (list[str]): Classifications to skip if provided by survey
    """

    if num_params not in (4, 5):
        raise ValueError('num_params argument must be 4 or 5.')

    # Get parameter names
    param_names = ['z', 't0', 'x0', 'x1', 'c']
    params_to_fit = param_names[5 - num_params:]

    # Create progress bar options
    pbar_txt = (f'{num_params} param '
                f'{model.source.name} - '
                f'{model.source.version} for '
                f'{module.survey_abbrev}')

    pbar_args = {'desc': pbar_txt, 'position': 1}

    _iter_fit_bands(
        out_dir=out_dir,
        module=module,
        model=model,
        params_to_fit=params_to_fit,
        kwargs=kwargs,
        verbose=pbar_args,
        time_out=time_out,
        skip_types=skip_types)


def iter_all_fits(out_dir, module, models, num_params, kwargs,
                  time_out=None, skip_types=()):
    """Iteratively fit data for a given survey

    Args:
        out_dir          (str): Directory to write fit results to
        module        (module): A data access module
        models   (list[Model]): List of models to fit
        num_params (list[int]): Number of params to fit
        kwargs          (dict): A dictionary of kwargs for nest_lc AND fit_lc
        time_out         (int): Seconds before nested sampling times out
        skip_types (list[str]): Classifications to skip if provided by survey
    """

    for model in models:
        for n in num_params:
            # Get kwargs
            model_key = f'{model.source.name}_{model.source.version}'
            kwargs_this = kwargs[model_key]
            kwargs_this['warn'] = kwargs_this.get('warn', False)

            _fit_n_params(
                out_dir=out_dir,
                num_params=n,
                module=module,
                model=model,
                kwargs=kwargs_this,
                time_out=time_out,
                skip_types=skip_types)
