#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module handles the iterative fitting of light-curves for multiple
models and band pass collections.
"""

from copy import deepcopy

from ._fit_funcs import create_results_table, fit_lc
from ._sampling import get_sampled_model
from .._paths import get_fit_result_paths
from ..utils import split_data


def iter_fit_bands(survey, model, params_to_fit, kwargs, verbose,
                   time_out, skip_types=()):
    """Fit light curves from a given data module in all bandpass collections

    Args:
        survey           (module): data access module for a particular survey
        model             (Model): SNCosmo model
        params_to_fit (list[str]): Parameters to fit for
        kwargs             (dict): Any arguments for nest_lc and fit_lc
        verbose            (dict): Optional arguments for tqdm progress bar
        time_out            (int): Seconds before nested sampling times out
        skip_types    (list[str]): Classifications to skip if provided by survey
    """

    # Create separate tables for each band's fit results
    out_tables = [create_results_table(model.param_names) for _ in range(3)]
    out_paths = get_fit_result_paths(model, survey, len(params_to_fit))

    # Create iterable of input data tables to be fitted
    filter_func = lambda x: x.meta.get('classification', '') not in skip_types
    iterable = survey.iter_data(
        verbose=verbose, format_sncosmo=True, filter_func=filter_func)

    for data in iterable:
        kwargs_this = deepcopy(kwargs)

        # If not fitting for redshift and it is available, set z
        if 'z' not in params_to_fit:
            z = data.meta.get('redshift', -9)
            if z <= 0:
                continue

            model.set(z=z)

        # Fit light-curves
        sampled_model, bounds = get_sampled_model(
            survey=survey,
            data=data,
            model=model,
            vparam_names=params_to_fit,
            time_out=time_out,
            **kwargs_this)

        kwargs_this['guess_amplitude'] = False
        kwargs_this['guess_t0'] = False
        kwargs_this['guess_z'] = False

        blue, red = split_data(data, survey.band_names, survey.lambda_effective)
        iter_data = zip(out_tables, out_paths, [data, blue, red])
        for table, path, input_table in iter_data:
            fit_results = fit_lc(
                data=input_table,
                model=sampled_model,
                vparam_names=params_to_fit,
                **kwargs_this)

            table.add_row(fit_results)
            table.write(path, overwrite=True)


def run_iter_fitting(survey, model_list, kwarg_list, fitz_list,
                     time_out, skip_types=()):
    """Iteratively fit data for a given survey with multiple models

    Args:
        survey          (module): A data access module
        model_list (list[Model]): List of models to fit
        fit_z       (list[bool]): Whether to fit for redshift
        kwarg_list   (list[dict]): A dictionary of kwargs for nest_lc AND fit_lc
        time_out            (int): Seconds before nested sampling times out
        skip_types    (list[str]): Classifications to skip if provided by survey
    """

    iter_args = zip(model_list, kwarg_list, fitz_list)
    for model, kwargs, fit_z in iter_args:

        # Get parameter names
        param_names = deepcopy(model.param_names)
        if not fit_z:
            param_names.pop(param_names.index('z'))

        # Create progress bar options
        pbar_txt = (f'{len(param_names)} param '
                    f'{model.source.name} - '
                    f'{model.source.version} for '
                    f'{survey.survey_abbrev}')

        pbar_args = {'desc': pbar_txt, 'position': 1}

        iter_fit_bands(
            survey=survey,
            model=model,
            params_to_fit=param_names,
            kwargs=kwargs,
            verbose=pbar_args,
            time_out=time_out,
            skip_types=skip_types)
