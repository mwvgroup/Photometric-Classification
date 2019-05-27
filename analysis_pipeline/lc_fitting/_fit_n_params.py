#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module handles the iterative fitting of light-curves for multiple
models, surveys, and band pass collections.
"""

from copy import deepcopy
from pathlib import Path

import numpy as np
from astropy.table import Table

from ._fit_funcs import create_empty_summary_table, fit_lc, nest_lc


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


def _split_data(data_table, band_names, lambda_eff):
    """Split a data table into blue and red data (by restframe)

    Split data by keeping filters that are redward or blueward of 5500 Ang.

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


def _create_table_paths(out_dir, model, survey):
    """Create a list of file paths for all, blue, and red band fit results

    Args:
        out_dir (str): Out directory of file paths
        model (Model): SNCosmo model
        survey  (str): Name of a survey

    Returns:
        ["<out_dir>/<model>_<model version><survey>_<bands>", ...]
    """

    def create_fname(band):
        return f'{model.source.name}_{model.source.version}_{survey}_{band}.csv'

    out_dir = Path(out_dir)
    return [out_dir / create_fname(band) for band in ('all', 'blue', 'red')]


def _iter_fit_bands(out_dir, module, model, params_to_fit, kwargs, verbose):
    """Fit light curves from a given data module in all bandpass collections

    Args:
        out_dir             (str): Directory of output files
        module           (module): data access module for a particular survey
        model             (Model): SNCosmo model
        params_to_fit (list[str]): Parameters to fit for
        kwargs             (dict): Any arguments for nest_lc and fit_lc
        verbose            (dict): Optional arguments for tqdm progress bar
    """

    Path(out_dir).mkdir(exist_ok=True)

    # Create separate tables for each band's fit results
    out_tables = [create_empty_summary_table() for _ in range(3)]
    out_paths = _create_table_paths(out_dir, model, module.survey_name)
    for data in module.iter_sncosmo_input(verbose=verbose):

        model_s = deepcopy(model)
        # If not fitting for redshift and it is available, set z
        if 'z' not in params_to_fit:
            z = data.meta.get('redshift', -9)
            if z <= 0:
                continue

            model_s.set(z=z)

        sampled_model = nest_lc(
            data, model_s, params_to_fit, **deepcopy(kwargs))

        # Fit light-curves
        blue, red = _split_data(data, module.band_names,
                                module.lambda_effective)
        iter_data = zip(out_tables, out_paths, [data, blue, red])
        for table, path, input_table in iter_data:
            kwargs['guess_amplitude'] = False
            kwargs['guess_t0'] = False
            kwargs['guess_z'] = False

            fit_results = fit_lc(
                data=input_table,
                model=deepcopy(sampled_model),
                vparam_names=params_to_fit,
                **deepcopy(kwargs))

            table.add_row(fit_results)
            table.write(path, overwrite=True)


def fit_n_params(out_dir, num_params, module, model, kwargs):
    """Fit light curves from in all bandpass collections with 4 or 5 parameters

    Args:
        out_dir    (str): Directory of output files
        num_params (int): Number of parameters to fit (4 or 5)
        module  (module): data access module for a particular survey
        model    (Model): SNCosmo model
        kwargs    (dict): Any arguments for nest_lc and fit_lc
    """

    if num_params not in (4, 5):
        raise ValueError('num_params argument must be 4 or 5.')

    # Get parameter names
    param_names = ['z', 't0', 'x0', 'x1', 'c']
    params_to_fit = param_names[5 - num_params:]

    # Create progress bar options
    pbar_txt = f'{num_params} param fit for {module.survey_name}'
    pbar_args = {'desc': pbar_txt, 'position': 1}

    _iter_fit_bands(out_dir, module, model, params_to_fit, kwargs, pbar_args)
