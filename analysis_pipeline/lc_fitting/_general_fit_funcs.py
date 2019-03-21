#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module provides functions for fitting light curves using sncosmo. It
provides the fit_n_params function which acts as a wrapper for the
sncosmo.fit_lc function. The wrapper adds the ability to iterate over data
provided by the data_access module and save results to file.
"""

import os
from copy import deepcopy

import numpy as np
import sncosmo
from astropy.table import Table
from sncosmo.fitting import DataQualityError


def _create_empty_summary_table(band_names):
    """Returns a table with columns:

         cid, *num_points_<band_names>, z, t0, x0, x1, z_err, t0_err, x0_err,
         x1_err, c_err, chi, dof, message.
    """

    names = ['cid']
    names.extend(['num_points_' + band for band in band_names])

    param_names = ['z', 't0', 'x0', 'x1', 'c']
    names.extend(param_names)
    names.extend((p + '_err' for p in param_names))
    names.extend(('chi', 'dof', 'message'))
    out_table = Table(names=names, dtype=[object for _ in names])

    return out_table


def _count_points_per_band(band_list, all_band_names):
    """Determine number of data points per band

    count the number of times each element in <all_band_names> appears in
    <band_list>.

    Returns:
        A list with the number of counts for each element in <all_band_names>
    """

    band_names, band_counts = np.unique(band_list, return_counts=True)
    count_dict = dict(zip(band_names, band_counts))
    return [count_dict.get(band, 0) for band in all_band_names]


def _run_fit(data, model, vparam_names, **kwargs):
    """Run a light curve fit

    Args:
        Any arguments for sncosmo.fit_lc

    Returns:
        A list of values for 'z', 't0', 'x0', 'x1', 'c', their respective
        errors, the fit chi-squared, number of DOF, and SNCosmo exit message.
    """

    out_data = []
    try:
        result, fitted_model = sncosmo.fit_lc(
            data=data, model=model, vparam_names=vparam_names, **kwargs)

    except (DataQualityError, RuntimeError, ValueError) as e:
        if 'z' in vparam_names:
            out_data.extend(np.full(12, np.NAN).tolist())
            out_data.append(str(e))

        else:
            z = data.meta['redshift']
            z_err = data.meta.get('redshift_err', 0)

            out_data.append(z)
            out_data.extend(np.full(4, np.NAN).tolist())
            out_data.append(z_err)
            out_data.extend(np.full(6, np.NAN).tolist())
            out_data.append(str(e))

    else:
        for param in ['z', 't0', 'x0', 'x1', 'c']:
            i = result.param_names.index(param)
            out_data.append(result.parameters[i])

        for param in ['z', 't0', 'x0', 'x1', 'c']:
            out_data.append(result.errors.get(param, 0))

        out_data.append(result.chisq)
        out_data.append(result.ndof)
        out_data.append(result.message)

    return out_data


def fit_n_params(out_path, num_params, inputs, bands, model, warn=False, **kwargs):
    """Fit light curves with a 4 parameter Salt2-like model using SNCosmo

    Redshift values are taken from the meta data of input tables using the
    'redshift' key. Inputs with negative, false, or missing redshift values
    are skipped.

    Args:
        out_path       (str): Where to write fit results
        num_params     (int): Number of parameters to fit. Either 4 or 5.
        inputs (iter[Table]): Iterable of SNCosmo input tables
        bands    (list[str]): The bands of the survey being fitted
        model        (model): SNCosmo model to use for fitting
        warn          (bool): Show sncosmo warnings (default = False)

        Additionally any arguments for sncosmo.fit_lc not mentioned above
    """

    if num_params not in (4, 5):
        raise ValueError("Parameter 'num_params' must be either 4 or 5")

    out_dir = os.path.dirname(out_path)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    kwargs['warn'] = warn
    out_table = _create_empty_summary_table(bands)

    # Run fit for each target
    for input_table in inputs:
        # Protect against mutating input args
        model_this = deepcopy(model)
        kwargs_this = deepcopy(kwargs)

        # Set redshift in model for 4 param fits
        if num_params == 4:
            z = input_table.meta.get('redshift', False)
            if z < 0 or not z:
                continue

            model_this.set(z=z)

        band_count = _count_points_per_band(input_table['band'], bands)
        fit_results = _run_fit(
            data=input_table,
            model=model_this,
            vparam_names=['t0', 'x0', 'x1', 'c'],
            **kwargs_this)

        new_row = [input_table.meta['cid']]
        new_row.extend(band_count)
        new_row.extend(fit_results)
        out_table.add_row(new_row)
        out_table.write(out_path, overwrite=True)
