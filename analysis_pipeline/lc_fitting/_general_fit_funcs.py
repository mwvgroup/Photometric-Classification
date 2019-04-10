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


def _create_empty_summary_table():
    """Returns a table with columns:

         cid, num_points, z, t0, x0, x1, z_err, t0_err, x0_err,
         x1_err, c_err, chi, dof, tmin, tmax, pre_max, post_max, message
    """

    names = ['cid', 'num_points']
    dtype = ['U20', int]

    param_names = ['z', 't0', 'x0', 'x1', 'c']
    names.extend(param_names)
    names.extend((p + '_err' for p in param_names))
    dtype.extend([float for _ in range(2 * len(param_names))])

    names.extend(('chi', 'dof', 'b_max', 'delta_15', 'tmin', 'tmax',
                  'pre_max', 'post_max', 'message'))

    dtype.extend([float, float, float, float, float, float, int, int, 'U250'])

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


def fit_lc(data, model, vparam_names, **kwargs):
    """A wrapper for sncosmo.fit_lc that returns results as a list

    Exceptions raised by sncosmo.fit_lc are caught and stored as the exit
    message.

    Args:
        Any arguments for sncosmo.fit_lc

    Returns:
        A list of values for 'z', 't0', 'x0', 'x1', 'c', their respective
        errors, the fit chi-squared, number of DOF, and SNCosmo exit message.
    """

    out_data = []

    # Try fitting the light-curve
    try:
        result, fitted_model = sncosmo.fit_lc(
            data=data, model=model, vparam_names=vparam_names, **kwargs)

    # If the fit fails fill out_data with place holder values (NANs and zeros)
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
            out_data.extend(np.full(10, np.NAN).tolist())
            out_data.extend((0, 0))
            out_data.append(str(e).replace('\n', ' '))

    # Finish populating out_data with fit results
    else:
        # Append fit parameter values
        for param in ['z', 't0', 'x0', 'x1', 'c']:
            i = result.param_names.index(param)
            value = result.parameters[i]
            out_data.append(value)

            if param == 't0':
                tmax = value

        # Append fit parameter errors
        for param in ['z', 't0', 'x0', 'x1', 'c']:
            out_data.append(result.errors.get(param, 0))

        # Determine peak magnitude and decline rate
        b_max = fitted_model.source_peakabsmag('bessellb', 'AB')
        delta_15 = np.nan

        out_data.append(result.chisq)
        out_data.append(result.ndof)
        out_data.append(b_max)
        out_data.append(delta_15)
        out_data.extend(_count_pre_and_post_max(data['time'], tmax))
        out_data.append(result.message.replace('\n', ' '))

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
    out_table = _create_empty_summary_table()
    out_table.meta = kwargs

    # Run fit for each target
    for input_table in inputs:
        # Protect against mutating input args
        model_this = deepcopy(model)
        kwargs_this = deepcopy(kwargs)

        # Set redshift in model for 4 param fits
        if num_params == 4:
            vparam_names = ['t0', 'x0', 'x1', 'c']
            z = input_table.meta.get('redshift', False)
            if z < 0 or not z:
                continue

            model_this.set(z=z)

        else:
            vparam_names = ['z', 't0', 'x0', 'x1', 'c']

        fit_results = fit_lc(
            data=input_table,
            model=model_this,
            vparam_names=vparam_names,
            **kwargs_this)

        new_row = [input_table.meta['cid'], len(input_table)]
        new_row.extend(fit_results)
        out_table.add_row(new_row)
        out_table.write(out_path, overwrite=True)
