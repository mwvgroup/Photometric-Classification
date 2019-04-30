#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Provides a function based interface for iteratively fitting light-curves for
a single model / survey.

The fit_lc function acts as a wrapper for sncosmo.fit_lc that returns fit
results as a list.

The fit_n_params function iteratively fits a collection of light-curves using a
4 or 5 parameter fit.
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


def _simplify_t0_bounds(bounds_dict, test_time):
    """Simplify user specified bounds into a form compatable with SNCsomo

    If bounds_dict['t0'] is a list of bounds for multiple observing seasons,
    return a copy of bounds_dict where 't0' is replaced with the bounds
    encompassing the value test_time. If no such bound exists, combine the list
    of bounds into a boundary spanning the entire survey.

    Args:
        bounds_dict (dict):
        test_time(dict)

    Returns:
         A Dictionary
    """

    if 't0' not in bounds_dict:
        return bounds_dict

    bounds_dict = bounds_dict.copy()
    if not isinstance(bounds_dict['t0'][0], (list, tuple)):
        return bounds_dict

    for low_bound, up_bound in bounds_dict['t0']:
        if low_bound <= test_time <= up_bound:
            bounds_dict['t0'] = [low_bound, up_bound]
            return bounds_dict

    bounds_dict['t0'] = [bounds_dict['t0'][0][0], bounds_dict['t0'][-1][-1]]
    return bounds_dict


def _fit_with_nesting(data, model, nest, vparam_names, **kwargs):
    kwargs['bounds']['t0'] = (min(data['time']) - 20, max(data['time']))

    if nest:
        nest_result, _ = sncosmo.nest_lc(
            data, model, vparam_names,
            bounds=kwargs['bounds'],
            verbose=True,
            maxiter=3000
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

    return result, fitted_model


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
        result, fitted_model = _fit_with_nesting(
            data, model, nest, vparam_names, **kwargs)

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
    """Iteratively fit light curves with a 4 or 5 parameter Salt2-like model

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
