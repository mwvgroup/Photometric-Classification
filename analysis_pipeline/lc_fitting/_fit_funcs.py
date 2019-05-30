#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module provides wrappers for SNCosmo fitting functions.

The nest_lc function acts as a wrapper for sncosmo.nest_lc but assumes
different default input parameters.

The fit_lc function acts as a wrapper for sncosmo.fit_lc that returns fit
results as a list.
"""

from copy import deepcopy

import numpy as np
import sncosmo
from astropy.table import Table
from sncosmo.fitting import DataQualityError


def create_empty_summary_table():
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


def nest_lc(data, model, vparam_names, **kwargs):
    """Set initial model params using nested sampling

    If not specified, bounds on t0 assume the first observation is less than
    15 days after peak and before the last data point.

    Arguments are protected against mutation except for ``bounds`` which may be
    modified as outlined above.

    Args:
        data             (Table): Table of light curve data for SNCosmo
        model            (Model): SNCosmo model
        vparam_names (list[str]): List of parameters to vary in the model
        bounds            (dict): Boundaries on fit parameters
        verbose           (bool): Whether to display progress (Default: True)
        maxiter            (int): Maximum sampling iterations (Default: 5000)
        maxcall            (int): Maximum function calls (Default: 20000)
        method             (str): Nested sampling method (Default: 'multi')
    """

    # Assume first observation < 15 days after peak and before last data point
    t0_start = min(data['time']) - 15
    t0_end = max(data['time'])
    kwargs['bounds']['t0'] = kwargs['bounds'].get('t0', (t0_start, t0_end))

    # Protect mutable arguments
    kwargs = deepcopy(kwargs)
    model = deepcopy(model)

    # Set other default values
    kwargs['verbose'] = kwargs.get('verbose', True)
    kwargs['maxiter'] = kwargs.get('maxiter', 10000)
    kwargs['maxcall'] = kwargs.get('maxcall', 20000)
    kwargs['method'] = kwargs.get('method', 'multi')

    # Set initial parameters in model
    nest_result, _ = sncosmo.nest_lc(data, model, vparam_names, **kwargs)
    model_args = dict()
    for vp, p in zip(nest_result.vparam_names, nest_result.parameters):
        model_args[vp] = p

    model.set(**model_args)
    return model


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
    message. Mutable arguments are not guaranteed to go unchanged.

    Args:
        Any arguments for sncosmo.fit_lc

    Returns:
        A list of values for the target id, number of photometric points,
        'z', 't0', 'x0', 'x1', 'c', their respective errors,
        the fit chi-squared, number of DOF, b_max, delta_15, t_min, t_max,
        the number of points pre / post max, and the SNCosmo exit message.
    """

    out_data = [data.meta['cid'], len(data)]

    try:
        if len(data) == 0:
            raise RuntimeError('Empty table passed to fit_lc.')

        result, fitted_model = sncosmo.fit_lc(
            data, deepcopy(model), vparam_names, **deepcopy(kwargs))

    except KeyboardInterrupt:
        raise

    # If the fit fails fill out_data with place holder values (NANs and zeros)
    except (DataQualityError, RuntimeError, ValueError) as e:
        out_data.extend(np.full(16, np.NAN).tolist())
        out_data.extend((0, 0, str(e).replace('\n', ' ')))
        if 'z' not in vparam_names:
            out_data[2] = data.meta['redshift']
            out_data[7] = data.meta.get('redshift_err', np.NAN)

    else:
        params = {p: v for p, v in zip(result.param_names, result.parameters)}

        # Add fit results to out_data
        out_data.extend(result.parameters)
        for param in result.param_names:
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
