#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module provides functions for fitting light curves with SNCosmo."""

from copy import deepcopy

import numpy as np
import sncosmo
from astropy.table import Table
from sncosmo.fitting import DataQualityError


def create_results_table(param_names):
    """Create an empty table for storing fit results

    Columns:
         obj_id, num_points, *param_names, *param_names + _err, chi, dof,
         tmin, tmax, pre_max, post_max, message

    Args:
        param_names (list[str]): List of parameter names

    Returns:
        An astropy Table
    """

    # Specify column names
    names = ['obj_id', 'num_points']
    names.extend(param_names)
    names.extend((name + '_err' for name in param_names))
    names.extend(
        ('chi', 'dof', 'b_max',
         'delta_15', 'tmin', 'tmax',
         'pre_max', 'post_max', 'message'))

    # Specify column data types
    dtype = ['U20', int]
    dtype.extend((float for _ in param_names))
    dtype.extend((float for _ in param_names))
    dtype.extend((float, float, float, float, float, float, int, int, 'U1000'))

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


def calc_chisq(data, model):
    """
    Calculate the chi-squared for a given data table and model
    
    Args:
        data  (Table): An SNCosmo input table
        model (Model): An SNCosmo Model
        
    Returns:
        The un-normalized chi-squared
        The number of data points used in the calculation
    """

    while True:
        try:
            # Model flux and keep only non-zero values
            model_flux = model.bandflux(data['band'], data['time'])
            data = data[model_flux > 0]
            model_flux = model_flux[model_flux > 0]

            residuals = model_flux - data['flux']
            chisq = np.sum((residuals / data['fluxerr']) ** 2)
            return chisq, len(data)

        except ValueError as err:
            # Remove bands that are out of model range
            data = data[data['band'] != err.args[0].split()[1][1:-1]]


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

    kwargs = deepcopy(kwargs)
    kwargs['warn'] = kwargs.get('warn', False)

    out_data = [data.meta['obj_id'], len(data)]

    try:
        if len(data) == 0:
            raise RuntimeError('Empty table passed to fit_lc.')

        result, fitted_model = sncosmo.fit_lc(
            data, deepcopy(model), vparam_names, **kwargs)

    except KeyboardInterrupt:
        raise

    # If the fit fails fill out_data with place holder values (NANs and zeros)
    except (DataQualityError, RuntimeError, ValueError) as e:
        num_cols = 2 * len(model.param_names) + 6
        out_data.extend(np.full(num_cols, np.NAN).tolist())
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

        # Calc chi-squared
        man_chi, num_points = calc_chisq(data, fitted_model)
        man_dof = num_points - len(vparam_names)

        # Add remaining data
        out_data.append(man_chi)
        out_data.append(man_dof)
        out_data.append(b_max)
        out_data.append(delta_15)
        out_data.extend(_count_pre_and_post_max(data['time'], params['t0']))
        out_data.append(result.message.replace('\n', ' '))

    return out_data
