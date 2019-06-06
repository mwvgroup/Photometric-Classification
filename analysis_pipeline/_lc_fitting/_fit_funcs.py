#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module provides wrappers for SNCosmo fitting functions.

The nest_lc function acts as a wrapper for sncosmo.nest_lc but assumes
different default input parameters and only returns a model object

The get_sampled_model function is the same as the nest_lc functions, but
uses cached values when available.

The fit_lc function acts as a wrapper for sncosmo.fit_lc that returns fit
results as a list.
"""

from copy import deepcopy
from pathlib import Path

import numpy as np
import sncosmo
from astropy.table import Table
from sncosmo.fitting import DataQualityError

PRIOR_DIR = Path(__file__).resolve().parent.parent / 'priors'
PRIOR_DIR.mkdir(exist_ok=True)
PRIORS = dict()  # For lazy loading priors


def create_empty_summary_table():
    """Create an empty table for storing fit results.

    Columns:
         obj_id, num_points, z, t0, x0, x1, z_err, t0_err, x0_err,
         x1_err, c_err, chi, dof, tmin, tmax, pre_max, post_max, message

    Returns:
        An astropy Table
    """

    param_names = (
        'z', 't0', 'x0', 'x1', 'c',
        'z_err', 't0_err', 'x0_err', 'x1_err', 'c_err'
    )

    # Specify column names
    names = ['obj_id', 'num_points']
    names.extend(param_names)
    names.extend(
        ('chi', 'dof', 'b_max',
         'delta_15', 'tmin', 'tmax',
         'pre_max', 'post_max', 'message'))

    # Specify colomn data types
    dtype = ['U20', int]
    dtype.extend((float for _ in param_names))
    dtype.extend((float, float, float, float, float, float, int, int, 'U1000'))

    return Table(names=names, dtype=dtype)


def create_empty_priors_table(model):
    """Create an empty table for storing fit priors.

    Columns:
         obj_id, *<model.param_names>, *<model.param_names>_min,
         *<model.param_names>_max

    Args:
        model (Model): An SNCosmo model

    Returns:
        An astropy Table
    """

    col_names = ['obj_id']
    dtype = ['U100']
    for param in model.param_names:
        col_names.extend((param, param + '_min', param + '_max'))
        dtype.extend((float, float, float))

    priors_table = Table(names=col_names, dtype=dtype)
    return priors_table


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
    # Intentionally mutate kwargs
    t0_start = min(data['time']) - 15
    t0_end = max(data['time'])
    kwargs['bounds']['t0'] = kwargs['bounds'].get('t0', (t0_start, t0_end))

    # Protect against **further** argument mutation
    kwargs = deepcopy(kwargs)
    model = deepcopy(model)

    # Set kwargs specific to sncosmo.nest_lc and get model with nested values
    kwargs['verbose'] = kwargs.get('verbose', True)
    kwargs['maxiter'] = kwargs.get('maxiter', 10000)
    kwargs['maxcall'] = kwargs.get('maxcall', 20000)
    _, nest_model = sncosmo.nest_lc(data, model, vparam_names, **kwargs)
    return nest_model


def get_sampled_model(survey_name, data, model, vparam_names, **kwargs):
    """Set initial model params using cached values

    If cached values are not available, determine them using nested sampling
    and save them to file.

    Args:
        survey_name        (str): The name of the survey being fit
        data             (Table): Table of light curve data for SNCosmo
        model            (Model): SNCosmo model
        vparam_names (list[str]): List of parameters to vary in the model
        bounds            (dict): Boundaries on fit parameters
        Any other arguments for analysis_pipeline.lc_fitting.nest_lc
    """

    # Get path of priors file
    model_name = f'{model.source.name}_{model.source.version}'
    file_name = f'{survey_name.lower()}_{model_name}.ecsv'
    file_path = PRIOR_DIR / file_name

    if file_path.exists():
        # Lazy load prior for object represented by ``data``
        priors_table = PRIORS.setdefault(file_path, Table.read(file_path))
        prior = priors_table[priors_table['obj_id'] == data.meta['obj_id']]

        if prior:
            sampled_model = deepcopy(model)
            for param in vparam_names:
                sampled_model.update({param: prior[param]})
                kwargs['bounds'][param] = \
                    (prior[f'{param}_min'], prior[f'{param}_max'])

            return sampled_model

    else:
        PRIORS[file_path] = create_empty_priors_table(model)

    # Calculate prior values
    try:
        sampled_model = nest_lc(data, model, vparam_names, **kwargs)

    except ValueError:
        sampled_model = deepcopy(model)
        param_dict = {p: np.median(kwargs['bounds'][p]) for p in sampled_model.param_names}
        sampled_model.update(param_dict)

    new_row = [data.meta['obj_id']]
    for param_name, param_val in zip(sampled_model.param_names,
                                     sampled_model.parameters):
        new_row.append(param_val)
        new_row.extend(kwargs['bounds'][param_name])

    PRIORS[file_path].add_row(new_row)
    PRIORS[file_path].write(file_path, overwrite=True)

    return sampled_model


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

    out_data = [data.meta['obj_id'], len(data)]

    try:
        if len(data) == 0:
            raise RuntimeError('Empty table passed to fit_lc.')

        result, fitted_model = sncosmo.fit_lc(
            data, deepcopy(model), vparam_names, **deepcopy(kwargs))

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

        # Add remaining data
        out_data.append(result.chisq)
        out_data.append(result.ndof)
        out_data.append(b_max)
        out_data.append(delta_15)
        out_data.extend(_count_pre_and_post_max(data['time'], params['t0']))
        out_data.append(result.message.replace('\n', ' '))

    return out_data