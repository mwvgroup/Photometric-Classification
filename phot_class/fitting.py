#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``fitting`` module runs a series of fits on individual light-curves and
tabulates the results. This includes the ability to fit each bandpass
independently (like in SiFTO) or to fit restframe blue/red band passes os
separate, collective sets.

Function Documentation
----------------------
"""

from copy import deepcopy

import numpy as np
import sncosmo
from astropy.table import Table
from matplotlib import pyplot

from . import utils

DUST = sncosmo.F99Dust()


def create_empty_table(parameters, **kwargs):
    """Create an empty table for storing fit results

    Columns:
        - obj_id
        - band
        - source
        - pre_max
        - post_max
        - num_params
        - *parameters
        - *parameters + _err
        - chisq
        - ndof
        - b_max
        - delta_15
        - message

    Args:
        parameters (iter): List of parameter names to add columns for
        Any arguments to pass ``astropy.Table``

    Returns:
        A masked astropy Table
    """

    # Specify column names
    names = ['obj_id', 'band', 'source', 'pre_max', 'post_max', 'vparams']
    names += list(parameters) + [param + '_err' for param in parameters]
    names += ['chisq', 'ndof', 'b_max', 'delta_15', 'message']

    # Specify column data types
    dtype = ['U20', 'U100', 'U100', int, int, 'U100']
    dtype += [float for _ in range(2 * len(parameters))]
    dtype += [float, float, float, float, 'U10000']

    # Unless otherwise specified, we default to returning a masked table
    kwargs = deepcopy(kwargs)
    kwargs.setdefault('masked', True)
    return Table(names=names, dtype=dtype, **kwargs)


def fit_results_to_dict(data, obj_id, band_set, results, fitted_model):
    """Format sncosmo fit results so they can be appended to an astropy table

    See the ``create_empty_table`` function for information on the assumed
    table format.

    Args:
        data         (Table): The data used in the fit
        obj_id         (str): The id of the object that was fit
        band_set       (str): The name of the band set ('all', 'blue', 'red')
        results     (Result): Fitting results returned by ``sncosmo``
        fitted_model (Model): A fitted ``sncosmo`` model

    Returns:
        Fit results as a dictionary
    """

    new_row = {
        'obj_id': obj_id,
        'band': band_set,
        'source': fitted_model.source.name,
        'vparams': ','.join(results.vparam_names)
    }

    # Determine number of points pre and post maximum
    t0 = results.parameters[results.param_names.index('t0')]
    new_row['pre_max'] = sum(data['time'] < t0)
    new_row['post_max'] = sum(data['time'] >= t0)

    # Add parameters and their errors
    params = {p: v for p, v in zip(results.param_names, results.parameters)}
    new_row.update(params)
    for param, error in results.errors.items():
        new_row[param + '_err'] = error

    # Calc chi-squared
    chisq, ndof = utils.calc_model_chisq(data, results, fitted_model)
    new_row['chisq'] = np.round(chisq, 2)
    new_row['ndof'] = ndof

    # Determine peak magnitude and decline rate
    b_max = fitted_model.source_peakabsmag('bessellb', 'ab')
    peak_phase = fitted_model.source.peakphase('bessellb')
    b_0 = fitted_model.source.bandmag('bessellb', 'ab', peak_phase)
    b_15 = fitted_model.source.bandmag('bessellb', 'ab', peak_phase + 15)
    delta_15 = b_15 - b_0
    new_row['b_max'] = np.round(b_max, 2)
    new_row['delta_15'] = np.round(delta_15, 3)

    # Add fitting exit status message. Not all fitting routines include
    # this attribute, so we assign a default value of 'NONE'.
    message = getattr(results, 'message', 'NONE')
    new_row['message'] = message
    return new_row


def _plot_lc(data, result, fitted_model, show=True):
    """Plot fit results

    Args:
        data         (Table): The data used in the fit
        result      (Result): The fit results
        fitted_model (Model): Model with params set to fitted values
    """

    fig = sncosmo.plot_lc(data, fitted_model, errors=result.errors)
    xs, d = utils.calc_model_chisq(data, result, fitted_model)
    print(f'chisq / ndof = {xs} / {d} = {xs / d}', flush=True)
    if show:
        pyplot.show()

    return fig


def _create_fit_data_iter(priors_hs, priors_bg, kwargs_hs, kwargs_bg):
    """Create an iterable of data used to run light-curve fits

    Args:
        priors_hs  (dict): Priors to use when fitting hsiao
        priors_bg  (dict): Priors to use when fitting sn91bg
        kwargs_hs  (dict): Kwargs to pass ``fit_func`` when fitting salt2
        kwargs_bg  (dict): Kwargs to pass ``fit_func`` when fitting sn91bg

    Returns:
        - An iterable of models, vparams, priors, and kwargs for fitting
        - A table to store output data
    """

    # Set default kwargs and protect against mutation
    priors_hs = deepcopy(priors_hs) or dict()
    priors_bg = deepcopy(priors_bg) or dict()
    kwargs_hs = deepcopy(kwargs_hs) or dict()
    kwargs_bg = deepcopy(kwargs_bg) or dict()

    # Define models for normal and 91bg SNe with host galaxy dust
    dust_kw = dict(effects=[DUST], effect_names=['mw'], effect_frames=['obs'])
    bg_source = sncosmo.get_source('sn91bg', version='hsiao_phase')
    sn91bg = sncosmo.Model(bg_source, **dust_kw)
    hsiao = sncosmo.Model('hsiao_x1', **dust_kw)

    # Determine what parameters to vary for each model
    # Hsiao does not have a c parameter. We don't vary mwebv
    vparams = {'z', 't0', 'amplitude', 'x1', 'c'}
    out_data = create_empty_table(vparams.union({'mwebv'}))
    if 'z' in priors_bg and 'z' in priors_hs:
        vparams -= {'z'}

    hsiao_vparams = set(hsiao.param_names).intersection(vparams)
    sn91bg_vparams = set(sn91bg.param_names).intersection(vparams)

    # Create iterators over the data we need to fit
    model_args = zip(
        (hsiao, sn91bg),  # The models
        (hsiao_vparams, sn91bg_vparams),  # The parameters to vary
        (priors_hs, priors_bg),  # The priors
        (kwargs_hs, kwargs_bg)  # The fitting kwargs
    )

    return model_args, out_data


def run_band_fits(
        obj_id, data, fit_func,
        priors_hs=None, priors_bg=None,
        kwargs_hs=None, kwargs_bg=None,
        show_plots=False):
    """Run light curve fits on a given target using the Hsiao and 91bg model

    Fits are run using both the ``hsiao_x1`` and ``sn91bg`` models for all
    available bands and then for each band individually.

    Varied parameters include ``z``, ``t0``, ``amplitude``, ``x1``, and ``c``.
    If the ``z`` is specified in the priors for both models, it is not varied
    in any fit. The parameters ``t0`` and ``z`` are not varied in the
    individual band fits.

    Args:
        obj_id      (str): Id of the object being fitted
        data      (Table): Table of photometric data
        fit_func   (func): Function to use to run fits (eg. ``fit_funcs.fit_lc``)
        priors_hs  (dict): Priors to use when fitting hsiao
        priors_bg  (dict): Priors to use when fitting sn91bg
        kwargs_hs  (dict): Kwargs to pass ``fit_func`` when fitting salt2
        kwargs_bg  (dict): Kwargs to pass ``fit_func`` when fitting sn91bg
        show_plots (bool): Plot and display each individual fit

    Returns:
       A table with results each model / dataset combination
    """

    model_args, out_data = _create_fit_data_iter(
        priors_hs, priors_bg, kwargs_hs, kwargs_bg)

    # Tabulate fit results for each band
    for model, vparams, prior, kwarg in model_args:
        model.update(prior)
        kwarg['bounds'] = \
            {p: v for p, v in kwarg.get('bounds', {}).items() if p in vparams}

        # Fit data in all bands
        result_all, fit_all = fit_func(data, model, vparams, **kwarg)
        new_row = fit_results_to_dict(data, obj_id, 'all', result_all,
                                      fit_all)
        out_data.add_row(new_row)

        if show_plots:
            _plot_lc(data, result_all, fit_all)

        # Fix t0 and z during individual band fits
        band_vparams = deepcopy(vparams) - {'t0', 'z'}
        kwarg['bounds'].pop('t0', None)
        kwarg['bounds'].pop('z', None)

        # Fit data in individual bands
        data = data.group_by('band')
        for band_name, band_data in zip(data.groups.keys['band'], data.groups):
            # Using amplitude from all data fit as initial guess works better
            kwarg['guess_amplitude'] = False
            result, fit = fit_func(band_data, fit_all, band_vparams, **kwarg)
            new_row = fit_results_to_dict(band_data, obj_id, band_name,
                                          result, fit)
            out_data.add_row(new_row)

            if show_plots:
                _plot_lc(band_data, result, fit)

    return out_data


def run_collective_fits(
        obj_id, data, fit_func,
        band_names, lambda_eff,
        priors_hs=None, priors_bg=None,
        kwargs_hs=None, kwargs_bg=None,
        show_plots=False):
    """Run light curve fits on a given target using the Hsiao and 91bg model

    Args:
        obj_id      (str): Id of the object being fitted
        data      (Table): Table of photometric data
        fit_func   (func): Function to use to run fits (eg. ``fit_funcs.fit_lc``)
        band_names (list): Name of bands included in ``data_iter``
        lambda_eff (list): Effective wavelength for bands in ``band_names``
        priors_hs  (dict): Priors to use when fitting hsiao
        priors_bg  (dict): Priors to use when fitting sn91bg
        kwargs_hs  (dict): Kwargs to pass ``fit_func`` when fitting salt2
        kwargs_bg  (dict): Kwargs to pass ``fit_func`` when fitting sn91bg
        show_plots (bool): Plot and display each individual fit

    Returns:
       A table with results each model / dataset combination
    """

    model_args, out_data = _create_fit_data_iter(
        priors_hs, priors_bg, kwargs_hs, kwargs_bg)

    # Tabulate fit results for each band
    for model, vparams, prior, kwarg in model_args:
        model.update(prior)
        kwarg['bounds'] = \
            {p: v for p, v in kwarg.get('bounds', {}).items() if p in vparams}

        # Fit data in all bands
        result_all, fit_all = fit_func(data, model, vparams, **kwarg)
        new_row = fit_results_to_dict(data, obj_id, 'all', result_all,
                                      fit_all)
        out_data.add_row(new_row)

        if show_plots:
            _plot_lc(data, result_all, fit_all)

        # Fix t0 during individual band fits
        band_vparams = deepcopy(vparams) - {'t0', 'z'}
        kwarg['bounds'].pop('t0', None)
        kwarg['bounds'].pop('z', None)

        # Get red and blue data
        z = fit_all.parameters[fit_all.param_names.index('z')]
        blue_data, red_data = utils.split_data(data, band_names, lambda_eff, z)

        for band_name, band_data in zip(('blue', 'red'),
                                        (blue_data, red_data)):
            # Using amplitude from all data fit as initial guess works better
            kwarg['guess_amplitude'] = False
            result, fit = fit_func(band_data, fit_all, band_vparams, **kwarg)
            new_row = fit_results_to_dict(band_data, obj_id, band_name,
                                          result, fit)
            out_data.add_row(new_row)

            if show_plots:
                _plot_lc(band_data, result, fit)

    return out_data
