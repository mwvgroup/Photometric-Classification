#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``classification`` module fits light curves, tabulates results, and
determines the corresponding fitting coordinates.
"""

from copy import deepcopy

import numpy as np
import sncosmo
from astropy.table import Table, vstack
from matplotlib import pyplot

from . import utils


def create_empty_table(parameters, **kwargs):
    """Create an empty table for storing fit results

    Columns:
        - obj_id
        - band
        - source
        - pre_max
        - post_max
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
    names = ['obj_id', 'band', 'source', 'pre_max', 'post_max']
    names += list(parameters) + [param + '_err' for param in parameters]
    names += ['chisq', 'ndof', 'b_max', 'delta_15', 'message']

    # Specify column data types
    dtype = ['U20', 'U100', 'U100', int, int]
    dtype += [float for _ in range(2 * len(parameters))]
    dtype += [float, float, float, float, 'U10000']

    # Unless otherwise specified, we default to returning a masked table
    kwargs = deepcopy(kwargs)
    kwargs.setdefault('masked', True)
    return Table(names=names, dtype=dtype, **kwargs)


def _fit_results_to_dict(data, obj_id, band_set, results, fitted_model):
    """Format sncosmo fit results so they can be appended to an astropy table

    See the ``create_empty_table`` function for information on the assumed
    table format.

    Args:
        data         (Table): The data used in the fit
        band_set       (str): The name of the band set ('all', 'blue', 'red')
        results     (Result): Fitting results returned by ``sncosmo``
        fitted_model (Model): A fitted ``sncosmo`` model

    Returns:
        Fit results as a list formatted for addition to an astropy table
    """

    new_row = {
        'obj_id': obj_id,
        'band': band_set,
        'source': fitted_model.source.name
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


def _raise_unspecified_params(fixed_params, *priors):
    """Raise RuntimeError if any fixed parameters are not set in a prior

    Args:
        fixed_params (iter): List of parameter names to assume as fixed
        *priors      (dict): Dictionary of values from a prior
    """

    for prior in priors:
        for p in fixed_params:
            if p not in prior:
                raise RuntimeError(
                    f'Encountered unspecified value for {p} in prior: {prior}')


def _plot_lc(data, result, fitted_model, show=True):
    """Plot fit results

    Args:
        data         (Table): The data used in the fit
        result      (Result): The fit results
        fitted_model (Model): Model with params set to fitted values
    """

    fig = sncosmo.plot_lc(data, fitted_model, color='blue')
    xs, d = utils.calc_model_chisq(data, result, fitted_model)
    print(f'chisq / ndof = {xs} / {d} = {xs / d}', flush=True)
    if show:
        pyplot.show()

    return fig


def run_band_fits(
        obj_id, data, vparams, fit_func,
        priors_hs=None, priors_bg=None,
        kwargs_hs=None, kwargs_bg=None,
        show_plots=False):
    """Run light curve fits on a given target using the Hsiao and 91bg model

    Fits are run using both models for all available bands and then for each
    band individually. All parameters specified in ``vparams`` are allowed to
    vary during all fits except for `t0` which is fixed to the value
    determined by fitting all available bands.

    Args:
        obj_id      (str): Id of the object being fitted
        data      (Table): Table of photometric data
        vparams    (iter): Iterable of param names to fit in any of the models
        fit_func   (func): Function to use to run fits (eg. ``fit_funcs.fit_lc``)
        priors_hs  (dict): Priors to use when fitting hsiao
        priors_bg  (dict): Priors to use when fitting sn91bg
        kwargs_hs  (dict): Kwargs to pass ``fit_func`` when fitting salt2
        kwargs_bg  (dict): Kwargs to pass ``fit_func`` when fitting sn91bg
        show_plots (bool): Plot and display each individual fit

    Returns:
       Fit results and the fitted model for each model / data combination
    """

    # Set default kwargs and protect against mutation
    priors_hs = deepcopy(priors_hs) or dict()
    priors_bg = deepcopy(priors_bg) or dict()
    kwargs_hs = deepcopy(kwargs_hs) or dict()
    kwargs_bg = deepcopy(kwargs_bg) or dict()

    # Define models for normal and 91bg SNe
    hsiao = sncosmo.Model('hsiao_x1')
    sn91bg = sncosmo.Model(sncosmo.get_source('sn91bg', version='hsiao_phase'))

    # Make sure any model parameters that are not being varied are
    # specified in the priors
    hsiao_params = set(hsiao.param_names)
    hsiao_vparams = hsiao_params.intersection(vparams)
    hsiao_fixed_params = hsiao_params - hsiao_vparams
    _raise_unspecified_params(hsiao_fixed_params, priors_hs)

    sn91bg_params = set(sn91bg.param_names)
    sn91bg_vparams = sn91bg_params.intersection(vparams)
    sn91bg_fixed_params = sn91bg_params - sn91bg_vparams
    _raise_unspecified_params(sn91bg_fixed_params, priors_bg)

    # Create iterators over the data we need to fit
    model_args = zip(
        (hsiao, sn91bg),  # The models
        (hsiao_vparams, sn91bg_vparams),  # The parameters to vary
        (priors_hs, priors_bg),  # The priors
        (kwargs_hs, kwargs_bg)  # The fitting kwargs
    )

    # Tabulate fit results for each band
    all_param_names = set(hsiao.param_names).union(sn91bg.param_names)
    out_data = create_empty_table(all_param_names)
    for model, vparams, prior, kwarg in model_args:
        model.update(prior)

        # Fit data in all bands
        result_all, fit_all = fit_func(data, model, vparams, **kwarg)
        new_row = _fit_results_to_dict(data, obj_id, 'all', result_all, fit_all)
        out_data.add_row(new_row)

        if show_plots:
            _plot_lc(data, result_all, fit_all)

        # Fix t0 during individual band fits
        band_vparams = deepcopy(vparams)
        if 't0' in band_vparams:
            band_vparams -= {'t0'}

        # Fit data in individual bands
        data = data.group_by('band')
        for band_name, band_data in zip(data.groups.keys['band'], data.groups):
            result, fit = fit_func(band_data, fit_all, band_vparams, guess_amplitude=False, **kwarg)
            new_row = _fit_results_to_dict(band_data, obj_id, band_name, result, fit)
            out_data.add_row(new_row)

            if show_plots:
                _plot_lc(band_data, result, fit)

    return out_data


def tabulate_fit_results(
        data_iter, band_names, lambda_eff, fit_func, vparams,
        config=None, out_path=None):
    """Tabulate fit results for a collection of data tables

    Args:
        data_iter  (iter): Iterable of photometric data for different SN
        band_names (list): Name of bands included in ``data_iter``
        lambda_eff (list): Effective wavelength for bands in ``band_names``
        vparams    (list): Name of parameters to vary
        fit_func   (func): Function to use to run fits
        config     (dict): Specifies priors / kwargs for fitting each model
        out_path    (str): Optionally cache progressive results to file

    Returns:
       An astropy table with fit results
    """

    # Set default kwargs
    config = deepcopy(config) or dict()

    # Add meta_data to output table meta data
    out_table = Table(names=['obj_id', 'message'], dtype=['U20', 'U10000'])
    out_table.meta['band_names'] = band_names
    out_table.meta['lambda_eff'] = lambda_eff
    out_table.meta['fit_func'] = fit_func.__name__
    out_table.meta['vparams'] = vparams
    out_table.meta['out_path'] = str(out_path)

    for data in data_iter:
        # Get fitting priors and kwargs
        obj_id = data.meta['obj_id']
        salt2_prior, salt2_kwargs, sn91bg_prior, sn91bg_kwargs = \
            utils.parse_config_dict(obj_id, config)

        try:
            fit_results = run_band_fits(
                obj_id=obj_id,
                data=data,
                vparams=vparams,
                fit_func=fit_func,
                priors_hs=salt2_prior,
                priors_bg=sn91bg_prior,
                kwargs_hs=salt2_kwargs,
                kwargs_bg=sn91bg_kwargs
            )

            out_table = vstack([out_table, fit_results])

        except KeyboardInterrupt:
            raise

        except Exception as e:
            out_table.add_row({
                'obj_id': obj_id,
                'message': str(e).replace('\n', '')
            })

        if out_path:
            out_table.write(out_path)

    return out_table


def classify_targets(
        fits_table, band_names=None, lambda_eff=None, out_path=None):
    """Tabulate fitting coordinates for SNe based on their fit results

    Assumed columns in ``fits_table`` include 'obj_id', 'source', 'band',
    'chisq', and 'ndof'.  Values in the 'source' column should be either
    'hsiao_x1' or 'sn91bg'.

    If ``band_names`` or ``lambda_eff`` are not given, they are taken from
    ``fits_table.meta``.

    Args:
        fits_table (Table): A table of fit results
        band_names  (list): List of band names used when fitting
        lambda_eff  (list): The effective wavelength of each band in angstroms
        out_path     (str): Optionally write results to file

    Returns:
        An astropy table of fitting coordinates
    """

    if band_names is None or lambda_eff is None:
        band_names = fits_table.meta['band_names']
        lambda_eff = fits_table.meta['lambda_eff']

    # Todo: this needs to take redshift into account
    blue_bands, red_bands = utils.split_bands(band_names, lambda_eff)

    # Convert input table to a DataFrame so we can leverage multi-indexing
    fits_df = fits_table.to_pandas()
    fits_df.set_index(['obj_id', 'source'], inplace=True)
    fits_df.dropna(subset=['band'])

    out_table = Table(names=['obj_id', 'x', 'y'], dtype=['U100', float, float])
    for obj_id in fits_df.index.unique(level='obj_id'):

        try:
            hsiao_data = fits_df.loc[obj_id, 'hsiao_x1']
            hsiao_blue = hsiao_data[hsiao_data['band'].isin(blue_bands)]
            hsiao_red = hsiao_data[hsiao_data['band'].isin(red_bands)]
            hsiao_blue_chisq = hsiao_blue['chisq'].sum() / hsiao_blue['ndof'].sum()
            hsiao_red_chisq = hsiao_red['chisq'].sum() / hsiao_red['ndof'].sum()

            sn91bg_data = fits_df.loc[obj_id, 'sn91bg']
            sn91bg_blue = sn91bg_data[sn91bg_data['band'].isin(blue_bands)]
            sn91bg_red = sn91bg_data[sn91bg_data['band'].isin(red_bands)]
            sn91bg_blue_chisq = sn91bg_blue['chisq'].sum() / sn91bg_blue['ndof'].sum()
            sn91bg_red_chisq = sn91bg_red['chisq'].sum() / sn91bg_red['ndof'].sum()

        except KeyError:
            continue

        else:
            x = hsiao_blue_chisq - sn91bg_blue_chisq
            y = hsiao_red_chisq - sn91bg_red_chisq
            out_table.add_row([obj_id, x, y])

            if out_path:
                out_table.write(out_path)

    return out_table
