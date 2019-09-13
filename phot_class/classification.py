#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``fitting`` module fits light curves, tabulates results, and
determines the corresponding fitting coordinates.
"""

from copy import deepcopy

import numpy as np
import sncosmo
from astropy.table import Table, vstack
from matplotlib import pyplot

from . import utils


def create_empty_table(**kwargs):
    """Create an empty table for storing fit results

    Columns:
        - obj_id
        - band_set
        - source
        - pre_max
        - post_max
        - z
        - t0
        - x0
        - x1
        - c
        - z_err
        - t0_err
        - x0_err
        - x1_err
        - c_err
        - chisq
        - ndof
        - b_max
        - delta_15
        - message

    Args:
        Any arguments to pass ``astropy.Table``

    Returns:
        A masked astropy Table
    """

    parameters = ['z', 't0', 'x0', 'x1', 'c']

    # Specify column names
    names = ['obj_id', 'band_set', 'source', 'pre_max', 'post_max']
    names += parameters
    names += [param + '_err' for param in parameters]
    names += ['chisq', 'ndof', 'b_max', 'delta_15', 'message']

    # Specify column data types
    dtype = ['U20', 'U4', 'U100', int, int]
    dtype += [float for _ in range(2 * len(parameters))]
    dtype += [float, float, float, float, 'U10000']

    # Unless otherwise specified, we default to returning a masked table
    kwargs = deepcopy(kwargs)
    kwargs.setdefault('masked', True)
    return Table(names=names, dtype=dtype, **kwargs)


def fit_results_to_table_row(data, band_set, results, fitted_model):
    """Format sncosmo fit results so they can be appended to an astropy table

    See the ``create_empty_table`` function for information on the assumed
    table format. ``fitted_model`` is assumed to have the parameters 'z', 't0',
    'x0', 'x1', and 'c' (in order).

    Args:
        data         (Table): The data used in the fit
        band_set       (str): The name of the band set ('all', 'blue', 'red')
        results     (Result): Fitting results returned by ``sncosmo``
        fitted_model (Model): A fitted ``sncosmo`` model

    Returns:
        Fit results as a list formatted for addition to an astropy table
    """

    new_row = [data.meta['obj_id'], band_set, fitted_model.source.name]

    # Determine number of points pre and post maximum
    t0 = results.parameters[results.param_names.index('t0')]
    pre_max = sum(data['time'] < t0)
    post_max = sum(data['time'] >= t0)
    new_row += [pre_max, post_max]

    # Add parameters and their errors to the new row
    new_row += list(results.parameters)
    new_row += [results.errors.get(p, 0) for p in results.param_names]

    # Calc chi-squared
    chisq, dof = utils.calc_model_chisq(data, results, fitted_model)
    new_row += [np.round(chisq, 2), dof]

    # Determine peak magnitude and decline rate
    b_max = fitted_model.source_peakabsmag('bessellb', 'ab')
    peak_phase = fitted_model.source.peakphase('bessellb')
    b_0 = fitted_model.source.bandmag('bessellb', 'ab', peak_phase)
    b_15 = fitted_model.source.bandmag('bessellb', 'ab', peak_phase + 15)
    delta_15 = b_15 - b_0
    new_row += [np.round(b_max, 2), np.round(delta_15, 3)]

    # Add fitting exit status message. Not all fitting routines include
    # this attribute, so we assign a default value of 'NONE'.
    message = getattr(results, 'message', 'NONE')
    new_row.append(message)
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


def run_classification_fits(
        blue_data, red_data, vparams, fit_func,
        priors_s2=None, priors_bg=None,
        kwargs_s2=None, kwargs_bg=None,
        show_plots=False):
    """Run light curve fits on a given target using the normal and 91bg model

    Fits are run for all data, ``red_data``, and then ``blue_data`` using salt2 and
    then the sn91bg model (i.e. three Salt2 fits followed by three 91bg fits).
    If a prior for t0 is not specified for a given value, then the salt2 t0
    value determined from the +/- 3 days is assumed.

    Args:
        blue_data (Table): Photometric data in just the 'red' bands
        red_data  (Table): Photometric data in just the 'blue' bands
        vparams    (iter): Iterable of param names to fit
        fit_func   (func): Function to use to run fits (eg. ``sncosmo.fit_lc``)
        priors_s2  (dict): Priors to use when fitting salt2
        priors_bg  (dict): Priors to use when fitting sn91bg
        kwargs_s2  (dict): Kwargs to pass ``fit_func`` when fitting salt2
        kwargs_bg  (dict): Kwargs to pass ``fit_func`` when fitting sn91bg
        show_plots (bool): Plot and display each individual fit

    Returns:
       Fit results and the fitted model for each model / data combination
    """

    fixed_params = {'z', 't0', 'x0', 'x1', 'c'} - set(vparams)
    _raise_unspecified_params(fixed_params, priors_s2, priors_bg)

    # Set default kwargs and protect against mutation
    priors_s2 = deepcopy(priors_s2) or dict()
    priors_bg = deepcopy(priors_bg) or dict()
    kwargs_s2 = deepcopy(kwargs_s2) if kwargs_s2 else dict()
    kwargs_bg = deepcopy(kwargs_bg) if kwargs_bg else dict()

    # Load models and set initial guess from priors
    salt2 = sncosmo.Model('salt2')
    sn91bg = sncosmo.Model('sn91bg')
    salt2.update(priors_s2)
    sn91bg.update(priors_bg)

    # Build iterator over data and models for each fit we want to run
    # We rely on the the salt2.4 model and full data table
    # being first in the iterator
    all_data = vstack([blue_data, red_data])
    data_tables = 2 * [all_data, blue_data, red_data]
    models = 3 * [salt2] + 3 * [sn91bg]
    kwargs = 3 * [kwargs_s2] + 3 * [kwargs_bg]
    fitting_data = zip(data_tables, models, kwargs)

    out_data = []
    salt2_t0 = None
    for data_table, model, kwarg in fitting_data:
        result, fitted_model = fit_func(data_table, model, vparams, **kwarg)
        out_data.append((result, fitted_model))

        # Update models to use salt2 t0 if not already specified in the priors
        if salt2_t0 is None:
            salt2_t0 = fitted_model.parameters[
                fitted_model.param_names.index('t0')]
            salt2.set(t0=priors_s2.get('t0', salt2_t0))
            sn91bg.set(t0=priors_bg.get('t0', salt2_t0))

            # We could also set default bounds on t0
            # kwargs_s2.setdefault('bounds', {})
            # kwargs_s2['bounds'].setdefault('t0', (s2_t0 - 3, s2_t0 + 3))
            # kwargs_bg.setdefault('bounds', {})
            # kwargs_bg['bounds'].setdefault('t0', (bg_t0 - 3, bg_t0 + 3))

        if show_plots:
            sncosmo.plot_lc(data_table, fitted_model)
            pyplot.show()

    return out_data


def tabulate_fit_results(
        data_iter, band_names, lambda_eff, fit_func, vparams, timeout_sec=90,
        config=None, out_path=None):
    """Tabulate fit results for a collection of data tables

    Args:
        data_iter   (iter): Iterable of photometric data for different SN
        band_names  (list): Name of bands included in ``data_iter``
        lambda_eff  (list): Effective wavelength for bands in ``band_names``
        vparams     (list): Name of parameters to vary
        fit_func    (func): Function to use to run fits
        timeout_sec  (int): Number of seconds before timeout fitting a SN
        config      (dict): Specifies priors / kwargs for fitting each model
        out_path     (str): Optionally cache progressive results to file

    Returns:
       An astropy table with fit results
    """

    # Set default kwargs
    config = deepcopy(config) or dict()

    # Add arguments to output table meta data
    out_table = create_empty_table()
    out_table.meta['band_names'] = band_names
    out_table.meta['lambda_eff'] = lambda_eff
    out_table.meta['fit_func'] = fit_func.__name__
    out_table.meta['vparams'] = vparams
    out_table.meta['timeout_sec'] = timeout_sec
    out_table.meta['out_path'] = str(out_path)

    for data in data_iter:
        blue_data, red_data = utils.split_data(
            data, band_names, lambda_eff, data.meta['redshift'])

        # Get fitting priors and kwargs
        obj_id = data.meta['obj_id']
        salt2_prior, salt2_kwargs, sn91bg_prior, sn91bg_kwargs = \
            utils.parse_config_dict(obj_id, config)

        try:
            with utils.timeout(timeout_sec):
                fit_results = run_classification_fits(
                    blue_data=blue_data,
                    red_data=red_data,
                    vparams=vparams,
                    fit_func=fit_func,
                    priors_s2=salt2_prior,
                    priors_bg=sn91bg_prior,
                    kwargs_s2=salt2_kwargs,
                    kwargs_bg=sn91bg_kwargs)

            # Create an iterator over data necessary to make each row
            band_sets = 2 * ['all', 'blue', 'red']
            data_tables = 2 * [data, red_data, blue_data]
            row_data_iter = zip(band_sets, data_tables, fit_results)

            # Populate the output table
            for bs, dt, (results, fitted_model) in row_data_iter:
                row = fit_results_to_table_row(dt, bs, results, fitted_model)
                out_table.add_row(row)

        except KeyboardInterrupt:
            raise

        except Exception as e:
            # Add a masked row so we have a record in the output table
            # indicating something went wrong.
            num_cols = len(out_table.colnames)
            row = np.full(num_cols, -99).tolist()
            mask = np.full(num_cols, True)

            row[0] = obj_id
            row[-1] = str(e).replace('\n', '')
            mask[0] = mask[-1] = False
            out_table.add_row(row, mask=mask)

        if out_path:
            out_table.write(out_path)

    return out_table


def classify_targets(fits_table, out_path=None):
    """Tabulate fitting coordinates for SNe based on their fit results

    See the ``create_empty_table`` function for information on the assumed
    input table format.

    Args:
        fits_table (Table): A table of fit results
        out_path     (str): Optionally write results to file

    Returns:
        An astropy table of fitting coordinates
    """

    # Convert input table to a DataFrame so we can leverage multi-indexing
    fits_df = fits_table.to_pandas()
    fits_df.set_index(['obj_id', 'source', 'band_set'], inplace=True)
    fits_df.dropna(inplace=True)

    out_table = Table(names=['obj_id', 'x', 'y'], dtype=['U100', float, float])
    for obj_id in fits_df.index.unique(level='obj_id'):

        try:

            salt2_blue_chisq = (
                    fits_df.loc[obj_id, 'salt2', 'blue']['chisq'] /
                    fits_df.loc[obj_id, 'salt2', 'blue']['ndof'])

            salt2_red_chisq = (
                    fits_df.loc[obj_id, 'salt2', 'red']['chisq'] /
                    fits_df.loc[obj_id, 'salt2', 'red']['ndof'])

            bg_blue_chisq = (
                    fits_df.loc[obj_id, 'sn91bg', 'blue']['chisq'] /
                    fits_df.loc[obj_id, 'sn91bg', 'blue']['ndof'])

            bg_red_chisq = (
                    fits_df.loc[obj_id, 'sn91bg', 'red']['chisq'] /
                    fits_df.loc[obj_id, 'sn91bg', 'red']['ndof'])

        except KeyError:
            continue

        else:
            x = salt2_blue_chisq - bg_blue_chisq
            y = salt2_red_chisq - bg_red_chisq
            out_table.add_row([obj_id, x, y])

            if out_path:
                out_table.write(out_path)

    return out_table
