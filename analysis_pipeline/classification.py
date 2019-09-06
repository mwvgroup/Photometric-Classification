#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``classification`` module fits light curves, tabulates results, and
determines the corresponding classification coordinates.
"""

from copy import deepcopy

import numpy as np
import sncosmo
from astropy.table import Table

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
        - dof
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
    dtype += [float, float, float, float, 'U1000']

    # Unless otherwise specified, we default to returning a masked table
    kwargs = kwargs.copy()
    kwargs.setdefault('masked', True)
    return Table(names=names, dtype=dtype, **kwargs)


def fit_results_to_table_row(data, band_set, results, fitted_model):
    """Format fit results so they can be appended to an astropy table

    See the ``create_empty_table`` function for information on the assumed
    table format. Fitted model is assumed to have the parameters 'z', 't0',
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
    t0 = fitted_model.parameters[fitted_model.param_names.index('t0')]
    pre_max = sum(data['time'] < t0)
    post_max = sum(data['time'] >= t0)
    new_row += [pre_max, post_max]

    # Add parameters and their errors to the new row
    new_row += list(fitted_model.parameters)
    new_row += [results.errors.get(p, 0) for p in fitted_model.param_names]

    # Calc chi-squared
    chisq, dof = utils.calc_model_chisq(data, fitted_model)
    new_row += [chisq, dof]

    # Determine peak magnitude and decline rate
    b_max = fitted_model.source_peakabsmag('bessellb', 'ab')
    peak_phase = fitted_model.source.peakphase('bessellb')
    b_0 = fitted_model.source.bandmag('bessellb', 'ab', peak_phase)
    b_15 = fitted_model.source.bandmag('bessellb', 'ab', peak_phase + 15)
    delta_15 = b_15 - b_0
    new_row += [b_max, delta_15]

    # Add fitting exit status message. Not all fitting routines include
    # this attribute, se we assign a default value.
    message = getattr(results, 'message', 'NONE')
    new_row.append(message)
    return new_row


def run_fits(all_data, blue_data, red_data, vparams, fit_func,
             kwargs_s2=None, kwargs_bg=None, show_plots=False):
    """Run light curve fits on a given target using the normal and 91bg model

    Fits are run for all_data, red_data, and then blue_data using salt2 and
    then the sn91bg model.

    Args:
        all_data  (Table): Photometric data in all bands
        blue_data (Table): Photometric data in just the 'red' bands
        red_data  (Table): Photometric data in just the 'blue' bands
        vparams    (iter): Iterable of param names to fit
        fit_func   (func): Function to use to run fits (eg. ``sncosmo.fit_lc``)
        kwargs_s2  (dict): Kwargs to pass ``fit_func`` when fitting salt2
        kwargs_bg  (dict): Kwargs to pass ``fit_func`` when fitting sn91bg
        show_plots (bool): Plot and display each individual fit

    Returns:
       Fit results and the fitted model for each model / data combination
    """

    # Copy kwargs if given to avoid mutation
    kwargs_s2 = deepcopy(kwargs_s2) if kwargs_s2 else dict()
    kwargs_bg = deepcopy(kwargs_bg) if kwargs_bg else dict()

    # Load models
    salt2 = sncosmo.Model('salt2')
    sn91bg = sncosmo.Model('sn91bg')

    # Set redshift if not fitting for it
    if 'z' not in vparams:
        z = all_data.meta['redshift']
        salt2.set(z=z)
        sn91bg.set(z=z)

    # Build iterator over data and models for each fit we want to run
    data_tables = 2 * [all_data, red_data, blue_data]
    models = 3 * [salt2] + 3 * [sn91bg]
    kwargs = 3 * [kwargs_s2] + 3 * [kwargs_bg]
    fitting_data = zip(data_tables, models, kwargs)

    out_data = []
    salt2_t0 = None
    for data_table, model, kwarg in fitting_data:

        # Default to bounding t0 to the salt2.4 t0 +/- 3 days
        if salt2_t0 is not None and 't0' not in kwargs_bg.get('bounds', {}):
            kwargs_bg.setdefault('bounds', {})
            kwargs_bg['bounds']['t0'] = (salt2_t0 - 3, salt2_t0 + 3)

        result, fitted_model = fit_func(data_table, model, vparams, **kwarg)
        out_data.append((result, fitted_model))

        # We rely on the fact that the salt2.4 model and full data table
        # are first in the iterator
        if salt2_t0 is None:
            salt2_t0 = fitted_model.parameters[fitted_model.param_names.index('t0')]
            salt2.set(t0=salt2_t0)
            sn91bg.set(t0=salt2_t0)

        if show_plots:
            sncosmo.plot_lc(data_table, fitted_model)

    return out_data


def tabulate_fit_results(
        data_iter, band_names, lambda_eff, fit_func, vparams,
        timeout_sec=90, kwargs_s2=None, kwargs_bg=None, out_path=None):
    """Tabulate fit results for a collection of data tables

    Args:
        data_iter  (iter): Iterable of photometric data for different SN
        band_names (list): Name of bands included in ``data_iter``
        lambda_eff (list): Effective wavelength for each band in ``band_names``
        vparams    (list): Name of parameters to vary
        fit_func   (func): Function to use to run fits
        timeout_sec (int): Number of seconds before timeout fitting for target
        kwargs_s2  (dict): Kwargs to pass ``fit_func`` when fitting salt2
        kwargs_bg  (dict): Kwargs to pass ``fit_func`` when fitting sn91bg
        out_path    (str): Optionally cache progressive results to file

    Returns:
       An astropy table with fit results
    """

    # Add arguments to output table meta data
    out_table = create_empty_table()
    out_table.meta['band_names'] = band_names
    out_table.meta['lambda_eff'] = lambda_eff
    out_table.meta['fit_func'] = fit_func.__name__
    out_table.meta['vparams'] = vparams
    out_table.meta['timeout_sec'] = timeout_sec
    out_table.meta['kwargs_s2'] = kwargs_s2
    out_table.meta['kwargs_bg'] = kwargs_bg
    out_table.meta['out_path'] = str(out_path)

    for data in data_iter:
        blue_data, red_data = utils.split_data(data, band_names, lambda_eff)

        try:
            with utils.timeout(timeout_sec):
                fit_results = run_fits(
                    data, blue_data, red_data,
                    vparams,
                    fit_func,
                    kwargs_s2,
                    kwargs_bg)

        except KeyboardInterrupt:
            raise

        except Exception as e:
            # Add a masked row so we have a record in the output table
            # indicating something went wrong.
            num_cols = len(out_table.colnames)
            row = np.full(num_cols, -99).tolist()
            mask = np.full(num_cols, True)

            row[0] = data.meta['obj_id']
            row[-1] = str(e).replace('\n', '')
            mask[0] = mask[-1] = False
            out_table.add_row(row, mask=mask)

        else:
            # Create an iterator over data necessary to make each row
            band_sets = 2 * ['all', 'blue', 'red']
            data_tables = 2 * [data, red_data, blue_data]
            row_data_iter = zip(band_sets, data_tables, fit_results)

            # Populate the output table
            for bs, dt, (results, fitted_model) in row_data_iter:
                row = fit_results_to_table_row(dt, bs, results, fitted_model)
                out_table.add_row(row)

        if out_path:
            out_table.write(out_path)

    return out_table


def classify_targets(fits_table, out_path=None):
    """Tabulate classification coordinates for SNe based on their fit results

    See the ``create_empty_table`` function for information on the assumed
    input table format.

    Args:
        fits_table (Table): A table of fit results
        out_path    (str): Optionally write results to file

    Returns:
        An astropy table of classification coordinates
    """

    # Convert input table to a DataFrame so we can leverage multi-indexing
    fits_df = fits_table.to_pandas()
    fits_df.set_index(['obj_id', 'source', 'band_set'], inplace=True)
    fits_df.dropna(inplace=True)

    out_table = Table(names=['obj_id', 'x', 'y'], dtype=['U100', float, float])
    for obj_id in fits_df.index.unique(level='obj_id'):

        # We expect fit for 3 band sets x 3 models = 6 fits total
        # Less than 6 means one of the fits failed and was dropped by the
        # `drop_na` call above
        if len(fits_df.loc[obj_id]) < 6:
            continue

        salt2_blue_chisq = (fits_df.loc[obj_id, 'salt2', 'blue']['chisq'] /
                            fits_df.loc[obj_id, 'salt2', 'blue']['ndof'])

        salt2_red_chisq = (fits_df.loc[obj_id, 'salt2', 'red']['chisq'] /
                           fits_df.loc[obj_id, 'salt2', 'red']['ndof'])

        bg_blue_chisq = (fits_df.loc[obj_id, 'sn91bg', 'blue']['chisq'] /
                         fits_df.loc[obj_id, 'sn91bg', 'blue']['ndof'])

        bg_red_chisq = (fits_df.loc[obj_id, 'sn91bg', 'red']['chisq'] /
                        fits_df.loc[obj_id, 'sn91bg', 'red']['ndof'])

        x = salt2_blue_chisq - bg_blue_chisq
        y = salt2_red_chisq - bg_red_chisq
        out_table.add_row([obj_id, x, y])

    if out_path:
        out_table.write(out_path)

    return out_table
