#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``classification`` module fits light curves, tabulates results, and
determines the corresponding classification coordinates.
"""

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
    new_row += fitted_model.parameters
    new_row += [results.errors[p] for p in fitted_model.param_names]

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


def run_fits(all_data, red_data, blue_data, vparams, fit_func,
             kwargs_s2=None, kwargs_bg=None):
    """Run light curve fits on a given target using the normal and 91bg model

    Fits are run for all_data, red_data, and then blue_data using salt2 and
    then the sn91bg model.

    Args:
        all_data  (Table): Photometric data in all bands
        red_data  (Table): Photometric data in just the 'blue' bands
        blue_data (Table): Photometric data in just the 'red' bands
        vparams    (iter): Iterable of param names to fit
        fit_func   (func): Function to use to run fits (eg. ``sncosmo.fit_lc``)
        kwargs_s2  (dict): Kwargs to pass ``fit_func`` when fitting salt2
        kwargs_bg  (dict): Kwargs to pass ``fit_func`` when fitting sn91bg

    Returns:
       Fit results and the fitted model for each model / data combination
    """

    # Copy kwargs if given to avoid mutation
    kwargs_s2 = kwargs_s2.copy() if kwargs_s2 else dict()
    kwargs_bg = kwargs_bg.copy() if kwargs_bg else dict()

    # Load models
    salt2 = sncosmo.Model('salt2')
    sn91bg = sncosmo.Model('sn91bg')

    # Set redshift if not fitting for it
    if 'z' not in vparams:
        z = all_data.meta['redshift']
        salt2.set(z=z)
        sn91bg.set(z=z)

    # Fit salt2 model to determine t0
    norm_result_all, norm_fit_all = fit_func(all_data, salt2, vparams,
                                             **kwargs_s2)
    t0 = norm_fit_all.parameters[norm_fit_all.param_names.index('t0')]

    # Set initial t0 value for remainder of fits. Unless a set of bounds is
    # already specified in the kwargs, we bound the 91bg t0 value since we
    # expect it to be close to the slat2 value
    salt2.set(t0=t0)
    sn91bg.set(t0=t0)
    if 't0' not in kwargs_bg.get('bounds', {}):
        kwargs_bg.setdefault('bounds', {})
        kwargs_bg['bounds']['t0'] = (t0 - 3, t0 + 3)

    # Run the remaining fits
    norm_all = (norm_result_all, norm_fit_all)
    norm_blue = fit_func(blue_data, salt2, vparams, **kwargs_s2)
    norm_red = fit_func(red_data, salt2, vparams, **kwargs_s2)
    bg_all = fit_func(all_data, sn91bg, vparams, **kwargs_bg)
    bg_blue = fit_func(blue_data, sn91bg, vparams, **kwargs_bg)
    bg_red = fit_func(red_data, sn91bg, vparams, **kwargs_bg)

    return norm_all, norm_blue, norm_red, bg_all, bg_blue, bg_red


def classify_data(data_iter, band_names, lambda_eff, fit_func, vparams,
                  timeout_seconds=90, kwargs_s2=None, kwargs_bg=None,
                  out_path=None):
    """Tabulate fit results for a collection of data tables

    Args:
        data_iter  (iter): Iterable of hhotometric data for different SN
        band_names (list): Name of bands included in ``data_iter``
        vparams    (list): Effective wavelength for each band in ``band_names``
        fit_func   (func): Function to use to run fits
        kwargs_s2  (dict): Kwargs to pass ``fit_func`` when fitting salt2
        kwargs_bg  (dict): Kwargs to pass ``fit_func`` when fitting sn91bg

    Returns:
       An astropy table with fit results
    """

    out_table = create_empty_table()
    for data in data_iter:
        blue_data, red_data = utils.split_data(data, band_names, lambda_eff)

        try:
            with utils.timeout(timeout_seconds):
                fit_results = run_fits(
                    data, red_data, blue_data,
                    vparams,
                    fit_func,
                    kwargs_s2,
                    kwargs_bg)

        except KeyboardInterrupt:
            raise

        except:
            pass  # Todo: add masked row with error message

        else:
            # Todo: Add fit results to table
            band_sets = ('all', 'blue', 'red', 'all', 'blue', 'red')
            data_tables = (data, red_data, blue_data, data, red_data, blue_data)
            for (results, fitted_model) in fit_results:
                row = fit_results_to_table_row(data, 'all', results, fitted_model)
                out_table.add_row(row)

        if out_path:
            out_table.write(out_path)

    return out_table
