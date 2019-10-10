#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``classification`` module tabulates light-curve fits for a given survey
and determines the corresponding classification coordinates:

.. math::

    x \def \chi^2_{blue}(Ia) - \chi^2_{blue}(91bg)
    y \def \chi^2_{red}(Ia) - \chi^2_{red}(91bg)
"""

from copy import deepcopy
from functools import partial
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

from . import fitting
from . import utils


def _get_fitting_func(fitting_method, band_names=None, lambda_eff=None):
    """Get the correct fitting function to match the described fittingmethod

    Args:
        fitting_method (str): Desired fitting method ('band' or 'collective')

    Returns:
        A callable from the ``fitting`` module
    """

    if fitting_method.lower() == 'band':
        return fitting.run_band_fits

    elif fitting_method.lower() == 'collective':
        return partial(fitting.run_collective_fits,
                       band_names=band_names,
                       lambda_eff=lambda_eff)

    else:
        raise ValueError(
            f"Fitting method should be 'band' "
            f"or 'collective' not {fitting_method}")


def tabulate_fit_results(
        data_iter, band_names, lambda_eff, fit_func, fitting_method='band',
        config=None, out_path=None):
    """Tabulate fit results for a collection of data tables

    Results already written to out_path are skipped.

    Args:
        data_iter  (iter): Iterable of photometric data for different SN
        band_names (list): Name of bands included in ``data_iter``
        lambda_eff (list): Effective wavelength for bands in ``band_names``
        fit_func   (func): Function to use to run fits
        config     (dict): Specifies priors / kwargs for fitting each model
        out_path    (str): Optionally cache results to file

    Returns:
       An astropy table with fit results
    """

    # Set default kwargs
    config = deepcopy(config) or dict()
    fitting_func = _get_fitting_func(fitting_method, band_names, lambda_eff)

    # Add meta_data to output table meta data
    if Path(out_path).exists():
        out_table = Table.read(out_path)

    else:
        out_table = Table(names=['obj_id', 'message'], dtype=['U20', 'U10000'])

    out_table.meta['band_names'] = band_names
    out_table.meta['lambda_eff'] = lambda_eff
    out_table.meta['fit_func'] = fit_func.__name__
    out_table.meta['out_path'] = str(out_path)

    for data in data_iter:
        # Get fitting priors and kwargs
        obj_id = data.meta['obj_id']
        if obj_id in out_table['obj_id']:
            continue

        salt2_prior, salt2_kwargs, sn91bg_prior, sn91bg_kwargs = \
            utils.parse_config_dict(obj_id, config)

        try:
            fit_results = fitting_func(
                obj_id=obj_id,
                data=data,
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


# noinspection PyPep8
def classify_targets(
        fits_table, band_names=None, lambda_eff=None, out_path=None):
    """Tabulate fitting coordinates for SNe based on their fit results

    See the ``create_empty_table`` function for the assumed input table format.
    If ``band_names`` or ``lambda_eff`` are not given, they are taken from
    ``fits_table.meta``. Any targets having one or more fits with the string
    'failed' in the message are skipped (case insensitive).

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

    # Keep only objects that don't have failed fits in any band
    fits_df = fits_table.to_pandas()
    failed_fits = fits_df['message'].str.lower().str.contains('failed')
    failed_ids = fits_df[failed_fits].obj_id.unique()
    good_fits = fits_df[~fits_df.obj_id.isin(failed_ids)]
    good_fits.set_index(['obj_id', 'source'], inplace=True)

    out_table = Table(names=['obj_id', 'x', 'y'], dtype=['U100', float, float])
    for obj_id in good_fits.index.unique(level='obj_id'):

        try:
            hsiao_data = good_fits.loc[obj_id, 'hsiao_x1']
            redshift = hsiao_data[hsiao_data['band'] == 'all']['z'][0]
            blue_bands, red_bands = utils.split_bands(
                band_names, lambda_eff, redshift)

            blue_bands = np.concatenate([blue_bands, ['blue']])
            red_bands = np.concatenate([red_bands, ['red']])

            hsiao_blue = hsiao_data[hsiao_data['band'].isin(blue_bands)]
            hsiao_red = hsiao_data[hsiao_data['band'].isin(red_bands)]
            hsiao_blue_chisq = hsiao_blue['chisq'].sum() / hsiao_blue[
                'ndof'].sum()
            hsiao_red_chisq = hsiao_red['chisq'].sum() / hsiao_red[
                'ndof'].sum()

            sn91bg_data = good_fits.loc[obj_id, 'sn91bg']
            sn91bg_blue = sn91bg_data[sn91bg_data['band'].isin(blue_bands)]
            sn91bg_red = sn91bg_data[sn91bg_data['band'].isin(red_bands)]
            sn91bg_blue_chisq = sn91bg_blue['chisq'].sum() / sn91bg_blue[
                'ndof'].sum()
            sn91bg_red_chisq = sn91bg_red['chisq'].sum() / sn91bg_red[
                'ndof'].sum()

        except KeyError:
            continue

        else:
            x = hsiao_blue_chisq - sn91bg_blue_chisq
            y = hsiao_red_chisq - sn91bg_red_chisq
            out_table.add_row([obj_id, x, y])

            if out_path:
                out_table.write(out_path)

    return out_table
