#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module fits SDSS light curves using sncosmo"""

import os
import logging
from warnings import warn

import numpy as np
import sncosmo
from astropy.table import Table
from matplotlib import pyplot as plt
from tqdm import tqdm

from parse_sn_data import get_cid_data, MASTER_TABLE


@np.vectorize
def map_index_to_band(i):
    return 'sdss' + 'ugriz'[i]


def iter_sncosmo_input():
    """Iterate through SDSS supernova and yield the SNCosmo input tables

    Yields:
        An astropy table formatted for use with SNCosmo
    """

    for cid in tqdm(MASTER_TABLE['CID']):
        sncosmo_table = Table()
        all_sn_data = get_cid_data(cid)

        sncosmo_table['time'] = all_sn_data['MJD']
        sncosmo_table['band'] = map_index_to_band(all_sn_data['FILT'])
        sncosmo_table['flux'] = all_sn_data['FLUX']
        sncosmo_table['fluxerr'] = all_sn_data['FLUXERR']
        sncosmo_table['zp'] = np.full(len(all_sn_data), 25)
        sncosmo_table['zpsys'] = np.full(len(all_sn_data), 'ab')
        sncosmo_table.meta = all_sn_data.meta
        sncosmo_table.meta['cid'] = cid

        yield sncosmo_table


def fit_cid(input_table, model_name, out_dir):

    cid = input_table.meta['cid']

    # Configure model for fitting
    model = sncosmo.Model(source=model_name)
    z = input_table.meta['redshift']
    if z == -9.:
        params_to_fit = ['z', 't0', 'x0', 'x1', 'c']
        bounds = {'z': (0.002, 1)}

    else:
        model.set(z=z)
        params_to_fit = ['t0', 'x0', 'x1', 'c']
        bounds = None

    # Run fit
    result, fitted_model = sncosmo.fit_lc(
        input_table, model, params_to_fit, bounds=bounds)

    # save plot
    try:
        fig_path = os.path.join(out_dir, '{}.pdf'.format(cid))
        sncosmo.plot_lc(input_table, model=fitted_model, errors=result.errors)
        plt.savefig(fig_path)
        plt.close()

    except ValueError as e:
        warn('Could not plot figure: {}'.format(e))

    return result


def fit_sdss_data(out_dir, model_name='salt2'):
    """Fit smp data with sncosmo

    Files are named as <out_dir>/<target cid>.txt

    Args:
        out_dir    (str): Where to write fit results
        model_name (str): Model to use for fitting. Default = salt2
    """

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Store warnings and errors
    logging.captureWarnings(True)
    log_path = os.path.join(out_dir, 'summary.log')
    logging.basicConfig(filename=log_path, level=logging.WARNING)

    # Table to store summary of fit results
    names = ['cid', 'class', 'z', 'fit_z', 't0', 'x0', 'x1', 'c', 'chi', 'dof', 'message']
    out_table = Table(names=names, dtype=[object for _ in names])
    out_path = os.path.join(out_dir, 'summary.csv')
    out_table.write(out_path, overwrite=True)

    for input_table in iter_sncosmo_input():
        try:
            result = fit_cid(input_table, model_name, out_dir)

        except (sncosmo.fitting.DataQualityError, RuntimeError) as e:
            logging.error('{}: {}\n'.format(input_table.meta['cid'], e))
            continue

        # Update summary file
        z_was_fit = int(input_table.meta['redshift'] == -9)
        new_row = [input_table.meta['cid'], input_table.meta['classification'], z_was_fit]
        new_row.extend(result.parameters)
        new_row.append(result.chisq)
        new_row.append(result.ndof)
        new_row.append(result.message)
        out_table.add_row(new_row)
        out_table.write(out_path, overwrite=True)


if __name__ == '__main__':
    fit_sdss_data('./sncosmo_results')
