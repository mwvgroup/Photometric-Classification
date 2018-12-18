#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module fits SDSS light curves using sncosmo"""

import numpy as np
import sncosmo
from astropy.table import Table
from matplotlib import pyplot as plt

from parse_sn_data import get_cid_data, MASTER_TABLE


@np.vectorize
def map_index_to_band(i):
    return 'sdss' + 'ugriz'[i]


def iter_sncosmo_input():
    """Iterate through SDSS supernova and yield the SNCosmo input tables

    Yields:
        An astropy table formatted for use with SNCosmo
    """

    for cid in MASTER_TABLE['CID']:
        sncosmo_table = Table()
        all_sn_data = get_cid_data(cid)

        sncosmo_table['time'] = all_sn_data['MJD']
        sncosmo_table['band'] = map_index_to_band(all_sn_data['FILT'])
        sncosmo_table['flux'] = all_sn_data['FLUX']
        sncosmo_table['fluxerr'] = all_sn_data['FLUXERR']
        sncosmo_table['zp'] = np.full(len(all_sn_data), 25)
        sncosmo_table['zpsys'] = np.full(len(all_sn_data), 'ab')
        sncosmo_table.meta = all_sn_data.meta

        yield sncosmo_table


def iter_fit_sdss_data():

    for input_table in iter_sncosmo_input():
        model = sncosmo.Model(source='salt2')
        z = input_table.meta['redshift']
        if z != -9.:
            model.set(z=z)
            params_to_fit = ['t0', 'x0', 'x1', 'c']
            bounds = None

        else:
            params_to_fit = ['z', 't0', 'x0', 'x1', 'c']
            bounds = {'z': (0.02, 1)}

        result, fitted_model = sncosmo.fit_lc(
            input_table, model, params_to_fit, bounds=bounds)

        yield input_table, result, fitted_model


if __name__ == '__main__':

    for input_table, result, fitted_model in iter_fit_sdss_data():
        sncosmo.plot_lc(input_table, model=fitted_model, errors=result.errors)
        plt.show()