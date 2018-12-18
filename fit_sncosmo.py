#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module fits SDSS light curves using sncosmo"""

import numpy as np
import sncosmo
from astropy.table import Table

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
        sncosmo_table['flux'] = all_sn_data['MAG']
        sncosmo_table['fluxerr'] = all_sn_data['MERR']
        sncosmo_table['zp'] = np.full(len(all_sn_data), 25)
        sncosmo_table['zpsys'] = np.full(len(all_sn_data), 'ab')

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

        yield result
