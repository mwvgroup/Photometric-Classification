#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module provides utility functions useful when fitting data with \
SNCosmo.
"""

import numpy as np
from astropy.table import Table


def create_empty_summary_table(band_names,  extra_cols=[]):
    """Returns a table with columns:

         cid, *num_points_<band_names>, z, t0, x0, x1, z_err, t0_err, x0_err,
         x1_err, c_err, chi, dof, message.
    """

    names = extra_cols + ['cid']
    names.extend(['num_points_' + band for band in band_names])

    param_names = ['z', 't0', 'x0', 'x1', 'c']
    names.extend(param_names)
    names.extend((p + '_err' for p in param_names))
    names.extend(('chi', 'dof', 'message'))
    out_table = Table(names=names, dtype=[object for _ in names])

    return out_table


def count_points_per_band(band_list, all_band_names):
    """Determine number of data points per band

    count the number of times each element in <all_band_names> appears in
    <band_list>.

    Returns:
        A list with the number of counts for each element in <all_band_names>
    """

    band_names, band_counts = np.unique(band_list, return_counts=True)
    count_dict = dict(zip(band_names, band_counts))
    return [count_dict.get(band, 0) for band in all_band_names]
