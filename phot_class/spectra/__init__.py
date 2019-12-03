#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module tabulates the properties of spectral features such as area,
velocity, and  equivalent width. All functions in this module are built to
support ``uarray`` objects from the ``uncertainties`` package as inputs.
"""

from astropy.table import Table as _Table

from .calc_properties import (bin_spectrum, correct_extinction, dust_map,
                              feature_area, feature_pew, feature_velocity,
                              find_peak_wavelength, guess_feature_bounds,
                              line_locations)
from .graphical_interface import SpectrumInspector


def tabulate_spectral_properties(
        data_iter, nstep=5, bin_size=3, method='avg', rv=3.1):
    """Tabulate spectral properties for multiple spectra of the same object

    Spectra are rest-framed and corrected for MW extinction using the
    Schlegel et al. 98 dust map and the Fitzpatrick et al. 99 extinction law.

    Args:
        data_iter (iter[Table]): Iterable of spectroscopic data tables
        nstep             (int): The number of sampling steps to take
        bin_size        (float): The width of the bins (Default: 5)
        method            (str): Either 'avg' or 'sum' the values of each bin
        rv              (float): Rv value to use for extinction

    Returns:
        A Table with measurements for each spectrum and feature
    """

    table_rows = []
    for spectrum in data_iter:
        inspector = SpectrumInspector(spectrum)
        table_rows = inspector.run(
            nstep=nstep, bin_size=bin_size, method=method, rv=rv)

    if not table_rows:
        table_rows = None

    # Format results as a table
    col_names = ['obj_id', 'sid', 'date', 'type', 'feat_name', 'feat_start', 'feat_end']
    dtype = ['U100', 'U100', 'U100', 'U100', 'U20', float, float]
    for value in ('vel', 'pew', 'area'):
        col_names.append(value)
        col_names.append(value + '_err')
        col_names.append(value + '_samperr')
        dtype += [float, float, float]

    col_names.append('msg')
    dtype.append('U1000')

    return _Table(rows=table_rows, names=col_names, dtype=dtype)
