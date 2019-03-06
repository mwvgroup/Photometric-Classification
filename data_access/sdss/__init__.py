#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module provides access to data from the SDSS-II SN Catalog Data
Release. See the package level docstring for usage instructions.

For more information on SDSS data products, see:
    https://data.sdss.org/sas/dr10/boss/papers/supernova/
"""

from ._data_access_funcs import get_data_for_id
from ._data_access_funcs import get_input_for_id
from ._data_access_funcs import iter_sncosmo_input
from ._data_access_funcs import master_table
from ._module_meta_data import band_names
from ._module_meta_data import lambda_effective
