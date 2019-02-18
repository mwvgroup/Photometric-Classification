#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module provides access to data from the SDSS-II SN Catalog Data
Release.
"""

from ._data_access_funcs import get_data_for_id
from ._data_access_funcs import get_input_for_id
from ._data_access_funcs import iter_sncosmo_input
from ._data_access_funcs import master_table
from ._module_meta_data import band_names
from ._module_meta_data import lambda_effective