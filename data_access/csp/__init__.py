#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module provides access to data from the third data release of the
Carnegie Supernova Project (CSP). See the package level docstring for usage
instructions.

For more information on CSP data products, see:
    https://csp.obs.carnegiescience.edu
"""

from ._data_access_funcs import get_data_for_id
from ._data_access_funcs import get_input_for_id
from ._data_access_funcs import iter_sncosmo_input
from ._data_access_funcs import master_table
from ._module_meta_data import band_names
from ._module_meta_data import lambda_effective
