#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module provides access to data from the third data release of the
Carnegie Supernova Project (CSP). Light curve data is provided from
Krisciunas et al. 2017. The master_table object provides data from Table 1
of the same paper.

No data cuts are applied by this module or its functions.

For more information on CSP data products, see:
    https://csp.obs.carnegiescience.edu
"""

from ._data_access_funcs import get_data_for_id
from ._data_access_funcs import get_input_for_id
from ._data_access_funcs import iter_sncosmo_input
from ._data_access_funcs import master_table
from ._module_meta_data import band_names
from ._module_meta_data import lambda_effective
