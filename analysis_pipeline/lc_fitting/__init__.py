#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module handles the nested sampling and fitting of light curve data."""

from ._fit_funcs import (create_empty_summary_table,
                         fit_lc,
                         get_sampled_model,
                         nest_lc)

from ._fit_n_params import split_data
from ._module_interface import LCFitting
