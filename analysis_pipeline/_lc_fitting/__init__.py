#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module handles the nested sampling and fitting of light curve data."""

from ._fit_funcs import create_empty_summary_table
from ._fit_funcs import fit_lc
from ._fit_funcs import get_sampled_model
from ._fit_funcs import nest_lc
from ._iter_fitting import iter_all_fits, split_data
