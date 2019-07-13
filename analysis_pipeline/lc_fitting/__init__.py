#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module handles the nested sampling and fitting of light curve data."""

from ._fit_funcs import calc_chisq, create_results_table, fit_lc
from ._iter_fitting import run_iter_fitting
from ._sampling import get_sampled_model, nest_lc
