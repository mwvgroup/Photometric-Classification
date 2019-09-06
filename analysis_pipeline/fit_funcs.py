#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``fit_funcs`` module provides wrappers that combine the the various
fitting functions in ``sncosmo``.
"""

from copy import deepcopy

import sncosmo as _sncosmo


def simple_fit(data, model, vparam_names, **kwargs):
    """Fit light curves using the basic ``fit_lc`` functionality in ``sncosmo``

    Args:
        data (Table): Table of photometric data
        model (Model): The model to fit
        vparam_names (iterable) Model parameters to vary

    Returns:
         - A dictionary like object with fit results
         - A fitted ``Model`` instance
    """

    kwargs = deepcopy(kwargs)
    return _sncosmo.fit_lc(data, model, vparam_names, **kwargs)


def nest_fit(data, model, vparam_names, **kwargs):
    """Fit light curves using nested sampling

    Args:
        data (Table): Table of photometric data
        model (Model): The model to fit
        vparam_names (iterable) Model parameters to vary

    Returns:
         - A dictionary like object with fit results
         - A fitted ``Model`` instance
    """

    kwargs = deepcopy(kwargs)
    return _sncosmo.nest_lc(data, model, vparam_names, **kwargs)


def mcmc_fit(data, model, vparam_names, **kwargs):
    """Fit light curves Monte Carlo sampling

    Args:
        data (Table): Table of photometric data
        model (Model): The model to fit
        vparam_names (iterable) Model parameters to vary

    Returns:
         - A dictionary like object with fit results
         - A fitted ``Model`` instance
    """

    kwargs = deepcopy(kwargs)
    return _sncosmo.mcmc_lc(data, model, vparam_names, **kwargs)


def nested_simple_fit(data, model, vparam_names, **kwargs):
    """Use nested sampling for initial parameters and then perform a simple fit

    Args:
        data (Table): Table of photometric data
        model (Model): The model to fit
        vparam_names (iterable) Model parameters to vary

    Returns:
         - A dictionary like object with fit results
         - A fitted ``Model`` instance
    """

    kwargs = deepcopy(kwargs)
    _, nested_model = _sncosmo.nest_lc(data, model, vparam_names, **kwargs)
    return _sncosmo.fit_lc(data, model, vparam_names, **kwargs)


def nested_mcmc_fit(data, model, vparam_names, **kwargs):
    """Use nested sampling for initial parameters and then perform a MCMC fit

    Args:
        data (Table): Table of photometric data
        model (Model): The model to fit
        vparam_names (iterable) Model parameters to vary

    Returns:
         - A dictionary like object with fit results
         - A fitted ``Model`` instance
    """

    kwargs = deepcopy(kwargs)
    _, nested_model = _sncosmo.nest_lc(data, model, vparam_names, **kwargs)
    return _sncosmo.mcmc_lc(data, model, vparam_names, **kwargs)
