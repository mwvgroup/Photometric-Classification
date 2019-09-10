#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``fit_funcs`` module provides wrappers that combine the the various
fitting functions in ``sncosmo``. Importantly, the fitting functions in this
module guarantee that arguments will not be mutated, which is not true for
``sncosmo`` in general.
"""

from copy import deepcopy

import sncosmo


def _copy_data(*args):
    """Return a copy of the passed arguments"""

    return tuple(deepcopy(a) for a in args)


def simple_fit(data, model, vparam_names, **kwargs):
    """Fit light curves using the basic ``fit_lc`` functionality in ``sncosmo``

    Args:
        data            (Table): Table of photometric data
        model           (Model): The model to fit
        vparam_names (iterable): Model parameters to vary
        Any other parameters for ``sncosmo.fit_lc``

    Returns:
         - A dictionary like object with fit results
         - A fitted ``Model`` instance
    """

    data, model, vparam_names, kwargs = \
        _copy_data(data, model, vparam_names, kwargs)

    return sncosmo.fit_lc(data, model, vparam_names, **kwargs)


def nest_fit(data, model, vparam_names, **kwargs):
    """Fit light curves using nested sampling

    Args:
        data            (Table): Table of photometric data
        model           (Model): The model to fit
        vparam_names (iterable): Model parameters to vary
        Any other parameters for ``sncosmo.nest_lc``

    Returns:
         - A dictionary like object with fit results
         - A fitted ``Model`` instance
    """

    data, model, vparam_names, kwargs = \
        _copy_data(data, model, vparam_names, kwargs)

    return sncosmo.nest_lc(data, model, vparam_names, **kwargs)


def mcmc_fit(data, model, vparam_names, **kwargs):
    """Fit light curves Monte Carlo sampling

    Args:
        data            (Table): Table of photometric data
        model           (Model): The model to fit
        vparam_names (iterable): Model parameters to vary
        Any other parameters for ``sncosmo.mcmc_lc``

    Returns:
         - A dictionary like object with fit results
         - A fitted ``Model`` instance
    """

    data, model, vparam_names, kwargs = \
        _copy_data(data, model, vparam_names, kwargs)

    return sncosmo.mcmc_lc(data, model, vparam_names, **kwargs)


def nested_simple_fit(data, model, vparam_names, **kwargs):
    """Use nested sampling for initial parameters and then perform a simple fit

    Args:
        data            (Table): Table of photometric data
        model           (Model): The model to fit
        vparam_names (iterable): Model parameters to vary
        Any other parameters for ``sncosmo.nest_lc`` and ``sncosmo.fit_lc``

    Returns:
         - A dictionary like object with fit results
         - A fitted ``Model`` instance
    """

    data, model, vparam_names, kwargs = \
        _copy_data(data, model, vparam_names, kwargs)

    _, nested_model = sncosmo.nest_lc(data, model, vparam_names, **kwargs)
    return sncosmo.fit_lc(data, model, vparam_names, **kwargs)


def nested_mcmc_fit(data, model, vparam_names, **kwargs):
    """Use nested sampling for initial parameters and then perform a MCMC fit

    Args:
        data            (Table): Table of photometric data
        model           (Model): The model to fit
        vparam_names (iterable): Model parameters to vary
        Any other parameters for ``sncosmo.nest_lc`` and ``sncosmo.mcmc_lc``

    Returns:
         - A dictionary like object with fit results
         - A fitted ``Model`` instance
    """

    data, model, vparam_names, kwargs = \
        _copy_data(data, model, vparam_names, kwargs)

    _, nested_model = sncosmo.nest_lc(data, model, vparam_names, **kwargs)
    return sncosmo.mcmc_lc(data, model, vparam_names, **kwargs)