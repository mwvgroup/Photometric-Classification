#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``models`` module defines a ``Source`` class for modeling
91bg-like supernovae with ``sncosmo``.

Available Model Versions
------------------------

There are multiple versions of the 91bg model included in this package. For
information on how each version determines the flux for a given set of
parameters, see the documentation for the given source class.

+--------+------------------+-------------------------------------------------+
| Name   | Version          | Description                                     |
+========+==================+=================================================+
| sn91bg | 'salt2_phase'    | 1991bg model restricted to a phase similar to   |
|        | (Default)        | Salt 2.4 (extends from -18 to 50 days).         |
+--------+------------------+-------------------------------------------------+
| sn91bg | 'hsiao_phase'    | 1991bg model restricted to a phase similar to   |
|        |                  | Hsiao 3.0 (extends from -18 to 85 days).        |
+--------+------------------+-------------------------------------------------+
| sn91bg | 'full_phase'     | 1991bg model extending over to full available   |
|        |                  | phase range (-18 to 100 days).                  |
+--------+------------------+-------------------------------------------------+

Usage Example
-------------

>>> import sncosmo
>>> from matplotlib import pyplot as plt
>>>
>>> from phot_class import models
>>>
>>> # Make sncosmo aware of the 1991bg models
>>> models.register_sources()
>>>
>>> # Initialize a model
>>> source = sncosmo.get_source('sn91bg', version='phase_limited')
>>> model = sncosmo.Model(source=source)
>>>
>>> # Get more information on how the source determines flux
>>> print(help(source))
>>>
>>> # run the fit
>>> data = sncosmo.load_example_data()
>>> result, fitted_model = sncosmo.fit_lc(
>>>     data, model,
>>>     ['z', 't0', 'x0'],  # parameters of model to vary
>>>     bounds={'z':(0.3, 0.7)})  # bounds on parameters (if any)
>>>
>>> # Plot results
>>> fig = sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)
>>> plt.show()

Function Documentation
----------------------
"""

import sncosmo as _sncosmo

from ._sources import SN91bg, load_template


# noinspection PyPep8Naming, PyUnusedLocal
def _load_sn91bg(name=None, version='salt2_phase'):
    """Return a SN 1991bg-like source class for SNCosmo

    If ``version='full_phase'`` then return the full model. If ``version``
    is the name of a model registered with sncosmo, return the model with the
    the phase limited to not extend past that model's.

    Args:
         name   (None): A dummy argument for compatibility with SNCosmo
         version (str): The version of the template to load
    """

    if version == 'full_phase':
        source = SN91bg(-float('inf'), float('inf'))
        source.version = version
        return source

    try:
        source_name = version.rstrip('_phase')
        s = _sncosmo.Model(source_name)
        source = SN91bg(s.source.minphase(), s.source.maxphase())
        source.version = version
        return source

    except Exception:
        raise ValueError(
            f'No secondary model found for given version: {source_name}')


def register_sources(force=False):
    """Register SN 1991bg-like models with SNCosmo

    Versions include: 'salt2_phase', 'hsiao_phase', and 'full_phase'

    Args:
        force (bool): Whether to overwrite an already registered source
    """

    _sncosmo.register_loader(
        data_class=_sncosmo.Source,
        name='sn91bg',
        func=_load_sn91bg,
        version='salt2_phase',
        force=force)

    _sncosmo.register_loader(
        data_class=_sncosmo.Source,
        name='sn91bg',
        func=_load_sn91bg,
        version='hsiao_phase',
        force=force)

    _sncosmo.register_loader(
        data_class=_sncosmo.Source,
        name='sn91bg',
        func=_load_sn91bg,
        version='full_phase',
        force=force)
