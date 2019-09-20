#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``models`` module defines a ``Source`` class for modeling
91bg-like supernovae with ``sncosmo``.

Available Model Versions
------------------------

There are multiple versions of the 91bg model included in this package. For
information on how each version determines the flux for a given set of
parameters, see the documentation for the given source class. It should **NOT**
be assumed that flux is predicted the same way for different model versions.

+--------+---------------------------+---------------------------------------------------------+
| Name   | Version                   | Description                                             |
+========+===========================+=========================================================+
| sn91bg | 'phase_limited' (Default) | 1991bg model restricted to a phase similar to Salt 2.4. |
|        |                           | This model version extends from -18 to 50 days.         |
+--------+---------------------------+---------------------------------------------------------+
| sn91bg | 'full_phase'              | 1991bg model extending over to full phase range.        |
|        |                           | This model version extends from -18 to 100 days.        |
+--------+---------------------------+---------------------------------------------------------+

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
>>> # Initialize a model where the version is the model mass
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

Module Level Functions
----------------------
"""

import sncosmo as _sncosmo

from ._sources import SN91bg, load_template


# noinspection PyPep8Naming, PyUnusedLocal
def _load_sn91bg(name=None, version='phase_limited'):
    """Return a SN 1991bg-like source class for SNCosmo

    Versions include: 'phase_limited', 'full_phase'

    Args:
         name   (None): A dummy argument for compatibility with SNCosmo
         version (str): The version of the template to load
    """

    if version in 'phase_limited':
        model = SN91bg()
        model.version = version
        return model

    if version == 'full_phase':
        model = SN91bg(-float('inf'), float('inf'))
        model.version = version
        return model


def register_sources(force=False):
    """Register SN 1991bg-like models with SNCosmo

    Versions include: 'phase_limited', 'full_phase'

    Args:
        force (bool): Whether to overwrite an already registered source
    """

    _sncosmo.register_loader(
        data_class=_sncosmo.Source,
        name='sn91bg',
        func=_load_sn91bg,
        version='phase_limited',
        force=force)

    _sncosmo.register_loader(
        data_class=_sncosmo.Source,
        name='sn91bg',
        func=_load_sn91bg,
        version='full_phase',
        force=force)
