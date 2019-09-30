#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``models`` module defines custom ``Source`` classes for modeling
supernovae with ``sncosmo``.

Available Model Versions
------------------------

There are multiple versions of a SN 1991bg-like model included in this package.
For information on how each version determines the flux for a given set of
parameters, see the documentation for the given source class. A modified
version of the Hsiao model is also available, which is the same as sncosmo's
built in ``hsiao`` model (version 3.0) but with an added stretch parameter.

+----------+---------------+-------------------------------------------------+
| Name     | Version       | Description                                     |
+==========+=================+===============================================+
| sn91bg   | 'salt2_phase' | 1991bg model restricted to a phase similar to   |
|          | (Default)     | Salt 2.4 (extends from -18 to 50 days).         |
+----------+---------------+-------------------------------------------------+
| sn91bg   | 'hsiao_phase' | 1991bg model restricted to a phase similar to   |
|          |               | Hsiao 3.0 (extends from -18 to 85 days).        |
+----------+---------------+-------------------------------------------------+
| sn91bg   | 'full_phase'  | 1991bg model extending over to full available   |
|          |               | phase range (-18 to 100 days).                  |
+----------+---------------+-------------------------------------------------+
| hsiao_x1 | '3.0.x1'      | The same as sncosmo's built in ``hsiao``        |
|          |               | model (version 3.0) but with an added stretch   |
|          |               | parameter and limited to -18 to 85 days.        |
+----------+---------------+-------------------------------------------------+

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
def _load_source(name, version):
    """Return an instantiated Source class

    Args:
         name    (str): The name of the source to load
         version (str): The version of the source to load
    """

    if name == 'sn91bg':
        if version == 'full_phase':
            source = SN91bg(-float('inf'), float('inf'))
            source.version = version
            return source

        source_name = version.rstrip('_phase')
        s = _sncosmo.Model(source_name)
        source = SN91bg(s.source.minphase(), s.source.maxphase())
        source.version = version
        return source

    elif name == 'hsiao_x1':
        return _sources.HsiaoStretch()

    else:
        raise ValueError(f'Unknown source {name} / {version}')


def register_sources(force=False):
    """Register custom models with SNCosmo

    See the documentation of the parent module for more information on the
    available models.

    Args:
        force (bool): Whether to overwrite an already registered source
    """

    _sncosmo.register_loader(
        data_class=_sncosmo.Source,
        name='sn91bg',
        func=_load_source,
        version='salt2_phase',
        force=force)

    _sncosmo.register_loader(
        data_class=_sncosmo.Source,
        name='sn91bg',
        func=_load_source,
        version='hsiao_phase',
        force=force)

    _sncosmo.register_loader(
        data_class=_sncosmo.Source,
        name='sn91bg',
        func=_load_source,
        version='full_phase',
        force=force)

    _sncosmo.register_loader(
        data_class=_sncosmo.Source,
        name='hsiao_x1',
        func=_load_source,
        version='3.0.x1',
        force=force)
