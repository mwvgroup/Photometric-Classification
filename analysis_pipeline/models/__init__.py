#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module provides a SNCosmo Source class for modeling 91bg-like
supernovae. The model is based on the 91bg template from Nugent et al.
2002 but is extended into the ultra-violet.

For more information on the Nugent template see:
    https://iopscience.iop.org/article/10.1086/341707
"""

from ._sources import load_template


# noinspection PyPep8Naming, PyUnusedLocal
def SN91bg(name=None, version='phase_limited'):
    """Return a version a SN 1991bg-like model for SNCosmo

    Versions include: 'phase_limited', 'color_interpolation'

    Args:
         name   (None): A dummy argument for compatibility with SNCosmo
         version (str): The version of the template to load
    """

    from . import _sources

    if version == 'phase_limited':
        return _sources.PhaseLimited()

    if version == 'color_interpolation':
        return _sources.ColorInterpolation()

    else:
        raise ValueError(f"Unidentified version: '{version}'.")


def register_sources(force=False):
    """Register SN 1991bg-like models with SNCosmo

    Versions include: 'phase_limited', 'color_interpolation'

    Args:
        force (bool): Whether to overwrite an already registered source
    """

    import sncosmo

    sncosmo.register_loader(
        data_class=sncosmo.Source,
        name='sn91bg',
        func=SN91bg,
        version='phase_limited',
        force=force)

    sncosmo.register_loader(
        data_class=sncosmo.Source,
        name='sn91bg',
        func=SN91bg,
        version='color_interpolation',
        force=force)
