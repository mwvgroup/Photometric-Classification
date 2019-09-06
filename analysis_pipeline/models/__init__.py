#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module provides a custom SNCosmo source for modeling 91bg-like
supernovae. The model is based on the 91bg template from Nugent et al. 2002
but is extended into the UV.

For more information on the Nugent template see:
    https://iopscience.iop.org/article/10.1086/341707
"""


# noinspection PyPep8Naming, PyUnusedLocal
def SN91bg(name=None, version='salt2_phase'):
    """Return a version a SN 1991bg-like model for SNCosmo

    Versions include: ``salt2_phase``, ``color_interpolation``

    Args:
         name (None): Dummy arg for compatibility with SNCosmo
         version (str): The version of the template to load
             (Default: salt2_phase)
    """

    from . import sources

    if version == 'salt2_phase':
        return sources.Salt2Phase()

    if version == 'color_interpolation':
        return sources.ColorInterpolation()

    else:
        raise ValueError(f"Unidentified version: '{version}'.")


def register_sources(force=False):
    """Register SN 1991bg-like models with SNCosmo

    Versions include: ``salt2_phase``, ``color_interpolation``

    Args:
        force (bool): Whether to overwrite an already registered source
    """

    import sncosmo

    sncosmo.register_loader(
        data_class=sncosmo.Source,
        name='sn91bg',
        func=SN91bg,
        version='salt2_phase',
        force=force)

    sncosmo.register_loader(
        data_class=sncosmo.Source,
        name='sn91bg',
        func=SN91bg,
        version='color_interpolation',
        force=force)
