#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module provides a custom SNCosmo source for modeling 91bg-like
supernovae. The model is based on the 91bg template from Nugent et al. 2002
but is extended into the UV.

For more information on the Nugent template see:
    https://iopscience.iop.org/article/10.1086/341707
"""

from . import _color_interpolation, _salt2_phase


def SN91bgSource(version='salt2_phase'):
    if version == 'salt2_phase':
        return _salt2_phase.SN91bgSource()

    if version == 'color_interpolation':
        return _color_interpolation.SN91bgSource()

    else:
        raise Exception(f"No 'SN91bg' with version='{version}' in registry. "
                        f"Registered versions: "
                        f"'salt2_phase', 'color_interpolation'")
