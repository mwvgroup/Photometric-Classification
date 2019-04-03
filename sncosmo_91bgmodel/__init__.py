#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This package provides a custom SNCosmo source for modeling 91bg-like
supernovae. The model is based on the 91bg template from Nugent et al. 2002
but is extended into the UV.

For more information on the Nugent template see:
    https://iopscience.iop.org/article/10.1086/341707

To use:
    from sncosmo_91bgmodel import SN91bgSource
"""

from ._sncosmo_91bgmodel import SN91bgSource