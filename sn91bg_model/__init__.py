#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module defines a custom 91bg source for use with SNCosmo. The base
model is based on Nugent's 91bg model but is extended to the UV with synthetic
spectra from Hachinger et al. 2008.
"""

import os as _os

from ._combine_template import combine_sed_templates

FILE_DIR = _os.path.abspath(_os.path.dirname(__file__))
MODEL_SOURCE_DIR = _os.path.join(FILE_DIR, 'sed_templates')
COMPILED_MODEL_PATH = _os.path.join(FILE_DIR, 'complete_template.npy')

if not _os.path.exists(COMPILED_MODEL_PATH):
    print('Compiled 91bg template not found. Creating a new template.')
    combine_sed_templates(out_path=COMPILED_MODEL_PATH, in_dir=MODEL_SOURCE_DIR)
    print('Done')


from ._sncosmo_ported_model import SN91bgSource
