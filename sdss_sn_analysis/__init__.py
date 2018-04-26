#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This package downloades supernova data from SDSS and plots the corresponding
light curves.
"""

from warnings import warn as _warn

from ._settings import *
from .download_data import download_sdss_data
from .lightcurves import cid_data
from .lightcurves import plot_cid_light_curve
from .lightcurves import plot_sn_light_curves

for path in (MASTER_PTH, SMP_DIR, SNANA_DIR):
    if not os.path.exists(path):
        _warn('Detected missing local data. '
              'Attempting auto download from SDSS...')
        download_sdss_data()
        break
