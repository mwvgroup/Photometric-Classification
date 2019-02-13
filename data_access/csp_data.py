#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module provides access to data from the third data release of the
Carnegie Supernova Project (CSP).

For more information on CSP see: https://csp.obs.carnegiescience.edu
"""

from os import path as _path

from astropy.table import Table

from ._download_data import download_data
from ._utils import parse_snoopy_data

CSP_URL = 'https://csp.obs.carnegiescience.edu/data/'

# Define local paths of CSP data
DATA_DIR = _path.join(_path.dirname(_path.realpath(__file__)), 'csp_data')
PHOT_DIR = _path.join(DATA_DIR, 'DR3')                 # DR3 Light Curves
FILT_DIR = _path.join(DATA_DIR, 'CSP_filter_package')  # DR3 Light Curves
MASTER_PTH = _path.join(PHOT_DIR, 'tab1.dat')          # Master table

# Download data if it does not exist
download_data(
    base_url=CSP_URL,
    out_dir=DATA_DIR,
    remote_name=['CSP_Photometry_DR3.tgz', 'CSP_filter_package.tgz'],
    check_local_name=[PHOT_DIR, FILT_DIR])

master_table = Table.read(MASTER_PTH, format='ascii')
