#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module provides access to data from the year three DES supernova
cosmology paper.
"""

from os import path as _path

from _download_data import download_data

DES_URL = 'http://desdr-server.ncsa.illinois.edu/despublic/sn_files/y3/tar_files/'
DATA_DIR = _path.join(_path.dirname(_path.realpath(__file__)), 'des_data')
FILT_DIR = _path.join(DATA_DIR, '01-FILTERS')  # Master table path
PHOT_DIR = _path.join(DATA_DIR, '02-DATA_PHOTOMETRY')  # SMP data files

# Download data if it does not exist
download_data(DES_URL, './des_data', ['01-FILTERS.tar.gz',
                                      '02-DATA_PHOTOMETRY.tar.gz'])
