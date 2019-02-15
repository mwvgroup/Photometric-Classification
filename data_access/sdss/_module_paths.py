#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""Download any data files that do not exist locally and define file paths to
the local data for use by the parent module.
"""

import os
from itertools import product

from data_access._utils import download_data, register_filter

# Define local paths of published SDSS data
_file_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(_file_dir, 'data')
filter_dir = os.path.join(data_dir, 'doi_2010_filters/')       # SDSS filters
snana_dir = os.path.join(data_dir, 'SDSS_dataRelease-snana/')  # SNANA files
master_table_path = os.path.join(data_dir, 'master_data.txt')  # Master table
smp_dir = os.path.join(data_dir, 'SMP_Data/')                  # SMP data files

# Download light curve data if it does not exist locally
_sdss_url = 'https://data.sdss.org/sas/dr10/boss/papers/supernova/'
_local_file_names = [master_table_path, smp_dir, snana_dir]
_remote_file_names = ['master_data.txt',
                      'SMP_Data.tar.gz',
                      'SDSS_dataRelease-snana.tar.gz']

download_data(
    base_url=_sdss_url,
    out_dir=data_dir,
    remote_name=_remote_file_names,
    check_local_name=_local_file_names
)

# Download light curve data if it does not exist locally
_filt_url = 'http://www.ioa.s.u-tokyo.ac.jp/~doi/sdss/'
_local_filt_names = ['{}{}.dat'.format(a, b) for a, b in
                     product('ugriz', '123456')]
download_data(
    base_url=_filt_url,
    out_dir=filter_dir,
    remote_name=_local_filt_names,
    check_local_name=_local_filt_names
)

# Register filters if not already registered
for _filter_path in _local_filt_names:
    fpath = os.path.join(filter_dir, _filter_path)
    register_filter(fpath, 'doi_2010_' + _filter_path[:2])
