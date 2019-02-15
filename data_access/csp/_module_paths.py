#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""Download any data files that do not exist locally and define file paths to
the local data for use by the parent module.
"""

import os

from data_access._utils import download_data

# Define local paths of published CSP data
_file_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(_file_dir, 'data')
photometry_dir = os.path.join(data_dir, 'DR3')  # DR3 Light Curves
filter_dir = os.path.join(data_dir, 'CSP_filter_package')  # DR3 Light Curves
master_path = os.path.join(photometry_dir, 'tab1.dat')  # Master table

# Download data if it does not exist
_csp_url = 'https://csp.obs.carnegiescience.edu/data/'
_local_file_names = [photometry_dir, filter_dir]
_remote_file_names = ['CSP_Photometry_DR3.tgz', 'CSP_filter_package.tgz']

download_data(
    base_url=_csp_url,
    out_dir=data_dir,
    remote_name=_remote_file_names,
    check_local_name=_local_file_names)
