#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This dscript downloads photometric data from the SDSS supernova survey to
the directory ./data .
"""

from __future__ import print_function

import os
import tarfile

import requests


def _download_file(url, out_path):
    """Download a specified file to a given output path

    Any top level .tar.gz archives will be automatically unzipped.

    Args:
        url      (str): URL of the file to download
        out_path (str): The path where the downloaded file should be written
    """

    # Make temporary file path
    if os.path.isdir(out_path):
        temp_path = os.path.join(out_path, '.temp')
        out_dir = out_path

    else:
        temp_path = out_path + '.temp'
        out_dir = os.path.dirname(out_path)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    # Download data
    response = requests.get(url)
    response.raise_for_status()
    with open(temp_path, 'wb') as ofile:
        ofile.write(response.content)

    # Unzip file if its an archive
    if url.endswith(".tar.gz"):
        with tarfile.open(temp_path, "r:gz") as data:
            data.extractall(out_dir)

        os.remove(temp_path)

    else:
        os.rename(temp_path, out_path)


def download_data(base_url, out_dir, data_files):
    """Downloads data files from a given url

    Args:
        base_url         (str): Url to download files from
        out_dir          (str): Director to save files into
        data_files (list[str]): Name of files to download
    """

    for f_name in data_files:
        out_path = os.path.join(out_dir, f_name)
        if not os.path.exists(out_path):
            print('downloading', f_name)
            url = requests.compat.urljoin(base_url, f_name)
            _download_file(url, out_path)
