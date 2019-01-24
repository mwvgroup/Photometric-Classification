#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module downloads photometric data from the SDSS supernova survey to
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


def download_data(base_url, out_dir, remote_name, check_local_name=None):
    """Downloads data files from a given url and unzip if it is a .tar.gz

    If check_local_names is provided, check if <out_dir>/<check_local_name[i]>
    exists first and don't download the file if it does.

    Args:
        base_url               (str): Url to download files from
        out_dir                (str): Directory to save files into
        remote_name      (list[str]): Name of files to download
        check_local_name (list[str]): Names of file to check for
    """

    for i, f_name in enumerate(remote_name):
        out_path = os.path.join(out_dir, f_name)
        if check_local_name is not None:
            check_path = os.path.join(out_dir, check_local_name[i])
            if os.path.exists(check_path):
                continue

        print('downloading', f_name)
        url = requests.compat.urljoin(base_url, f_name)
        _download_file(url, out_path)
