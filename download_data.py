#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This document defines functions for downloading SDSS photometric data."""

import os

import requests
import tarfile

SDSS_URL = 'https://data.sdss.org/sas/dr10/boss/papers/supernova/'


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


def download_sdss_data(out_dir):
    """Downloads supernova data from SDSS

    Downloaded files include:
        master_data.txt
        SMP_Data.tar.gz
        SDSS_dataRelease-snana.tar.gz
    """

    # Define file_names and output paths of files to download
    master_path = os.path.join(out_dir, 'master_table.txt')  # Master table
    data_files = [('master_data.txt', master_path),
                  ('SMP_Data.tar.gz', out_dir),
                  ('SDSS_dataRelease-snana.tar.gz', out_dir)]

    # Download each file
    for f_name, out_path in data_files:
        print('downloading ', f_name)
        url = requests.compat.urljoin(SDSS_URL, f_name)
        _download_file(url, out_path)


if __name__ == '__main__':
    download_sdss_data('./data')
