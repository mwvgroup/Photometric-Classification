#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This document defines functions for downloading SDSS photometric data."""

import os
import shutil

import requests
import tarfile

from . import _settings


def _download_file(f_name, out_path):
    """Download a specified file to a given output path

    Any top level .tar.gz archives will be automatically unzipped.

    Args:
        f_name   (str): The name of the file to download
        out_path (str): The path where the downloaded file should be written
    """

    url = requests.compat.urljoin(_settings.SDSS_URL, f_name)
    response = requests.get(url)
    response.raise_for_status()

    # Write data to a temporary path
    temp_path = os.path.join(_settings.TEMP_DIR, f_name)
    with open(temp_path, 'wb') as ofile:
        ofile.write(response.content)

    # Unzip file if it is an archive
    if temp_path.endswith(".tar.gz"):
        with tarfile.open(temp_path, "r:gz") as data:
            data.extractall(_settings.TEMP_DIR)

        os.remove(temp_path)
        new_file = os.listdir(_settings.TEMP_DIR)[0]
        temp_path = os.path.join(_settings.TEMP_DIR, new_file)

    shutil.move(temp_path, out_path)


def download_sdss_data():
    """Downloads supernova data from SDSS

    Downloaded files include:
        master_data.txt
        SMP_Data.tar.gz
        SDSS_dataRelease-snana.tar.gz
    """

    if os.path.exists(_settings.TEMP_DIR):
        shutil.rmtree(_settings.TEMP_DIR)

    os.mkdir(_settings.TEMP_DIR)

    # Define file_names and output paths of files to download
    data_files = [('master_data.txt', _settings.MASTER_PTH),
                  ('SMP_Data.tar.gz', _settings.SMP_DIR),
                  ('SDSS_dataRelease-snana.tar.gz', _settings.SNANA_DIR)]

    for f_name, out_path in data_files:
        _download_file(f_name, out_path)

    shutil.rmtree(_settings.TEMP_DIR)
