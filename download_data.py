#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This document defines functions for downloading SDSS photometric data."""

import os
import shutil

import requests
import tarfile

from _path_settings import *


def _download_file(f_name, out_path):
    """Download a specified file to a given output path

    Any top level .tar.gz archives will be automatically unzipped.

    Args:
        f_name   (str): The name of the file to download
        out_path (str): The path where the downloaded file should be written
    """

    url = requests.compat.urljoin(SDSS_URL, f_name)
    response = requests.get(url)
    response.raise_for_status()

    # Write data to a temporary path
    temp_path = os.path.join(TEMP_DIR, f_name)
    with open(temp_path, 'wb') as ofile:
        ofile.write(response.content)

    # Unzip file if it is an archive
    if temp_path.endswith(".tar.gz"):
        with tarfile.open(temp_path, "r:gz") as data:
            data.extractall(TEMP_DIR)

        os.remove(temp_path)
        new_file = os.listdir(TEMP_DIR)[0]
        temp_path = os.path.join(TEMP_DIR, new_file)

    shutil.move(temp_path, out_path)


def download_sdss_data():
    """Downloads supernova data from SDSS

    Downloaded files include:
        master_data.txt
        SMP_Data.tar.gz
        SDSS_dataRelease-snana.tar.gz
    """

    if os.path.exists(DATA_DIR):
        shutil.rmtree(DATA_DIR)

    os.mkdir(DATA_DIR)
    os.mkdir(TEMP_DIR)

    # Define file_names and output paths of files to download
    data_files = [('master_data.txt', MASTER_PTH),
                  ('SMP_Data.tar.gz', SMP_DIR),
                  ('SDSS_dataRelease-snana.tar.gz', SNANA_DIR)]

    for f_name, out_path in data_files:
        _download_file(f_name, out_path)

    shutil.rmtree(TEMP_DIR)


def gen_cid_lists():
    """Generates input CID lists for SNANA

    Files are written to DATA_DIR/cids.list .
    """

    cids = sorted([int(x[4:10]) for x in os.listdir(SMP_DIR)])
    out_str = ' '.join([str(x) for x in cids])
    out_path = os.path.join(DATA_DIR, 'cids.list')
    with open(out_path, 'w') as ofile:
        ofile.write(out_str)


if __name__ == '__main__':
    download_sdss_data()
    gen_cid_lists()
