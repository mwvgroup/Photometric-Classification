#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This document downloads and sorts photometric data from SDSS"""

import os

import requests
import tarfile

SDSS_URL = 'https://data.sdss.org/sas/dr10/boss/papers/supernova/'
FILE_DIR = os.path.dirname(os.path.realpath(__file__))


def download_data(out_dir):
    """Downloads the SDSS supernova data

    Downloaded files:
        master_data.txt
        SMP_Data.tar.gz

    Args:
        out_dir (str): The directory where downloaded files are written
    """

    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    data_files = ['master_data.txt', 'SMP_Data.tar.gz']
    for fname in data_files:
        url = requests.compat.urljoin(SDSS_URL, fname)
        response = requests.get(url)
        response.raise_for_status()

        path = os.path.join(out_dir, fname)
        with open(path, 'wb') as ofile:
            ofile.write(response.content)

        if path.endswith("tar.gz"):
            with tarfile.open(path, "r:gz") as data:
                data.extractall(out_dir)

            os.remove(path)


if __name__ == '__main__':
    download_data('.')
