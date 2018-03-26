#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This document downloads and sorts photometric data from SDSS"""

import os

from astropy.table import Table
import numpy as np
import requests
import tarfile

SDSS_URL = 'https://data.sdss.org/sas/dr10/boss/papers/supernova/'
FILE_DIR = os.path.dirname(os.path.realpath(__file__))
DEFAULT_SMP = os.path.join(FILE_DIR, 'SMP_Data')
FILT_INDICES = {'u': 0, 'g': 1, 'r': 2, 'i': 3, 'z': 4}


def download_data(out_dir=FILE_DIR):
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


def read_smp_table(CID, filter, smp_dir=DEFAULT_SMP):
    """Returns photometric data for a supernova candidate in a given band

    Args:
        CID     (str): The Candidate ID of the desired object
        filter  (str): The desired filter (u, g, r, i, z)
        smp_dir (str): The directory of SMP data if not ./SMP_Data

    Returns:
        An astropy table of photometric data for the given candidate ID
    """

    file_name = 'SMP_{:06d}.dat'.format(CID)
    file_path = os.path.join(smp_dir, file_name)
    all_data = Table.read(file_path, format='ascii')

    col_names = all_data.meta['comments'][-1].split()
    for i, name in enumerate(col_names):
        all_data['col{}'.format(i + 1)].name = name

    filt_index = FILT_INDICES[filter]
    indices = np.where(all_data['FILT'] == filt_index)
    indexed_data = all_data[indices]
    indexed_data.remove_column('FILT')
    return indexed_data


if __name__ == '__main__':
    download_data()
