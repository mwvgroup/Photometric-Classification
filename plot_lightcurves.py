#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This document downloads and sorts photometric data from SDSS"""

import os
import shutil

from astropy.table import Table
from matplotlib import pyplot as plt
import numpy as np
import requests
import tarfile
from tqdm import tqdm

SDSS_URL = 'https://data.sdss.org/sas/dr10/boss/papers/supernova/'
FILE_DIR = os.path.dirname(os.path.realpath(__file__))
MASTR_PTH = os.path.join(FILE_DIR, 'master_table.txt')  # Path of master table
SMP_DIR = os.path.join(FILE_DIR, 'smp_data/')           # For SMP data files
SNANA_DIR = os.path.join(FILE_DIR, 'snana_data/')       # For SNANA data files

# Indices used by SDSS to label different photometric filters
FILT_INDICES = {'u': 0, 'g': 1, 'r': 2, 'i': 3, 'z': 4}


def download_data(verbose=False):
    """Downloads supernova data from SDSS

    Downloaded files:
        master_data.txt
        SMP_Data.tar.gz
        SDSS_dataRelease-snana.tar.gz

    Args:
        verbose (bool): Whether to output progress to console
    """

    # Define file_names and output paths of files to download
    data_files = [('master_data.txt', MASTR_PTH),
                  ('SMP_Data.tar.gz', SMP_DIR),
                  ('SDSS_dataRelease-snana.tar.gz', SNANA_DIR)]

    if os.path.exists('./.temp'):
        shutil.rmtree('./.temp')
    os.mkdir('./.temp')

    iterable = data_files
    if verbose:
        print('Downloading SDSS data...', flush=True)
        iterable = tqdm(iterable)

    for f_name, out_path in iterable:
        url = requests.compat.urljoin(SDSS_URL, f_name)
        response = requests.get(url)
        response.raise_for_status()

        # Download data to a temporary path
        temp_path = os.path.join('./.temp', f_name)
        with open(temp_path, 'wb') as ofile:
            ofile.write(response.content)

        # Unzip file if it is an archive
        if temp_path.endswith(".tar.gz"):
            with tarfile.open(temp_path,"r:gz") as data:
                data.extractall('./.temp')

            os.remove(temp_path)
            file_info = next(os.walk('./.temp'))
            temp_path = os.path.join(file_info[0], file_info[1][0])

        shutil.move(temp_path, out_path)

    os.rmdir('./.temp')


def read_smp_table(cid, filt_name):
    """Returns photometric data for a supernova candidate in a given filter

    Args:
        cid       (str): The Candidate ID of the desired object
        filt_name (str): The desired filter (u, g, r, i, z)

    Returns:
        An astropy table of photometric data for the given candidate ID
    """

    file_name = 'SMP_{:06d}.dat'.format(cid)
    file_path = os.path.join(SMP_DIR, file_name)
    all_data = Table.read(file_path, format='ascii')

    col_names = all_data.meta['comments'][-1].split()
    for i, name in enumerate(col_names):
        all_data['col{}'.format(i + 1)].name = name

    filt_index = FILT_INDICES[filt_name]
    indices = np.where(all_data['FILT'] == filt_index)
    indexed_data = all_data[indices]
    indexed_data.remove_column('FILT')

    return indexed_data


def plot_all_cid_lightcurves(out_dir, verbose=False):
    """Create plots of light curves for all objects with locally available data

    Args:
        out_dir (str): The directory to write output plots to
        verbose (bool): Whether to output progress to console
    """

    out_path = os.path.join(out_dir, '{:06d}.pdf')
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    master_table = Table.read(MASTR_PTH, format='ascii')
    iterable = master_table['CID']
    if verbose:
        print('Plotting light curves:', flush=True)
        iterable = tqdm(iterable)

    for cid in iterable:
        for filter in ['u', 'g', 'r', 'i', 'z']:
            data_table = read_smp_table(cid, filter)
            plt.scatter(data_table['MJD'], data_table['MAG'], label=None, s=8)
            plt.plot(data_table['MJD'], data_table['MAG'],
                     label=filter, alpha=0.4)

        plt.title('CID: {}'.format(cid))
        plt.ylabel('Magnitude (mag)')
        plt.xlabel('Date (MJD)')
        plt.legend()
        plt.tight_layout()

        plt.savefig(out_path.format(cid), format='pdf')
        plt.clf()


def main():
    """Download SDSS supernova data and plot all light curves

    If SDSS data is not found locally, downloaded it to the destinations
    specified by global variables. Using this data, generate light curve plots
    for each supernova and save them to './light_curves'.
    """

    for path in (MASTR_PTH, SMP_DIR, SNANA_DIR):
        if not os.path.exists(path):
            download_data(True)
            break

    plot_all_cid_lightcurves('./light_curves', verbose=True)


if __name__ == '__main__':
    main()
