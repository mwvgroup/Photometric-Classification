#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This document defines function for plotting the light curves of objects from
the SDSS-II SN Catalog Data Release.
"""

import os

from astropy.table import Table
from matplotlib import pyplot as plt
import numpy as np
from tqdm import tqdm

from _path_settings import *

# Indices used by SDSS to label different photometric filters
FILT_INDICES = {'u': 0, 'g': 1, 'r': 2, 'i': 3, 'z': 4}
MASTER_TABLE = Table.read(MASTER_PTH, format='ascii')


def cid_data(cid, filt_name):
    """Returns photometric data for a supernova candidate in a given filter

    Args:
        cid       (int): The Candidate ID of the desired object
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


def _plot_master_row(cid_row, out_dir='./'):
    """Plot the light curve for a given entry in the master table

    Args:
        cid_row (row): A row from the SMP master table
        out_dir (str): The directory where output plots are written
    """

    for filter_name in ['u', 'g', 'r', 'i', 'z']:
        data_table = cid_data(cid_row['CID'], filter_name)
        plt.scatter(data_table['MJD'], data_table['MAG'], label=None, s=8)
        plt.plot(data_table['MJD'], data_table['MAG'],
                 label=filter_name, alpha=0.4)

    if cid_row['IAUName']:
        title = 'SN {} (CID {}, Type {})'.format(
            cid_row['IAUName'], cid_row['CID'], cid_row['Classification'])

    else:
        title = ' CID: {}'.format(cid_row['CID'])

    plt.title(title)
    plt.ylabel('Magnitude (mag)')
    plt.xlabel('Date (MJD)')
    plt.legend()
    plt.tight_layout()

    out_path = os.path.join(out_dir, '{:06d}.pdf')
    plt.savefig(out_path.format(cid_row['CID']), format='pdf')
    plt.clf()


def plot_cid_light_curve(cid, out_dir):
    """Plot the light curve for a given CID""

    Args:
        cid     (int): The id number of the desired SMP object
        out_dir (str): The directory where output plots are written
    """

    cid_str = '{:06d}'.format(cid)
    index = np.where(MASTER_TABLE['CID'] == cid_str)[0][0]
    _plot_master_row(MASTER_TABLE[index], out_dir)


def plot_sn_light_curves(out_dir, verbose):
    """Create plots of light curves for all objects with locally available data

    Args:
        out_dir  (str): The directory where output plots are written
        verbose (bool): Whether to output progress to console
    """

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    is_classified = np.logical_not(MASTER_TABLE['IAUName'].mask)
    classified_objects = MASTER_TABLE[is_classified]

    if verbose:
        print('Plotting SN light curves ...', flush=True)
        classified_objects = tqdm(classified_objects)

    for row in classified_objects:
        _plot_master_row(row, out_dir)


if __name__ == '__main__':
    plot_sn_light_curves('./light_curves', True)
