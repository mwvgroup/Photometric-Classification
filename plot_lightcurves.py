#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This document defines function for plotting the light curves of objects from
the SDSS-II SN Catalog Data Release.
"""

import os

from matplotlib import pyplot as plt
import numpy as np
from tqdm import tqdm

from parse_sn_data import MASTER_TABLE, get_cid_data


def _plot_master_row(cid_row, out_dir='./'):
    """Plot the light curve for a given entry in the master table

    Args:
        cid_row (row): A row from the SMP master table
        out_dir (str): The directory where output plots are written
    """

    for filter_name in ['u', 'g', 'r', 'i', 'z']:
        data_table = get_cid_data(cid_row['CID'], filter_name)
        plt.scatter(data_table['MJD'], data_table['MAG'], label=None, s=8)
        plt.plot(data_table['MJD'],
                 data_table['MAG'],
                 label=filter_name,
                 alpha=0.4)

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

    index = np.where(MASTER_TABLE['CID'] == f'{cid:06d}')[0][0]
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
