#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module formats SDSS data for use with snpy and fits the SDSS light
curves"""

import os
import logging

import snpy
from astropy.table import Table
from matplotlib import pyplot as plt

from parse_sn_data import MASTER_TABLE, get_cid_data


def write_snoopy_sdss_file(out_path, **kwargs):
    """Creates a snoopy input file <out_path>/<name>.txt

    Args:
        out_path         (str): Where to write the input file
        name             (str): The name of the target
        redshift       (float): The redshift of the target
        ra             (float): The ra of the target
        dec            (float): The dec of the target
        <u, g, r, i, z> (dict): Dict with date ('MJD'), magnitude ('MAG'),
                                  and magnitude error ('MERR')
    """

    file_text = '{} {} {} {}\n'.format(kwargs["name"],
                                       kwargs["redshift"],
                                       kwargs["ra"],
                                       kwargs["dec"])

    for band in 'ugriz':
        filter_name = band + '_s'
        if band not in kwargs:
            continue

        file_text += 'filter {}\n'.format(filter_name)
        band_data = kwargs[band]
        for mjd, mag, error in \
                zip(band_data['MJD'], band_data['MAG'], band_data['MERR']):
            file_text += '{} {} {}\n'.format(mjd, mag, error)

    with open(out_path, 'w') as ofile:
        ofile.write(file_text)


def create_snoopy_inputs(out_dir, iter_paths=False):
    """Create snoopy input files for all SMP targets

    Files are named as <out_dir>/<target cid>.txt

    Args:
        out_dir     (str): The directory to write files to
        iter_paths (bool): Whether to yield output file paths
    """

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    for cid in MASTER_TABLE['CID']:
        out_path = os.path.join(out_dir, '{}.txt'.format(cid))

        cid_data = get_cid_data(cid)
        u = cid_data[cid_data['FILT'] == 0]
        g = cid_data[cid_data['FILT'] == 1]
        r = cid_data[cid_data['FILT'] == 2]
        i = cid_data[cid_data['FILT'] == 3]
        z = cid_data[cid_data['FILT'] == 4]

        redshift = u.meta['redshift']
        ra = u.meta['ra']
        dec = u.meta['dec']

        write_snoopy_sdss_file(
            out_path,
            name=cid,
            redshift=redshift,
            ra=ra,
            dec=dec,
            u=u,
            g=g,
            r=r,
            i=i,
            z=z
        )

        if iter_paths:
            yield out_path


def fit_sdss_data(out_dir, model_name='color_model'):

    # Store warnings and errors
    logging.captureWarnings(True)
    log_path = os.path.join(out_dir, 'summary.log')

    if os.path.exists(log_path):
        os.remove(log_path)

    logging.basicConfig(filename=log_path, level=logging.WARNING)

    # Table to store summary of fit results
    names = ['cid', 'class', 'z', 'fit_z', 't0', 'x0', 'x1', 'c', 'chi', 'dof', 'message']
    out_table = Table(names=names, dtype=[object for _ in names])
    out_path = os.path.join(out_dir, 'summary.csv')
    out_table.write(out_path, overwrite=True)

    # Make output directories
    inputs_dir = os.path.join(out_dir, 'inputs')
    fig_dir = os.path.join(out_dir, 'figs')
    for directory in (inputs_dir, fig_dir):
        if not os.path.exists(directory):
            os.makedirs(directory)

    # Fit smp data with snoopy
    for path in create_snoopy_inputs(inputs_dir, iter_paths=True):
        fig_path = path.replace('inputs', 'figs').replace('.txt', '.pdf')

        try:
            s = snpy.get_sn(path)
            s.choose_model(model_name)
            s.fit()
            s.plot()

            plt.title(path)
            plt.savefig(fig_path)

        except (ValueError, RuntimeError, TypeError) as e:
            logging.error('{}: {}\n'.format(path, e))


if __name__ == '__main__':
    fit_sdss_data('./snoopy_results')
