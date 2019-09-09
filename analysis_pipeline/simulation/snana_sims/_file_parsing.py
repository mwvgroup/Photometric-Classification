# !/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

""" This module helps read SNANA simulation results data, transform into SNcosmo-format
    data table with metadata: z,t0,x0,st,c;
    Write all tables and metadata into csv files.
"""

import os
from pathlib import Path

import numpy as np
from astropy.io import ascii
from astropy.table import Table, vstack


def create_file_path(cid, data_dir):
    """Construct the file path for an SNANA simulation result

    Args:
        cid      (int): cid of the SN
        data_dir (str): path of SNANA simulation light curve data directory

    return:
        A Path object representing simulation results for the given cid
    """

    return Path(data_dir) / f'91bg_SN{cid:06d}.DAT'


def read_lc(path):
    """Read SNANA simulation results and transform to sncosmo format
    
    Args:
        path (str): Path of the file to read
        DATA_DIR   (str): path of SNANA simulation light curve data directory
       
    return:
        A SNcosmo-format data table with metadata: z,t0,x0,st,c
    """

    h_start = 0
    d_start1 = 0
    d_end1 = 0
    d_start2 = 0
    d_end2 = 0
    section_2 = False
    header = True
    DETECT = False
    meta_data = {}

    with open(path) as infile:
        for line in infile:  # Find the header_start, data_start and data_end
            if header:
                if line[:7] == 'VARLIST':
                    d_start1 = h_start + 1
                    d_end1 = d_start1
                    header = False

                elif line[:1] not in ('#' + '\n' + ' '):
                    h_start = h_start + 1

                # 5 fit parameters: [z,t0,x0,st,c]
                if line[:16] == 'SIM_REDSHIFT_CMB':
                    meta_data['z'] = line.split()[1]

                if line[:11] == 'SIM_PEAKMJD':
                    meta_data['t0'] = line.split()[1]

                if line[:11] == 'SIM_SALT2x0':
                    meta_data['x0'] = line.split()[1]

                if line[13:20] == 'stretch':
                    meta_data['st'] = line.split()[1]

                if line[13:18] == 'color':
                    meta_data['c'] = line.split()[1]

            if line[:4] == 'OBS:':
                if not section_2:
                    d_end1 = d_end1 + 1

                else:
                    d_end2 = d_end2 + 1

            if line[:9] == 'DETECTION':
                DETECT = True
                d_start2 = d_end1 + 1
                d_end2 = d_start2
                section_2 = True

    if DETECT:
        table1 = ascii.read(
            path, guess=False, delimiter=' ', header_start=h_start,
            data_start=d_start1, data_end=d_end1)

        table2 = ascii.read(
            path, guess=False, delimiter=' ', header_start=h_start,
            data_start=d_start2, data_end=d_end2)

        table = vstack([table1, table2])

    else:
        table = ascii.read(path, guess=False, delimiter=' ',
                           header_start=h_start, data_start=d_start1,
                           data_end=d_end1)

    # creaet sncosmo formated table
    mjd = table['MJD']
    bands = ['sdss' + b for b in table['FLT']]
    flux = table['FLUXCAL']  # * 10 ** (-11 + 0.4 * table['ZPT'])
    fluxerr = table['FLUXCALERR']  # * 10 ** (-11 + 0.4 * table['ZPT'])
    zp = np.full(len(mjd), 27.5)
    zp_sys = np.full(len(mjd), 'ab')

    sncosmo_table = Table(
        data=[mjd, bands, flux, fluxerr, zp, zp_sys],
        names=['time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'])

    sncosmo_table.meta.update(meta_data)
    return sncosmo_table


def write_lc(num, out_dir, path):
    """Parse SNANA data table into sncosmo table and write to a file
    
    Args:
        num     (int): how many light curves to parse
        out_dir (str): path of SNANA simulation light curve data directory
        path    (str): path to write output files

    """

    meta = {}
    cid = []
    for i in range(num):
        table = read_lc(i + 1, out_dir)

        if i == 0:
            for par in table.meta:
                meta[par] = []

        for par in meta:
            meta[par].append(table.meta[par])

        cid.append(i + 1)
        table.write(os.path.join(path, 'sn91bg_{:05d}.csv'.format(i + 1)),
                    overwrite=True)

    metadata = Table(
        [cid, meta['z'], meta['t0'], meta['x0'], meta['st'], meta['c']],
        names=['cid', 'z', 't0', 'x0', 'st', 'c'])

    metadata.write(os.path.join(path, 'meta.csv'), overwrite=True)
