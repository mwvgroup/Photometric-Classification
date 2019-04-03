#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module reads SNANA simulation results data, transforms them into
an SNcosmo format data table with metadata: z, t0, x0, st, c;
It also writes the resulting tables and metadata into csv files.
"""

import os

import numpy as np
from astropy.io import ascii
from astropy.table import Table, vstack


def read_lc(number, DATA_DIR):
    """
    Read SNANA simulation results and transform to SNcosmo format
    
    Args:
        number     (int): cid of the SN
        DATA_DIR   (str): path of SNANA simulation light curve data directory
       
    return:
        A SNcosmo-format data table with metadata: z,t0,x0,st,c 
    
    """

    path = os.path.join(DATA_DIR, '91bg_SN{:06d}.DAT'.format(number))
    file = open(path)
    h_start = 0
    d_start1 = 0
    d_end1 = 0
    d_start2 = 0
    d_end2 = 0
    section_2 = False
    header = True
    DETECT = False
    meta_data = {}

    # This for-loop find the header_start, data_start and data_end
    for line in file:
        if header == True:
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
                meta_data['x1'] = line.split()[1]

            if line[13:18] == 'color':
                meta_data['c'] = line.split()[1]

        if line[:4] == 'OBS:':
            if section_2 == False:
                d_end1 = d_end1 + 1

            elif section_2 == True:
                d_end2 = d_end2 + 1

        if line[:9] == 'DETECTION':
            DETECT = True
            d_start2 = d_end1 + 1
            d_end2 = d_start2
            section_2 = True

    file.close()
    if DETECT == True:
        table1 = ascii.read(path, guess=False,
                            delimiter=' ',
                            header_start=h_start,
                            data_start=d_start1,
                            data_end=d_end1)

        table2 = ascii.read(path, guess=False,
                            delimiter=' ',
                            header_start=h_start,
                            data_start=d_start2,
                            data_end=d_end2)

        table = vstack([table1, table2])

    else:
        table = ascii.read(
            path,
            guess=False,
            delimiter=' ',
            header_start=h_start,
            data_start=d_start1,
            data_end=d_end1)

    # creat SNcosmo-format table
    flux = table['FLUXCAL']  # *10**(-11+0.4*table['ZPT'])
    fluxerr = table['FLUXCALERR']  # *10**(-11+0.4*table['ZPT'])
    SNcosmo_table = Table(
        [table['MJD'], ['sdss' + b for b in table['FLT']], flux,
         fluxerr, np.full(len(table), 27.5), np.full(len(table), 'ab')],
        names=['time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'])

    for s in meta_data:
        SNcosmo_table.meta[s] = meta_data[s]
    return SNcosmo_table


def write_lc(num, DATA_DIR, path):
    """
    Parse SNANA data table into SNcosmo table and write to a file
    
    Args:
        num       (int): how many light curves to parse
        DATA_DIR  (str): path of SNANA simulation light curve data directory
        path      (str): path to write output files

    """

    meta = {}
    cid = []
    for i in range(num):
        table = read_lc(i + 1, DATA_DIR)
        if i == 0:
            for par in table.meta:
                meta[par] = []

        for par in meta:
            meta[par].append(table.meta[par])

        cid.append(i + 1)
        path_this = os.path.join(path, 'sn91bg_{:05d}.csv'.format(i + 1))
        table.write(path_this, overwrite=True)

    metadata = Table(
        [cid, meta['z'], meta['t0'], meta['x0'], meta['x1'], meta['c']],
        names=['cid', 'z', 't0', 'x0', 'x1', 'c'])

    metadata.write(os.path.join(path, 'meta.csv'), overwrite=True)
