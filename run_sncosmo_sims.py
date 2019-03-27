# !/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""
This module generates light curves using SNcosmo and write simulation results
to csv files
"""

import os
import sys

import numpy as np
import sncosmo
from astropy.table import Table
from numpy.random import multivariate_normal, uniform

sys.path.insert(0, '../')
from analysis_pipeline import SN91bgSource

FILE_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(FILE_DIR, 'sncosmo_sims')
if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR, exist_ok=True)


def generate_lc(num, model):
    """
    Generate a series of light curves with specified model and write to
    csv files

    Args:
        num   (int): Number of light curves to generate
        model (sncosmo model): The model to be used
    """

    # generate redshifts
    area = 1.
    tmin = 53615.0
    tmax = 53710.0
    zmax = 0.45
    redshifts = list(sncosmo.zdist(0., zmax, time=(tmax - tmin), area=area))
    while len(redshifts) < num:
        area = area * 2
        redshifts = list(
            sncosmo.zdist(0., zmax, time=(tmax - tmin), area=area))

    redshifts = redshifts[0:num]

    # randomly generate peakmjd
    peakmjd = [uniform(tmin, tmax) for i in range(num)]

    # randomly generate stretch and color in 2D normal distribution
    # (covariance matrix included)
    x1 = np.zeros(num)
    c = np.zeros(num)
    for i in range(num):
        stretch, color = multivariate_normal(
            [0.975, 0.557],
            [[0.096 ** 2, -0.0110263], [-0.0110263, 0.175 ** 2]])

        while not (stretch < 1.25 and stretch > 0.65 and color < 1. and color > 0.):
            stretch, color = multivariate_normal(
                [0.975, 0.557],
                [[0.096 ** 2, -0.0110263], [-0.0110263, 0.175 ** 2]])

        x1[i] = stretch
        c[i] = color

    # parameters to generate light curves
    params = [{'z': redshifts[i],
               'x0': 1e-14,
               't0': peakmjd[i],
               'x1': x1[i],
               'c': c[i]} for i in range(num)]

    date = [(t + np.linspace(-13, 100, 20)).tolist() * 5 for t in peakmjd]
    flt = []
    for b in 'ugriz':
        flt = flt + ['sdss' + b] * (int(len(date[0]) / 5))

    obs = [Table({'time': date[i],
                  'band': flt,
                  'gain': np.full(len(flt), 1.0),
                  'skynoise': np.zeros(len(flt)),
                  'zp': np.full(len(flt), 27.5),
                  'zpsys': np.full(len(flt), 'ab')}) for i in range(num)]

    # generate light curves
    lcs = [sncosmo.realize_lcs(obs[i], model, [params[i]])[0]
           for i in range(num)]

    # write light curve tables into csv files
    cid = []
    for i in range(len(lcs)):
        cid.append('{:05d}'.format(i + 1))
        lcs_path = os.path.join(OUT_DIR + f'/sn91bg_{i + 1:05d}.csv')
        lcs[i].write(lcs_path, overwrite=True)

    lc_dump = Table([cid, redshifts, peakmjd, x1, c],
                    names=['cid', 'z', 't0', 'stretch', 'color'])

    lc_dump_path = os.path.join(OUT_DIR, 'sn91bg_dump.csv')
    lc_dump.write(lc_dump_path, overwrite=True)


if __name__ == '__main__':
    model = sncosmo.Model(source=SN91bgSource())
    generate_lc(1000, model)
