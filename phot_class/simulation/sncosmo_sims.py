# !/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``simulation.sncsomo_sims`` module generates light-curves using sncosmo
and write models results to file.
"""

from pathlib import Path

import numpy as np
import sncosmo
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
from numpy.random import multivariate_normal

__all__ = ['AVG_COLOR', 'AVG_STRETCH', 'COVARIANCE', 'bg_stretch_color',
           'sim_bg_params', 'generate_lc']

AVG_COLOR = 0.557
AVG_STRETCH = 0.975
COVARIANCE = ((0.096 ** 2, -0.0110263), (-0.0110263, 0.175 ** 2))


def bg_stretch_color(
        size, avg_color=AVG_COLOR, avg_stretch=AVG_STRETCH,
        covariance=COVARIANCE, min_stretch=0.65, max_stretch=1.25,
        min_color=0, max_color=1):
    """Return random arrays for 91bg stretch and color

    Args:
        size          (int): Number of stretches / colors to return
        avg_color   (float): Average color value
        avg_stretch (float): Average stretch value
        covariance  (float): Covariance matrix for stretch / color
        min_stretch (float): Minimum stretch to return
        max_stretch (float): Maximum stretch to return
        min_color   (float): Minimum color to return
        max_color   (float): Maximum color to return

    Returns:
        An array of stretch values
        An array of color values
    """

    size = int(size)

    # Simulate stretch and color. Then check if any are out of our model bounds
    # and replace out of bound values with new values
    x1, c, out_of_bounds = np.ones(size), np.ones(size), np.full(size, True)
    while sum(out_of_bounds):
        x1[out_of_bounds], c[out_of_bounds] = multivariate_normal(
            mean=[avg_stretch, avg_color],
            cov=covariance,
            size=sum(out_of_bounds)).T

        out_of_bounds = ~(
                (min_stretch < x1)
                & (x1 < max_stretch)
                & (min_color < c)
                & (c < max_color)
        )

    return x1, c


def sim_bg_params(
        zmin, zmax, tmin, tmax, area=1., ratefunc=None, cosmo=None, **kwargs):
    """Simulate parameters for "observed" 91bg light curves

    Args:
        zmin (float): Minimum redshift value
        zmax (float): Maximum redshift value
        tmin (float): Minimum observation time in days
        tmax (float): Maximum observation time in days
        area (float): Sky area in deg2 to simulate. Determines number of SNe.
        ratefunc (callable):
            A callable that accepts a single float (redshift) and returns the
            comoving volumetric rate at each redshift in units of yr^-1 Mpc^-3.
            The default is a function that returns ``1.e-4``.
        cosmo (Cosmology):
            The cosmology used to determine volume. The default is a
            FlatLambdaCDM cosmology with ``Om0=0.3``, ``H0=70.0``.
        Any other kwargs for ``bg_stretch_color`` except ``size``

    Returns:
        A list of dictionaries with model parameters for each light-curve
    """

    ratefunc = ratefunc or (lambda z: 1.e-4)
    cosmo = cosmo or FlatLambdaCDM(H0=70.0, Om0=0.3)

    # Generated redshifts includes SN rate so we base the number of SNe off z
    redshifts = list(
        sncosmo.zdist(zmin, zmax, (tmax - tmin), area, ratefunc, cosmo))
    peakmjd = np.random.uniform(tmin, tmax, len(redshifts))
    stretch, color = bg_stretch_color(len(redshifts), **kwargs)

    # Following the design of sncosmo, we return data into a list of dicts
    param_iter = zip(redshifts, peakmjd, stretch, color)
    return [{'z': z, 'x0': 1e-14, 't0': t0, 'x1': x1, 'c': c} for
            z, t0, x1, c in param_iter]


def generate_lc(model, phase_range, model_params, out_dir):
    """Generate 91bg light curves in SDSS band passes

    This is a very simple simulation that assumes observations are performed
    simultaneously at each epoch with a gain of one and zero sky noise.
    Simulations use SDSS bandpasses

    Args:
        model       (Model): The model to be used
        phase_range (tuple): Phase range to simulate light-curves over
        model_params (list): List of model parameters for each light-curve
        out_dir      (Path): The output directory to write light-curves to
    """

    out_dir = Path(out_dir)  # In case out_dir is a string

    for i, params_dict in enumerate(model_params):
        time = 5 * (params_dict['t0'] + np.arange(*phase_range).tolist())
        bands = np.concatenate(
            [np.full(len(time) // 5, 'sdss' + band) for band in 'ugriz']
        )

        table_length = len(time)
        observations = Table(
            {'time': time,
             'band': bands,
             'gain': np.full(table_length, 1.0),
             'skynoise': np.zeros(table_length),
             'zp': np.full(table_length, 27.5),
             'zpsys': np.full(table_length, 'ab')})

        light_curve = sncosmo.realize_lcs(observations, model, [params_dict])[0]
        light_curve.meta.update(params_dict)
        light_curve.write(out_dir / f'sn91bg_{i}.ecsv', overwrite=True)
