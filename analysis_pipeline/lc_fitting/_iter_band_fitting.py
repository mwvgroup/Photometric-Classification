#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module handles the iterative fitting of light-curves for multiple
models, surveys, and band pass collections.
"""

from pathlib import Path

import numpy as np
import sncosmo
from astropy.table import Table

from ._fit_funcs import fit_lc, nest_lc
from .._sn91bg_model import SN91bgSource


def _split_bands(bands, lambda_eff):
    """Split band-passes into collections of blue and red bands

    Blue bands have an effective wavelength < 5500 Ang. Red bands have an
    effective wavelength >= 5500 Ang.

    Args:
        bands        (array[str]): Name of band-passes
        lambda_eff (array[float]): Effective wavelength of band-passes
    """

    is_blue = np.array(lambda_eff) < 5500
    b_array = np.array(bands)
    return b_array[is_blue], b_array[~is_blue]


def _split_data(data_table, band_names, lambda_eff):
    """Split a data table into blue and red data (by restframe)

    Split data by keeping filters that are redward or blueward of 5500 Ang.

    Args:
        data_table (Table): An SNCosmo input table with column 'band'
        band_names  (iter): List of all bands available in the survey
        lambda_eff  (iter): The effective wavelength of each band in band_names

    Returns:
        A new input table for SNCosmo only containing select rest frame bands
    """

    # Type cast to allow numpy indexing
    band_names = np.array(band_names)
    lambda_eff = np.array(lambda_eff)

    @np.vectorize
    def lambda_for_band(band):
        return lambda_eff[band_names == band]

    # Rest frame effective wavelengths for each observation
    z = data_table.meta['redshift']
    observed_lambda = lambda_for_band(data_table['band'])
    rest_frame_lambda = observed_lambda / (1 + z)

    # Get the name of the observer frame band with the smallest distance
    # to each rest frame lambda
    delta_lambda = np.array([
        np.abs(rest_frame_lambda - l_eff) for l_eff in lambda_eff])

    min_indx = np.argmin(delta_lambda, axis=0)
    rest_frame_filters = np.array(band_names)[min_indx]

    # Keep only the specified filters that are within 700 Angstroms of the
    # rest frame effective wavelength
    is_ok_diff = delta_lambda[min_indx, np.arange(delta_lambda.shape[1])] < 700

    # Split into blue and red band passes
    out_list = []
    for bands in _split_bands(band_names, lambda_eff):
        is_in_bands = np.isin(rest_frame_filters, bands)
        indices = np.logical_and(is_in_bands, is_ok_diff)
        out_list.append(data_table[indices])

    return out_list


def _create_empty_summary_table():
    """Returns a table with columns:

         cid, num_points, z, t0, x0, x1, z_err, t0_err, x0_err,
         x1_err, c_err, chi, dof, tmin, tmax, pre_max, post_max, message
    """

    names = ['cid', 'num_points']
    dtype = ['U20', int]

    param_names = ('z', 't0', 'x0', 'x1', 'c')
    names.extend(param_names)
    names.extend((p + '_err' for p in param_names))
    dtype.extend((float for _ in range(len(names) - 2)))

    names.extend(
        ('chi', 'dof', 'b_max',
         'delta_15', 'tmin', 'tmax',
         'pre_max', 'post_max', 'message'))

    dtype.extend([float, float, float, float, float, float, int, int, 'U1000'])
    return Table(names=names, dtype=dtype)


def _create_table_paths(out_dir, survey):
    """Create a list of file paths

    Args:
        out_dir (str): Out directory of file paths
        survey  (str): Name of a survey

    Returns:
        [">out_dir>/<survey>_<bands>" for bands in ('all', 'blue', 'red')]
    """
    out_dir = Path(out_dir)
    return [out_dir / f'{survey}_{bands}' for bands in ('all', 'blue', 'red')]


def iter_all_fits(out_dir, module, nkwargs=dict(), fkwargs=dict()):
    """Iteratively fit data for a given survey

    Args:
        out_dir   (str): Directory to write fit results to
        module (module): A data access module
        nkwargs  (dict): A dictionary of kwargs for sncosmo.nest_lc
        fkwargs  (dict): A dictionary of kwargs for sncosmo.fit_lc
    """

    salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))
    salt_2_0 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.0'))
    sn_91bg = sncosmo.Model(source=SN91bgSource())
    models_dict = dict(salt_2_4=salt_2_4, salt_2_0=salt_2_0, sn_91bg=sn_91bg)

    for model_name, model in models_dict.items():

        # Iteratively perform a 4 and 5 parameter fit
        param_names = ['z', 't0', 'x0', 'x1', 'c']
        for num_params in (4, 5):
            params_to_fit = param_names[5 - num_params:]

            # Create separate tables for each band and populate with fit results
            out_tables = [_create_empty_summary_table() for _ in range(3)]
            out_paths = _create_table_paths(out_dir, module.survey_name)
            for data in module.iter_sncosmo_inputs():
                sampled_model = nest_lc(data, model, params_to_fit, **nkwargs)

                band_data = [data]
                band_data.extend(_split_data(
                    data, module.band_names, module.lambda_effective))

                for table, path, input_table in zip(out_tables, out_paths,
                                                    band_data):
                    fkwargs['guess_amplitude'] = False
                    fkwargs['guess_t0'] = False
                    fkwargs['guess_z'] = False

                    fit_results = fit_lc(
                        input_table, sampled_model, params_to_fit, **fkwargs)

                    table.add_row(fit_results)
                    table.write(path)
