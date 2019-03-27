#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Functions for reading in fit results from the analysis pipeline."""

import os

from astropy.table import Table


def get_fit_results(survey, model, params, out_dir='./'):
    """Get lightcurve fits for a given survey and model

    Args:
        survey  (str): The name of the survey
        model (Model): The SNCosmo model used to fit lightcurves
        params  (int): The number of Salt2 params that were fit (4 or 5)
        out_dir (str): The directory specified to the analysis pipeline for
                    outputting results (Default is current working directory)

    Returns:
        A DataFrame of fits in all bands
        A DataFrame of fits in blue bands
        A DataFrame of fits in red bands
    """

    model_name = model.source.name + '_' + model.source.version
    fname = f'{survey}/{model_name}_{params}param_{{}}.ecsv'
    path_pattern = os.path.join(out_dir, fname)

    all_data = None
    all_data_path = path_pattern.format('all')
    if os.path.exists(all_data_path):
        all_data = Table.read(path_pattern.format('all'))
        all_data = all_data.to_pandas()
        all_data.set_index('cid', inplace=True)

    blue_data = None
    blue_data_path = path_pattern.format('blue')
    if os.path.exists(blue_data_path):
        blue_data = Table.read(path_pattern.format('blue'))
        blue_data = blue_data.to_pandas()
        blue_data.set_index('cid', inplace=True)

    red_data = None
    red_data_path = path_pattern.format('red')
    if os.path.exists(red_data_path):
        red_data = Table.read(path_pattern.format('red'))
        red_data = red_data.to_pandas()
        red_data.set_index('cid', inplace=True)

    return all_data, blue_data, red_data
