#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Functions for reading in fit results from the analysis pipeline."""

import os

import pandas as _pd


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

    index_col = 0
    model_name = model.source.name + '_' + model.source.version
    fname = f'{survey}/{model_name}_{params}param_{{}}.csv'
    path_pattern = os.path.join(out_dir, fname)

    all_data = _pd.read_csv(path_pattern.format('all'), index_col=index_col)

    blue_data = _pd.read_csv(path_pattern.format('blue'),
                             index_col=index_col)

    red_data = _pd.read_csv(path_pattern.format('red'),
                            index_col=index_col)

    return all_data, blue_data, red_data
