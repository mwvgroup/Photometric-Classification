#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

import os as _os
from pathlib import Path as _Path

from astropy.table import Table as _Table

from ._sn91bg_model._model import SN91bgSource

_FIT_DIR = _Path(__file__).resolve().parent.parent / 'fit_results'


def get_fit_results(survey, model, params):
    """Get light-curve fits for a given survey and model

    Args:
        survey  (str): The name of the survey
        model (Model): The SNCosmo model used to fit light-curves
        params  (int): The number of Salt2 params that were fit (4 or 5)

    Returns:
        A DataFrame of fits in all bands or None
        A DataFrame of fits in blue bands or None
        A DataFrame of fits in red bands or None
    """

    survey = survey.lower()
    model_name = model.source.name + '_' + model.source.version
    fname = f'{survey}/{survey}_{params}_{model_name}_{{}}.ecsv'

    all_data = None
    all_data_path = _FIT_DIR / fname.format('all')
    if all_data_path.exists():
        all_data = _Table.read(all_data_path)
        all_data = all_data.to_pandas()
        all_data.set_index('obj_id', inplace=True)

    blue_data = None
    blue_data_path = _FIT_DIR / fname.format('blue')
    if _os.path.exists(blue_data_path):
        blue_data = _Table.read(blue_data_path)
        blue_data = blue_data.to_pandas()
        blue_data.set_index('obj_id', inplace=True)

    red_data = None
    red_data_path = _FIT_DIR / fname.format('red')
    if _os.path.exists(red_data_path):
        red_data = _Table.read(red_data_path)
        red_data = red_data.to_pandas()
        red_data.set_index('obj_id', inplace=True)

    return all_data, blue_data, red_data
