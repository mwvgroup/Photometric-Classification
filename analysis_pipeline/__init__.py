#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

import os as _os
from pathlib import Path as _Path

from astropy.table import Table as _Table

from ._lc_fitting import fit_lc, split_data
from ._lc_fitting import get_sampled_model
from ._lc_fitting import iter_all_fits
from ._lc_fitting import nest_lc
from ._sn91bg_model._model import SN91bgSource

_FIT_DIR = _Path(__file__).resolve().parent / 'fit_results'
_PRIOR_DIR = _Path(__file__).resolve().parent / 'priors'


def get_fit_results(survey, model, params):
    """Get light-curve fits for a given survey, model, and number of parameters

    Args:
        survey  (str): The name of the survey
        model (Model): The SNCosmo model used to fit light-curves
        params  (int): The number of Salt2 params that were fit (4 or 5)

    Returns:
        A DataFrame of fits in all bands or None
        A DataFrame of fits in blue bands or None
        A DataFrame of fits in red bands or None
    """

    survey = survey.survey_abbrev.lower()
    model_name = model.source.name + '_' + model.source.version
    fname = f'{survey}_{params}_{model_name}_{{}}.ecsv'

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


def get_priors_paths(survey, model):
    """Return the path to the light-curve priors for a given survey and model

    Args:
        survey  (str): The name of the survey
        model (Model): The SNCosmo model used to fit light-curves

    Returns:
        A Path object for the auto generated priors
        A Path object for the manually set priors
    """

    model_name = f'{model.source.name}_{model.source.version}'
    auto_priors_name = f'{survey.lower()}_{model_name}.ecsv'
    auto_priors_path = _PRIOR_DIR / auto_priors_name
    manual_priors_path = auto_priors_path.with_suffix('.man.ecsv')
    return auto_priors_path, manual_priors_path
