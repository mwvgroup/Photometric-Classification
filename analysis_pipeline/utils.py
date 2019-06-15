#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""A collection of general utilites"""

import os as _os
import signal
from pathlib import Path

from astropy.table import Column, Table, unique

FIT_DIR = Path(__file__).resolve().parent / 'fit_results'
PRIOR_DIR = Path(__file__).resolve().parent / 'priors'
PRIOR_DIR.mkdir(exist_ok=True)


class timeout:
    """A timeout context manager"""

    def __init__(self, seconds=1, error_message='Timeout'):
        """A timeout context manager
        Args:
            seconds       (int): The number of seconds until timeout
            error_message (str): The TimeOutError message on timeout
        """

        self.seconds = seconds
        self.error_message = error_message

    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)

    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)

    def __exit__(self, type, value, traceback):
        signal.alarm(0)


def get_fit_results(survey, model, params):
    """Get light-curve fits for a given survey, model, and number of parameters

    Args:
        survey (module): An SNData submodule for a particular data release
        model   (Model): The SNCosmo model used to fit light-curves
        params    (int): The number of Salt2 params that were fit (4 or 5)

    Returns:
        A DataFrame of fits in all bands or None
        A DataFrame of fits in blue bands or None
        A DataFrame of fits in red bands or None
    """

    survey = survey.survey_abbrev.lower()
    model_name = model.source.name + '_' + model.source.version
    fname = f'{survey}_{params}_{model_name}_{{}}.ecsv'

    all_data = None
    all_data_path = FIT_DIR / fname.format('all')
    if all_data_path.exists():
        all_data = Table.read(all_data_path)
        all_data = all_data.to_pandas()
        all_data.set_index('obj_id', inplace=True)

    blue_data = None
    blue_data_path = FIT_DIR / fname.format('blue')
    if _os.path.exists(blue_data_path):
        blue_data = Table.read(blue_data_path)
        blue_data = blue_data.to_pandas()
        blue_data.set_index('obj_id', inplace=True)

    red_data = None
    red_data_path = FIT_DIR / fname.format('red')
    if _os.path.exists(red_data_path):
        red_data = Table.read(red_data_path)
        red_data = red_data.to_pandas()
        red_data.set_index('obj_id', inplace=True)

    return all_data, blue_data, red_data


def _get_priors_paths(survey, model):
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
    auto_priors_path = PRIOR_DIR / auto_priors_name
    manual_priors_path = auto_priors_path.with_suffix('.man.ecsv')
    return auto_priors_path, manual_priors_path


def save_priors(obj_id, module, model, priors_dict, message='-'):
    """Save priors to the analysis pipeline's internal file structure

    Args:
        obj_id       (str): The ID of the object to save priors for
        module       (str): An SNData module
        model      (Model): The SNCosmo model of the priors
        priors_dict (dict): Dictionary of prior values
        message (str): Message to include with new priors (Default: '-')
    """

    auto_priors_path, manual_priors_path = \
        _get_priors_paths(module.survey_abbrev.lower(), model)

    try:
        existing_data = Table.read(manual_priors_path)
        existing_data['message'] = Column(existing_data['message'],
                                          dtype='U100')

    except FileNotFoundError:
        auto_data = Table.read(auto_priors_path)
        names = auto_data.colnames
        dtype = [float for _ in names]
        dtype[0] = dtype[-1] = 'U1000'
        existing_data = Table(names=names, dtype=dtype)

    priors_dict['obj_id'] = obj_id
    priors_dict['message'] = message
    existing_data.add_row(priors_dict)
    new_data = unique(existing_data, keep='last', keys=['obj_id'])
    new_data.write(manual_priors_path, overwrite=True)
