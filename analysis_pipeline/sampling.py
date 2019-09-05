#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This module determines and caches prior values specific to a given
light-curve and model. Prior values are determined using nested sampling.

Todo: This module is currently a dump of old logic. It's saved as a reference
for future work only
"""

from copy import deepcopy
from pathlib import Path

import sncosmo
from astropy.table import Column, Table, unique, vstack


PRIORS = dict()  # For lazy loading priors

# !/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Define and provide package wide access to file paths of priors
and fit results.
"""

import os
from pathlib import Path

if 'FIT_RES_DIR' in os.environ:
    OUT_DIR = Path(os.environ['FIT_RES_DIR'])

else:
    OUT_DIR = Path(__file__).resolve().parent / 'results'

FIT_DIR = OUT_DIR / 'fit_results'
PRIOR_DIR = OUT_DIR / 'priors'

FIT_DIR.mkdir(exist_ok=True, parents=True)
PRIOR_DIR.mkdir(exist_ok=True, parents=True)


def get_priors_path(model, survey, manual=False):
    """Construct the path of a priors table for a given model and survey

    Args:
        model (sncosmo.Model): The model to return priors for
        survey       (module): An sndata data access module
        manual         (bool): Return path of manual priors (Default: False)
    """

    survey_abbrev = survey.survey_abbrev.lower()
    release = survey.release.lower()

    # Get path of priors file
    model_name = f'{model.source.name.lower()}_{model.source.version}'
    file_name = f'{survey_abbrev}_{release}_{model_name}.ecsv'
    path = PRIOR_DIR / file_name
    if manual:
        path = path.with_suffix('.man.ecsv')

    return path


def get_fit_result_paths(model, survey, num_params):
    """Create a list of file paths for all, blue, and red band fit results

    Args:
        model (sncosmo.Model): The model to return priors for
        survey       (module): An sndata data access module
        num_params      (int): The number of params being fit

    Returns:
        Path objects for all three tables
    """

    survey_abbrev = survey.survey_abbrev.lower()
    release = survey.release.lower()

    paths = []
    for band in ('all', 'blue', 'red'):
        file_name = (f'{survey_abbrev}_{release}_'
                     f'{model.source.name.lower()}_{model.source.version}_'
                     f'{num_params}_{band}.ecsv')

        paths.append(FIT_DIR / file_name)

    return paths


def create_priors_table(model):
    """Return an empty astropy table for storing fitting priors

    Columns of returned table:
         obj_id, *<model.param_names>, *<model.param_names>_min,
         *<model.param_names>_max

    Args:
        model (Model): An SNCosmo model

    Returns:
        An astropy Table
    """

    col_names = ['obj_id']
    dtype = ['U100']
    for param in model.param_names:
        col_names.extend((param, param + '_min', param + '_max'))
        dtype.extend((float, float, float))

    col_names.append('message')
    dtype.append('U100')
    return Table(names=col_names, dtype=dtype)


def read_priors_table(file_path):
    """Read fitting priors from file

    Read a table of automatically generated priors from file. If manual prior
    values have been set, overwrite the automatic values with manual values in
    the returned file. Assumes table has the same column names as created by
    ``create_priors_table``.

    Args:
        file_path (Path): Path of the table with prior values

    Returns:
        An astropy Table
    """

    priors = Table.read(file_path)

    # Typecast to avoid truncating values when appending to table
    priors['obj_id'] = Column(priors['obj_id'], dtype='U100')
    priors['message'] = Column(priors['message'], dtype='U100')

    # Prefer manually specified priors over automatic ones
    manual_priors_path = file_path.with_suffix('.man.ecsv')
    if manual_priors_path.exists():
        manual_priors = Table.read(manual_priors_path)
        manual_priors['obj_id'] = Column(manual_priors['obj_id'], dtype='U100')
        priors = vstack([manual_priors, priors])
        priors = unique(priors, keys='obj_id', keep='first')

    return priors


def get_cached_priors_table(model, survey):
    """Cache and return tabulated priors for a given model and survey

    Results are read from file and cached in memory.

    Args:
        model   (Model): An SNCosmo model
        survey (module): An sndata data access module

    Returns:
        An astropy Table
    """

    if model in PRIORS:
        return PRIORS[model]

    file_path = get_priors_path(model, survey)
    if file_path.exists():
        priors_table = read_priors_table(file_path)

    else:
        priors_table = create_priors_table(model)

    return PRIORS.setdefault(model, priors_table)


# noinspection PyIncorrectDocstring
def nest_lc(data, model, vparam_names, **kwargs):
    """Wrapper for sncosmo.nest_lc

    Use nested sampling to determine initial model parameters

    Args:
        data             (Table): Table of light curve data for SNCosmo
        model            (Model): SNCosmo model
        vparam_names (list[str]): List of parameters to vary in the model
        bounds            (dict): Boundaries on fit parameters
        verbose           (bool): Whether to display progress (Default: True)
        maxiter            (int): Maximum sampling iterations (Default: 10000)
        maxcall            (int): Maximum function calls (Default: 20000)
        method             (str): Nested sampling method (Default: 'multi')

    Returns:
        A sampled SNCosmo Model
        A dictionary of parameter bounds used to generate the sampled model
    """

    # Protect against argument mutation
    kwargs = deepcopy(kwargs)

    # Define our own default values and get model with nested values
    kwargs['verbose'] = kwargs.get('verbose', True)
    kwargs['maxiter'] = kwargs.get('maxiter', 10000)
    kwargs['maxcall'] = kwargs.get('maxcall', 20000)
    kwargs['warn'] = kwargs.get('warn', False)

    _, sampled_model = sncosmo.nest_lc(data, model, vparam_names, **kwargs)

    return sampled_model, kwargs['bounds']


def calc_new_priors(data, model, vparam_names, survey, time_out=120, **kwargs):
    """Use nested sampling update cached prior values

    If not specified, bounds on t0 assume the first observation is less than
    15 days after peak and before the last data point. If sampling fails or
    times out, return the original model values.

    Args:
        data             (Table): Table of light curve data for SNCosmo
        model            (Model): SNCosmo model
        vparam_names (list[str]): List of parameters to vary in the model
        survey          (module): An sndata data access module
        time_out           (int): Seconds before timing out (Default: 120)
        Any other arguments for analysis_pipeline.lc_fitting.nest_lc

    Returns:
        The sampled sncosmo model
        The bounds used to sample the model
    """

    kwargs = deepcopy(kwargs)

    # Assume first observation < 15 days after peak and before last data point
    t0_start = min(data['time']) - 15
    t0_end = min(max(data['time']), min(data['time']) + 25)
    kwargs['bounds'].setdefault('t0', (t0_start, t0_end))

    try:
        with timeout(seconds=time_out):
            sampled_model, bounds = nest_lc(data, model, vparam_names,
                                            **kwargs)
            msg = 'Success'

    except (ValueError, TimeoutError, RuntimeError) as e:
        sampled_model = model
        bounds = kwargs['bounds']
        msg = f'Fail {e} - defaulting to original model'

    # Format prior data for appending to table
    prior_data = [data.meta['obj_id']]
    for param_name, param_val in zip(sampled_model.param_names,
                                     sampled_model.parameters):
        prior_data.append(param_val)
        prior_data.extend(kwargs['bounds'][param_name])

    prior_data.append(msg)

    # Update cached values
    PRIORS[model].add_row(prior_data)
    PRIORS[model].write(get_priors_path(model, survey), overwrite=True)

    return sampled_model, bounds


def get_sampled_model(data, model, vparam_names, survey, time_out=120,
                      **kwargs):
    """Set initial model params using cached values

    If cached values are not available, determine them using nested sampling
    and save them to file.

    Args:
        data             (Table): Table of light curve data for SNCosmo
        model            (Model): SNCosmo model
        vparam_names (list[str]): List of parameters to vary in the model
        survey          (module): An sndata data access module
        time_out           (int): Seconds before sampling times out (Default: 120)
        Any other arguments for analysis_pipeline.lc_fitting.nest_lc
    """

    priors_table = get_cached_priors_table(model, survey)
    cached_prior = priors_table[priors_table['obj_id'] == data.meta['obj_id']]

    if cached_prior:
        sampled_model = deepcopy(model)
        sampled_model.update({k: cached_prior[k] for k in vparam_names})
        bounds = {k: (cached_prior[f'{k}_min'], cached_prior[f'{k}_max']) for k
                  in vparam_names}

    else:
        sampled_model, bounds = calc_new_priors(
            data, model, vparam_names, survey, time_out, **kwargs)

    return sampled_model, bounds


def get_priors(survey, model):
    """Get light-curve priors for a given survey and model

    Args:
        survey  (module): An SNData submodule for a particular data release
        model    (Model): The SNCosmo model used to fit light-curves

    Returns:
        An astropy table
    """

    path = _paths.get_priors_path(model, survey)
    data = Table.read(path).to_pandas()
    data.set_index('obj_id', inplace=True)
    return data


def save_manual_priors(obj_id, survey, model, priors_dict, message='-'):
    """Save priors to the analysis pipeline's internal file structure

    Args:
        obj_id       (str): The ID of the object to save priors for
        survey    (module): An sndata data access module
        model      (Model): The SNCosmo model of the priors
        priors_dict (dict): Dictionary of prior values
        message (str): Message to include with new priors (Default: '-')
    """

    auto_priors_path = _paths.get_priors_path(model, survey, manual=True)
    manual_priors_path = _paths.get_priors_path(model, survey, manual=True)

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
