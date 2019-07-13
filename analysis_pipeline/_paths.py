#!/usr/bin/env python3.7
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
