#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Command line interface for the analysis_pipeline package."""

import argparse

import sncosmo
import sndata
import yaml

from analysis_pipeline import models
from analysis_pipeline.lc_fitting import run_iter_fitting

models.register_sources()


def read_yaml(file_path):
    """A yaml file reader compatible with Python 3.6 and 3.7

    Args:
        file_path (str): The yaml file path to read

    Returns:
        A dict of the file's contents
    """

    with open(file_path) as ofile:
        try:
            return yaml.load(ofile, Loader=yaml.FullLoader)

        except AttributeError:
            return yaml.load(ofile)


def get_models(cli_args):
    """Return a list of SNCosmo models specified by command line arguments

    Args:
        cli_args (argparse.Namespace): Command line arguments

    Returns:
        A list of SNCosmo models
    """

    kwargs = read_yaml(cli_args.args_path)[cli_args.survey]

    models_list = []
    kwargs_list = []
    for name, version in zip(cli_args.models, cli_args.versions):
        source = sncosmo.get_source(name, version=version)
        model = sncosmo.Model(source=source)
        models_list.append(model)
        kwargs_list.append(kwargs[f'{name}_{version}'])

    return models_list, kwargs_list


def run(cli_args):
    """Run light curve fits using command line args

    Args:
        cli_args (argparse.Namespace): Command line arguments
    """

    # Download and register data for fitting
    data_module = getattr(getattr(sndata, cli_args.survey), cli_args.release)
    data_module.download_module_data()
    data_module.register_filters()

    # Get list of specified models and fit them all
    model_lists, kwarg_list = get_models(cli_args)
    run_iter_fitting(
        survey=data_module,
        model_list=model_lists,
        kwarg_list=kwarg_list,
        fitz_list=cli_args.fit_z,
        time_out=cli_args.time_out,
        skip_types=cli_args.skip_types)


# Parse command line input
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Fit light-curves for a given survey.')

    parser.add_argument(
        '-s', '--survey',
        type=str,
        required=True,
        help='Survey name (e.g. csp)')

    parser.add_argument(
        '-r', '--release',
        type=str,
        required=True,
        help='Release name (e.g. dr3)')

    parser.add_argument(
        '-a', '--args_path',
        type=str,
        required=True,
        help='Path of fitting arguments')

    parser.add_argument(
        '-m', '--models',
        type=str,
        nargs='+',
        default=['salt2'],
        help='Models to fit (salt2, sn91bg)')

    parser.add_argument(
        '-v', '--versions',
        type=str,
        nargs='+',
        default=['2.4'],
        help='Version of each model to fit')

    parser.add_argument(
        '-z', '--fit_z',
        nargs='+',
        type=int,
        default=False,
        help='Whether to fit for redshift')

    parser.add_argument(
        '-k', '--skip_types',
        type=str,
        nargs='+',
        default=['AGN', 'SLSN', 'SNII', 'Variable'],
        help='Object classifications to skip. Only supported for SDSS.'
    )

    parser.add_argument(
        '-t', '--time_out',
        type=int,
        default=120,
        help='Seconds before nested sampling times out.'
    )

    cli_args = parser.parse_args()
    if not (len(cli_args.models) ==
            len(cli_args.versions) ==
            len(cli_args.fit_z)):

        raise ValueError(
            'Number of models, version, and redshift flags must match.')

    run(cli_args)
