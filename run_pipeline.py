#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Command line interface for the analysis_pipeline package."""

import argparse

import sncosmo
import sndata
import yaml

import analysis_pipeline

analysis_pipeline.models.register_sources()
out_dir = analysis_pipeline.FIT_DIR


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

    models_list = list()
    for name, version in zip(cli_args.models, cli_args.versions):
        source = sncosmo.get_source(name, version=version)
        model = sncosmo.Model(source=source)
        models_list.append(model)

    return models_list


def run(cli_args):
    """Run light curve fits using command line args

    Args:
        cli_args (argparse.Namespace): Command line arguments
    """

    models = get_models(cli_args)
    data_module = getattr(getattr(sndata, cli_args.survey), cli_args.release)
    data_module.download_module_data()
    data_module.register_filters()

    # Run fitting
    kwargs = read_yaml(cli_args.args_path)[cli_args.survey]
    analysis_pipeline.lc_fitting.iter_all_fits(
        out_dir=cli_args.out_dir,
        module=data_module,
        models=models,
        num_params=cli_args.num_params,
        time_out=cli_args.time_out,
        kwargs=kwargs,
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
        '-n', '--num_params',
        type=int,
        nargs='+',
        default=[4, 5],
        help='Number of params to fit (4, 5)')

    parser.add_argument(
        '-k', '--skip_types',
        type=str,
        nargs='+',
        default=['AGN', 'SLSN', 'SNII', 'Variable'],
        help='Object classifications to skip. Only supported for SDSS.'
    )

    parser.add_argument(
        '-o', '--out_dir',
        type=str,
        default=out_dir,
        help='Output directory for fit results.'
    )

    parser.add_argument(
        '-t', '--time_out',
        type=int,
        default=90,
        help='Seconds before nested sampling times out.'
    )

    cli_args = parser.parse_args()
    if cli_args.survey not in ('csp', 'des', 'sdss'):
        raise ValueError(f"Survey name '{cli_args.survey}' not recognized")

    run(cli_args)
