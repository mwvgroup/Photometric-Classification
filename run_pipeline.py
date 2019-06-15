#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Command line interface for the analysis_pipeline package."""

import argparse

import yaml
from SNData.csp import dr3
from SNData.des import sn3yr
from SNData.sdss import sako18

import analysis_pipeline

out_dir = analysis_pipeline.FIT_DIR
for data in (dr3, sn3yr, sako18):
    data.download_module_data()
    data.register_filters()


def read_yaml(file_path):
    """A yaml file reader compatible with Python 3.6 adn 3.7

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


def run(args):
    """Run light curve fits using command line args"""

    import sncosmo

    from analysis_pipeline import SN91bgSource
    from analysis_pipeline import iter_all_fits

    # Define surveys and models for fitting
    models_dict = dict(
        salt_2_4=sncosmo.Model(
            source=sncosmo.get_source('salt2', version='2.4')),
        salt_2_0=sncosmo.Model(
            source=sncosmo.get_source('salt2', version='2.0')),
        sn_91bg=sncosmo.Model(source=SN91bgSource())
    )
    models = [models_dict[model_name] for model_name in args.models]
    survey = {'csp': dr3, 'des': sn3yr, 'sdss': sako18}[args.survey]

    # Run fitting
    kwargs = read_yaml(args.args_path)[args.survey]
    iter_all_fits(
        out_dir=out_dir,
        module=survey,
        models=models,
        num_params=args.num_params,
        time_out=args.time_out,
        kwargs=kwargs,
        skip_types=args.skip_types)


# Parse command line input
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Fit light-curves for a given survey.')

    parser.add_argument(
        '-s', '--survey',
        type=str,
        required=True,
        help='Survey name (csp, des, sdss)')

    parser.add_argument(
        '-a', '--args_path',
        type=str,
        required=True,
        help='Path of fitting arguments')

    parser.add_argument(
        '-m', '--models',
        type=str,
        nargs='+',
        default=['salt_2_4'],
        help='Models to fit (salt_2_0, salt_2_4, sn_91bg)')

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
        default=None,
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
