#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Command line interface for the analysis_pipeline package."""

import argparse
import warnings
from pathlib import Path

import sndata

from analysis_pipeline import classification
from analysis_pipeline import fit_funcs
from analysis_pipeline import models

warnings.simplefilter('ignore')
models.register_sources()


def run(cli_args):
    """Run light curve fits using command line args

    Args:
        cli_args (argparse.Namespace): Command line arguments
    """

    # Create output file paths
    out_dir = Path(out_path=cli_args.out_dir)
    out_dir.mkdir(exist_ok=True, parents=True)
    fit_path = out_dir/ 'fits.ecsv'
    classification_path = out_dir / 'classification.ecsv'

    # Download and register data for fitting
    data_module = getattr(getattr(sndata, cli_args.survey), cli_args.release)
    data_module.download_module_data()
    data_module.register_filters()

    # specify arguments for classification.tabulate_fit_results
    data_iter = data_module.iter_data(format_sncosmo=True, verbose=True)
    band_names = data_module.band_names
    lambda_eff = data_module.lambda_effective
    fit_func = getattr(fit_funcs, cli_args.fit_func)
    vparams = cli_args.vparams
    timeout_sec = cli_args.timeout

    # Todo: this should be specified externally from the CLI somehow
    kwargs_bg = {'bounds': {
        'z': [0.01, 0.8],
        'x0': [0, 0.02],
        'x1': [0.65, 1.25],
        'c': [0, 1]}
    }

    fit_results = classification.tabulate_fit_results(
        data_iter, band_names, lambda_eff, fit_func, vparams,
        timeout_sec=timeout_sec, kwargs_bg=kwargs_bg, out_path=fit_path)

    classification.classify_targets(fit_results, out_path=classification_path)


def create_cli_parser():
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
        '-f', '--fit_func',
        type=str,
        default='simple_fit',
        help='Which fitting function to use')

    parser.add_argument(
        '-v', '--vparams',
        type=str,
        nargs='+',
        required=True,
        help='What parameters to vary with the fit')

    parser.add_argument(
        '-t', '--timeout',
        type=int,
        default=120,
        help='Seconds before fitting times out.'
    )

    parser.add_argument(
        '-o', '--out_dir',
        type=str,
        required=True,
        help='Directory to write output files to.'
    )

    return parser


# Parse command line input
if __name__ == '__main__':
    parser = create_cli_parser()
    cli_args = parser.parse_args()
    run(cli_args)
