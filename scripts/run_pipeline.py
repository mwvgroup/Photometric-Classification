#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Command line interface for the ``phot_class`` package."""

import argparse
import math
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import sndata
import yaml
from sndata.sdss import sako18spec

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from phot_class import classification
from phot_class import fit_func_wraps
from phot_class import models
from phot_class import utils

warnings.simplefilter('ignore')


def load_yaml(path):
    """Return data from a yaml file using the full loader

    Args:
        path (str): File path to read from

    Returns:
        File contents as a dictionary
    """

    with open(path) as infile:
        return yaml.load(infile, Loader=yaml.FullLoader)


def get_phot_data_iter(data_module):
    """Iterate over photometric data tables while removing Y band observations

    This function is a wrapper around ``data_module.iter_data``

    Args:
        data_module (module): An sndata module

    Yields:
        Astropy tables
    """

    survey = data_module.survey_abbrev
    release = data_module.release
    if data_module.data_type != 'photometric':
        raise RuntimeError(
            f'{survey} - {release} is not a photometric data release')

    # Other classifications:
    # 'SNIa', 'SNIa?', 'pSNIa', 'zSNIa'
    # 'SLSN', 'SNIb', 'SNIc', 'pSNIbc', 'zSNIbc', 'Unknown', 'AGN', 'Variable'
    filter_func = utils.classification_filter_factory(
        ['AGN', 'Variable']
    )

    data_iter = data_module.iter_data(
        verbose={'desc': f'{survey} - {release} Objects'},
        filter_func=filter_func)

    for data in data_iter:
        data = data[data['band'] != 'csp_dr3_Ydw']
        data = data[data['band'] != 'csp_dr3_Y']
        data = data[data['band'] != 'csp_dr3_Jrc2']
        data = data[data['band'] != 'csp_dr3_Jdw']
        data = data[data['band'] != 'csp_dr3_J']
        data = data[data['band'] != 'csp_dr3_Hdw']
        data = data[data['band'] != 'csp_dr3_H']
        yield data


def run_photometric_classification(cli_args):
    """Run photometric classification of SNe

    Args:
        cli_args (argparse.Namespace): Command line arguments
    """

    # Create output file paths
    out_dir = Path(cli_args.out_dir).resolve()
    out_dir.mkdir(exist_ok=True, parents=True)

    file_prefix = f'{cli_args.survey}_{cli_args.release}_{cli_args.fit_func}_'
    fit_path = out_dir / (file_prefix + 'fits.ecsv')
    classification_path = out_dir / (file_prefix + 'class.ecsv')

    # Get the sndata module - download and register data for fitting
    data_module = getattr(getattr(sndata, cli_args.survey), cli_args.release)
    data_module.download_module_data()
    data_module.register_filters()

    # specify arguments for fitting.tabulate_fit_results
    data_iter = get_phot_data_iter(data_module)
    band_names = data_module.band_names
    lambda_eff = data_module.lambda_effective
    fit_func = getattr(fit_func_wraps, cli_args.fit_func)

    # Read in priors and fitting arguments from file
    config = load_yaml(cli_args.config) if cli_args.config else None

    # Run fits
    fit_results = classification.tabulate_fit_results(
        data_iter=data_iter,
        band_names=band_names,
        lambda_eff=lambda_eff,
        fit_func=fit_func,
        fitting_method=cli_args.method,
        config=config,
        out_path=fit_path
    )

    classification.classify_targets(fit_results, out_path=classification_path)


@np.vectorize
def calc_julian_date(date):
    """
    Convert a datetime object into julian float.

    Args:
        date (str): The date to convert in %Y-%m-%d format

    Returns:
        The Julian date as a float
    """

    date = datetime.strptime(date, '%Y-%m-%d')
    julian_datetime = (
            367 * date.year -
            int((7 * (date.year + int((date.month + 9) / 12.0))) / 4.0) +
            int((275 * date.month) / 9.0) + date.day +
            1721013.5 +
            (date.hour + date.minute / 60.0 + date.second / math.pow(60, 2)) / 24.0 -
            0.5 * math.copysign(1, 100 * date.year + date.month - 190002.5) + 0.5
    )

    return julian_datetime


def get_sdss_phase(obj_id, date_str):
    peak = sako18spec.load_table('master').to_pandas(index='CID')['MJDatPeakrmag']
    peak += 2400000.5

    if obj_id in peak.index:
        return calc_julian_date(date_str) - peak.loc[obj_id]

    else:
        return None


def create_cli_parser():
    """Return a command line argument parser"""

    parser = argparse.ArgumentParser(
        description='Arguments for the command line interface are as follows:')
    parser.set_defaults(
        func=run_photometric_classification,
        help='Classify targets photometrically'
    )

    parser.add_argument(
        '-s', '--survey',
        type=str,
        required=True,
        help='The name of the survey to analyze. This should be the name of a survey in the sndata package (e.g. csp).'
    )

    parser.add_argument(
        '-r', '--release',
        type=str,
        required=True,
        help='The name of the survey\'s data release. This should also match the sndata package (e.g. dr3).'
    )

    parser.add_argument(
        '-f', '--fit_func',
        type=str,
        default='simple_fit',
        help='The name of the fitting routine to use (simple_fit, nest_fit, mcmc_fit).'
    )

    parser.add_argument(
        '-m', '--method',
        type=str,
        default='band',
        help="Whether to fit bands independently ('band') or as red and blue sets ('collective')"
    )

    parser.add_argument(
        '-c', '--config',
        type=str,
        required=False,
        help='Path of the yaml config file.'
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
    models.register_sources()

    parser = create_cli_parser()
    cli_args = parser.parse_args()
    cli_args.func(cli_args)
