#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Command line interface for the ``phot_class`` package."""

import argparse
import warnings
from pathlib import Path

import numpy as np
import sndata
import yaml

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
    # 'SLSN', 'SNIb', 'SNIc', 'pSNIbc', 'zSNIbc', 'Unknown', 'AGN', 'Variable'
    filter_func = utils.classification_filter_factory(
        ['SNIa', 'SNIa?', 'pSNIa', 'zSNIa']
    )

    data_iter = data_module.iter_data(
        verbose={'desc': f'{survey} - {release} Objects'},
        filter_func=filter_func)

    for data in data_iter:
        data = data[data['band'] != 'csp_dr3_Ydw']
        data = data[data['band'] != 'csp_dr3_Y']
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


def get_spec_data_iter(data_module):
    """Iterate over spectroscopic data tables while removing galaxy observations

    This function is a wrapper around ``data_module.iter_data``

    Args:
        data_module (module): An sndata module

    Yields:
        Astropy tables
    """

    survey = data_module.survey_abbrev
    release = data_module.release
    if data_module.data_type != 'spectroscopic':
        raise RuntimeError(
            f'{survey} - {release} is not a spectroscopic data release')

    for table in data_module.iter_data(verbose={'desc': 'Objects'}):
        # Skip sdss galaxy spectra
        if survey.lower() == 'sdss':
            table = table[np.isin(table['type'], ['Ia', 'Ia-pec', 'Ia?'])]

        if not table:
            continue

        for spectrum in table.group_by('date').groups:
            if spectrum.meta['ra'] is None:
                continue

            yield spectrum


def run_spectroscopic_classification(cli_args):
    """Run spectroscopic classification of SNe

    Args:
        cli_args (argparse.Namespace): Command line arguments
    """

    # Create output file path
    out_dir = Path(cli_args.out_dir).resolve() / 'spec_class'
    out_dir.mkdir(exist_ok=True, parents=True)

    rv_str = str(cli_args.rv).replace('.', '_')
    file_name = f'{cli_args.survey}_{cli_args.release}_{rv_str}_{cli_args.nstep}.ecsv'

    data_module = getattr(getattr(sndata, cli_args.survey), cli_args.release)
    data_module.download_module_data()

    out_table = spectra.tabulate_spectral_properties(
        get_spec_data_iter(data_module),
        rv=cli_args.rv,
        nstep=cli_args.nstep,
        plot=cli_args.plot)

    out_table.write(out_dir / file_name, overwrite=True)


def create_cli_parser():
    """Return a command line argument parser"""

    parser = argparse.ArgumentParser(
        description='Arguments for the command line interface are as follows:')
    subparsers = parser.add_subparsers()

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

    photometric_parser = subparsers.add_parser('photometric')
    photometric_parser.set_defaults(
        func=run_photometric_classification,
        help='Classify targets photometrically'
    )

    photometric_parser.add_argument(
        '-f', '--fit_func',
        type=str,
        default='simple_fit',
        help='The name of the fitting routine to use (simple_fit, nest_fit, mcmc_fit).'
    )

    photometric_parser.add_argument(
        '-m', '--method',
        type=str,
        default='band',
        help="Whether to fit bands independently ('band') or as red and blue sets ('collective')"
    )

    photometric_parser.add_argument(
        '-c', '--config',
        type=str,
        required=False,
        help='Path of the yaml config file.'
    )

    photometric_parser.add_argument(
        '-o', '--out_dir',
        type=str,
        required=True,
        help='Directory to write output files to.'
    )

    spectroscopic_parser = subparsers.add_parser('spectroscopic')
    spectroscopic_parser.set_defaults(
        func=run_spectroscopic_classification,
        help='Classify targets spectroscopically'
    )

    spectroscopic_parser.add_argument(
        '-r', '--rv',
        type=float,
        default=3.1,
        help='Rv value to use for extinction correction'
    )

    spectroscopic_parser.add_argument(
        '-n', '--nstep',
        type=int,
        default=5,
        help='Number of steps used in resampling'
    )

    spectroscopic_parser.add_argument(
        '-b', '--bin_size',
        type=int,
        default=5,
        help='Size of bins in angstroms'
    )

    spectroscopic_parser.add_argument(
        '-m', '--method',
        type=str,
        default='avg',
        help='Either "avg" or "sum" each bin'
    )

    spectroscopic_parser.add_argument(
        '-o', '--out_dir',
        type=str,
        required=True,
        help='Directory to write output files to.'
    )

    spectroscopic_parser.add_argument(
        '--plot',
        help='Display live plots of fitting results for the velocity calculation.',
        action='store_true')

    return parser


# Parse command line input
if __name__ == '__main__':
    from phot_class import classification
    from phot_class import fit_func_wraps
    from phot_class import models
    from phot_class import utils
    from phot_class import spectra

    models.register_sources()

    parser = create_cli_parser()
    cli_args = parser.parse_args()
    cli_args.func(cli_args)
