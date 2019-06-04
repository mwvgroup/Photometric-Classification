#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Command line interface for the analysis_pipeline package."""

import argparse


def run(args):
    """Run light curve fits using command line args"""

    import sncosmo

    from analysis_pipeline import SN91bgSource
    from analysis_pipeline.lc_fitting import LCFitting

    # Define models for fitting
    models = dict(
        salt_2_4=sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4')),
        salt_2_0=sncosmo.Model(source=sncosmo.get_source('salt2', version='2.0')),
        sn_91bg=sncosmo.Model(source=SN91bgSource())
    )

    # Run fitting
    lc_fitting = LCFitting(args.args_path, args.out_dir)
    fit_func = getattr(lc_fitting, f'fit_{args.survey}')
    fit_func(models=[models[s] for s in args.models],
             num_params=args.num_params,
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
        '-t', '--skip_types',
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

    args = parser.parse_args()
    if args.survey not in ('csp', 'des', 'sdss'):
        raise ValueError(f"Survey name '{args.survey}' not recognized")

    run(args)
