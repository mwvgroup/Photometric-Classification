#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Run the analysis pipeline for all data (CSP, DES, and SDSS) and all models
(SALT 2.4 and 91bg).
"""

import argparse


def run(args):
    """Run light curve fits using command line args"""

    from pathlib import Path

    import sncosmo

    from analysis_pipeline import SN91bgSource
    from analysis_pipeline.lc_fitting import LCFitting

    # Define models for fitting
    salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))
    salt_2_0 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.0'))
    sn_91bg = sncosmo.Model(source=SN91bgSource())

    # Create output directories
    out_dir = Path(__file__).resolve().parent / f'fit_results/{args.survey}'
    out_dir.mkdir(parents=True, exist_ok=True)

    # Run fitting
    lc_fitting = LCFitting(args.args_path)
    fit_func = getattr(lc_fitting, f'fit_{args.survey}')
    fit_func(out_dir, models=[salt_2_4, sn_91bg], num_params=args.num_params)


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
        '-n', '--num_params',
        type=int,
        nargs='+',
        default=[4, 5],
        help='Number of params to fit (4, 5)')

    args = parser.parse_args()
    if args.survey not in ('csp', 'des', 'sdss'):
        raise ValueError(f"Survey name '{args.survey}' not recognized")

    run(args)
