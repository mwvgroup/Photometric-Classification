#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This script runs SNID on SDSS Sako et. al 2018 spectra"""

import logging
import sys
from pathlib import Path

import pandas as pd
from astropy.table import Table
from run_snid_typing import run_snid_on_spectrum, sdss_data_iter

log = logging.getLogger()
log.setLevel(logging.DEBUG)

# Parent types for SNID templates used in the classification
# E.g. 'Ia' is a parent type, while 'Ia-norm' and 'Ia-91bg' are subtypes
TYPES = ['Ia', 'Ib', 'Ic', 'II', 'NotSN']


def read_peak_type(path):
    """Return the type summary from an SNID output file

    Args:
        path (str, Path): Path to read

    Returns:
         An astropy Table
    """

    names = ['type', 'ntemp', 'fraction', 'slope', 'redshift',
             'redshift_error', 'age', 'age_error']

    data = Table.read(
        str(path), header_start=4, data_start=4,
        data_end=28, format='ascii.basic', names=names
    ).to_pandas(index='type')

    # Calculate percentage of templates used for each type
    # Total matched templates equals the sum of matches for the parent types
    # (see TYPES global)
    peak_type = data.loc[TYPES].ntemp.idxmax()
    total_templates = data.loc[TYPES].ntemp.sum()
    percent_templates = data.loc[peak_type].ntemp / total_templates

    return peak_type, percent_templates


def compile_peak_types(results_dir, perc_cutoff=0.5):
    """Get peak types from previous SNID run

    Args:
        results_dir (Path): Directory of SNID outputs
        perc_cutoff (float): Only return results where most likley type has >=
            perc_cutoff of the template matches

    Returns:
        A DataFrame indexed by object ID
    """

    rows = []
    for path in results_dir.glob('*snid.output'):
        obj_id, phase, *_ = path.name.split('_')
        peak_type, percent_templates = read_peak_type(path)
        rows.append([obj_id, float(phase), peak_type, percent_templates])

    type_data = pd.DataFrame(
        rows,
        columns=['obj_id', 'phase', 'peak_type', 'perc_temp'])

    # Keep only the spectra nearest peak
    type_data['abs_phase'] = type_data.phase.abs()
    type_data = type_data.sort_values('abs_phase', ascending=True)
    type_data = type_data.drop_duplicates(keep='first', subset='obj_id')

    # Apply percentage cutoff
    type_data = type_data[type_data.perc_temp >= perc_cutoff]

    type_data['obj_id'] = type_data['obj_id'].astype('str')
    return type_data.set_index('obj_id')


def combine_typing_results(*paths, perc_cutoff=0.5):
    """Combining tables of peak types from SNID

    Values from the later tables take precedence over earlier tables.

    Args:
        *paths       (Path): Path of the types table
        perc_cutoff (float): Drop values with ``perctemp`` < this value

    Returns:
         A pandas data frame with types for each object
    """

    combined_data = None
    for path in paths:
        if combined_data is None:
            combined_data = compile_peak_types(path, perc_cutoff=perc_cutoff)
            continue

        new_data = compile_peak_types(path, perc_cutoff=perc_cutoff)
        combined_data.update(new_data)

    return combined_data


def run_snid_subtyping_on_sdss(out_dir, object_types, **kwargs):
    """Run SNID on SDSS using a single type per object

    ``object_types`` must have a column ``type`` and be indexed by object id.

    Args:
        out_dir (Path): Directory to write results into
        object_types (DataFrame): The types to use for each object

    Returns:
         A pandas data frame with types for each object
    """

    out_dir.mkdir(exist_ok=True, parents=True)
    for spec in sdss_data_iter():
        obj_id = spec.meta['obj_id']

        try:
            type_record = object_types.loc[obj_id]

        except KeyError:
            log.info(f'No recorded type for <{obj_id}>')
            continue

        template_type = type_record['peak_type']
        log.info(f'Using template {template_type}')
        run_snid_on_spectrum(out_dir, spec, usetype=template_type, **kwargs)

    # The parameter file is overwritten for each target so is of little use
    # We remove it to avoid confusion.
    (out_dir / 'snid.param').unlink()


if __name__ == '__main__':
    formatter = logging.Formatter(f'%(levelname)8s (%(asctime)s): %(name)s - %(message)s')

    file_handler = logging.FileHandler('./snid_subtyping.log')
    file_handler.setFormatter(formatter)
    log.addHandler(file_handler)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    log.addHandler(stream_handler)

    results_dir = Path(__file__).resolve().parent.parent / 'results' / 'snid'
    rlap5_type_paths = results_dir / f'type_rlap_5'
    rlap10_type_paths = results_dir / f'type_rlap_10'

    log.info('Selecting best fit object types')
    types = combine_typing_results(rlap5_type_paths, rlap10_type_paths, perc_cutoff=.5)

    for rlap in (5, 10):
        log.info(f'Subtyping Spectra rlapmin={rlap}')
        subtype_path = results_dir / f'subtype_rlap_{rlap}'
        run_snid_subtyping_on_sdss(subtype_path, types, rlapmin=rlap, inter=0, plot=0, verbose=0)

    log.info(f'Subtyping Spectra with exlusively rlapmin=10 results')
    strict_dir = results_dir / 'subtype_rlap_strict'
    types = combine_typing_results(rlap10_type_paths, perc_cutoff=.5)
    run_snid_subtyping_on_sdss(strict_dir, types, rlapmin=10, inter=0, plot=0, verbose=0)
