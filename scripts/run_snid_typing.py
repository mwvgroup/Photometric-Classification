#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This script runs SNID on SDSS Sako et. al 2018 spectra

Important: For asserted default SNID arguments, see the
``run_snid_with_defaults`` function.
"""

import logging
import subprocess
import sys
import warnings
from datetime import datetime
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
from astropy.table import Table
from sndata._utils import convert_to_jd
from sndata.sdss import sako18spec

log = logging.getLogger()

# Add custom project code to python path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Make sure data is downloaded to the local machine
sako18spec.download_module_data()

# Load some data tables from the Sako et al. 2018 publication
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    sdss_master_table = sako18spec.load_table('master').to_pandas(index='CID')

# Specify default minimum and maximum phase of spectra to consider
min_phase = -15
max_phase = 15

# Specify default minimum and maximum wavelength range (observer frame) to consider
min_wave = 4000
max_wave = 9000

# Default age range to consider
dage = 10


###############################################################################
# Create an iterator over SDSS spectra formatted for use with SNID
###############################################################################

def get_sdss_t0(obj_id):
    """Get the t0 value for SDSS targets

    Args:
        obj_id (str): The object identifier

    Returns:
        The time of B-band maximum in units of
    """

    log.debug(f'Fetching t0 for <{obj_id}>')

    # Unknown object ID
    if obj_id not in sdss_master_table.index:
        log.error(f't0 not available - unknown object Id <{obj_id}>')
        raise ValueError(f't0 not available for <{obj_id}>')

    t0_mjd = sdss_master_table.loc[obj_id]['PeakMJDSALT2zspec']

    # Known object Id with unknown peak time
    if np.isnan(t0_mjd):
        log.info(f't0 not available for <{obj_id}>')
        raise ValueError(f't0 not available for <{obj_id}>')

    to_jd = convert_to_jd(t0_mjd)
    log.info(f'Found t0 for <{obj_id}>. t0 = {to_jd}')
    return to_jd


@np.vectorize
def convert_sdss_date_to_jd(observed_date):
    """Convert SDSS Spectra observation dates from string format to JD

    Args:
        observed_date (str): Date string with format ``%Y-%m-%d``

    Returns:
        Observed date in JD as a float
    """

    date_with_timezone = observed_date + '+0000'
    date = datetime.strptime(date_with_timezone, '%Y-%m-%d%z')

    unix_time = date.timestamp()
    january_1_1970_in_julian = 2440587.5
    day_in_seconds = 24 * 60 * 60
    date_in_jd = (unix_time / day_in_seconds) + january_1_1970_in_julian

    return date_in_jd


def format_table(table):
    """Formats data tables for use with SNID

    Changes:
        - Removes galaxy spectra
        - Adds Phase column
        - Drops data outside phase range (See globals)
        - Groups data by phase

    Args:
        table (Table): The data to format

    Returns:
        A new table with formatted data
    """

    obj_id = table.meta['obj_id']
    log.debug(f'Formatting table for <{obj_id}>')

    # Get tmax for object
    try:
        t0 = get_sdss_t0(obj_id)

    except ValueError:
        return Table()

    table['time'] = convert_sdss_date_to_jd(table['date'])
    table['phase'] = table['time'] - t0

    # Remove galaxy spectra and spectra outside phase range
    table = table[table['type'] != 'Gal']
    table = table[(table['phase'] >= min_phase) & (table['phase'] <= max_phase)]

    if table:
        table = table.group_by('phase')
        log.info(f'Formatted table for <{obj_id}> has {len(table.groups)} spectra')

    else:
        log.info(f'Formatted table for <{obj_id}> has 0 spectra')

    return table


def sdss_data_iter():
    """Iterate over SDSS spectra

    Includes all spectroscopically observed targets listed in Sako et al.
    2018 Table 9.

    Yields:
        Individual spectra as an astropy Table
    """

    spec_summary = sako18spec.load_table(9)
    obj_ids = sorted(spec_summary['CID'], key=int)

    # The spec_summary table multiple entries per object. We only want the unique ID's
    unique_ids = set(obj_ids)
    total_ids = len(unique_ids)
    for i, obj_id in enumerate(unique_ids):
        log.info(f'Fetching data for <{obj_id}> ({i} / {total_ids})')

        data = sako18spec.get_data_for_id(obj_id)
        processed_data = format_table(data)
        if processed_data:
            for spec in processed_data.groups:
                yield spec

        log.debug('Finished with target\n\n')


###############################################################################
# Read output files from SNID (sub)typing
###############################################################################


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
    type_names = ['Ia', 'Ib', 'Ic', 'II', 'NotSN']
    peak_type = data.loc[type_names].ntemp.idxmax()
    ntemp = data.loc[peak_type].ntemp
    total_templates = data.loc[type_names].ntemp.sum()
    percent_templates = ntemp / total_templates * 100

    return peak_type, ntemp, percent_templates


def compile_peak_types(snid_dir, perc_cutoff=0.5):
    """Get peak types from previous SNID run

    Args:
        snid_dir     (Path): Directory of SNID outputs
        perc_cutoff (float): Only return results best type has >=
            perc_cutoff of the template matches

    Returns:
        A DataFrame indexed by object ID
    """

    if not snid_dir.exists():
        raise FileNotFoundError(f'Results directory DNE: {snid_dir}')

    rows = []
    for path in snid_dir.glob('*snid.output'):
        obj_id, phase, *_ = path.name.split('_')
        peak_type, ntemp, percent_templates = read_peak_type(path)
        rows.append([int(obj_id), float(phase), peak_type, ntemp, percent_templates])

    type_data = pd.DataFrame(
        rows,
        columns=['objId', 'Phase', 'Type', 'nType', 'percType', ])

    # Keep only the spectra nearest peak
    type_data['abs_phase'] = type_data.Phase.abs()
    type_data = type_data.sort_values('abs_phase', ascending=True)
    type_data = type_data.drop_duplicates(keep='first', subset='objId')
    type_data = type_data.sort_values('objId', ascending=True)

    # Apply percentage cutoff
    type_data = type_data[type_data.percType >= perc_cutoff]
    return type_data.set_index('objId')


def read_peak_subtype(path):
    """Return the subtype summary from an SNID output file

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

    # Get the subtype with the most matches. Make sure the subtype
    # is not equally as "good" a match as the second best subtype
    sn_type, subtype, second_subtype = data.ntemp.nlargest(3).index
    assert data.loc[subtype].ntemp != second_subtype

    subtype_record = data.loc[subtype]
    type_record = data.loc[sn_type]
    perc_temp = subtype_record.ntemp / type_record.ntemp * 100

    return subtype_record.redshift, subtype, subtype_record.ntemp, perc_temp


def compile_peak_subtypes(snid_dir, perc_cutoff=0.5):
    """Get peak subtypes from all output files from a previous SNID run

    Args:
        snid_dir (Path): Directory of SNID outputs
        perc_cutoff (float): Only return results best subtype has >=
            perc_cutoff of the template matches

    Returns:
        A DataFrame indexed by object ID
    """

    if not snid_dir.exists():
        raise FileNotFoundError(f'Results directory DNE: {snid_dir}')

    rows = []
    for path in snid_dir.glob('*snid.output'):
        obj_id, phase, *_ = path.name.split('_')
        rows.append([int(obj_id), float(phase), *read_peak_subtype(path)])

    type_data = pd.DataFrame(
        rows,
        columns=['objId', 'Phase', 'Redshift', 'subType', 'nSubType', 'percSubType'])

    # Keep only the spectra nearest peak
    type_data['abs_phase'] = type_data['Phase'].abs()
    type_data = type_data.sort_values('abs_phase', ascending=True)
    type_data = type_data.drop_duplicates(keep='first', subset='objId')
    type_data = type_data.drop('abs_phase', axis='columns')

    type_data = type_data.sort_values('objId')

    type_data = type_data[type_data.percSubType >= perc_cutoff]
    return type_data.set_index('objId')


###############################################################################
# Python wrappers for running SNID
###############################################################################


def run_snid_on_spectrum(out_dir, spectrum, **kwargs):
    """Run SNID on a spectrum

    Args:
        out_dir   (Path): Directory to write results into
        spectrum (Table): Table with phase, wavelength, and flux columns
        Any additional SNID kwargs
    """

    obj_id = spectrum.meta['obj_id']
    phase = spectrum['phase'][0]
    if not len(set(spectrum['phase'])) == 1:
        raise ValueError('SNID passed multiple spectra')

    # Create input file in temporary directory
    with TemporaryDirectory() as temp_dir:
        snid_input_path = Path(temp_dir) / f'{obj_id}_{phase:.2f}.dat'
        np.savetxt(snid_input_path, spectrum['wavelength', 'flux'])

        # Run SNID
        fortran_kwargs = ' '.join(f'{key}={val}' for key, val in kwargs.items())
        bash_command = f'snid {fortran_kwargs} {snid_input_path}'
        log.info(bash_command)
        subprocess.Popen(bash_command.split(), cwd=str(out_dir)).communicate()


def run_snid_with_defaults(out_dir, spec, **kwargs):
    """Run SNID on a spectrum with customized defaults

    Asserts SNID arguments:
        - forcez = SDSS estimated redshift
        - emclip = SDSS estimated redshift
        - age = SDSS estimated phase
        - wmin = ``min_wave`` global
        - wmax = ``max_wave`` global
        - dage = ``dage`` global

    Args:
        out_dir   (Path): Directory to write results into
        spectrum (Table): Table with phase, wavelength, and flux columns
        Any additional SNID kwargs
    """

    z = spec.meta['z']
    phase = spec['phase'][0]
    run_snid_on_spectrum(
        out_dir,
        spec,
        forcez=z,
        emclip=z,
        age=phase,
        wmin=min_wave,
        wmax=max_wave,
        dage=dage,
        **kwargs)


def run_snid_typing_on_sdss_with_defaults(out_dir, **kwargs):
    """Run SNID for all SDSS spectra

    For a list of customized default arguments see
    ``run_snid_with_defaults``.

    Args:
        out_dir (Path): Directory to write results into
        Any additional SNID kwargs

    Returns:
        A dataframe with typing results for all targets with >50% temlate matches
    """

    out_dir.mkdir(exist_ok=True, parents=True)
    for spec in sdss_data_iter():
        run_snid_with_defaults(out_dir, spec, **kwargs)

    # The parameter file is overwritten for each target so is of little use
    # We remove it to avoid confusion.
    (out_dir / 'snid.param').unlink()

    return compile_peak_types(out_dir)


def run_snid_subtyping_on_sdss(out_dir, object_types, **kwargs):
    """Run SNID on SDSS using a single type per object

    ``object_types`` must have a column ``Type`` and be indexed by object id.

    Args:
        out_dir (Path): Directory to write results into
        object_types (DataFrame): The types to use for each object

    Returns:
         A pandas data frame with types for each object
    """

    out_dir.mkdir(exist_ok=True, parents=True)
    for spec in sdss_data_iter():
        obj_id = int(spec.meta['obj_id'])

        try:
            type_record = object_types.loc[obj_id]

        except KeyError:
            log.info(f'No recorded type for <{obj_id}>')
            continue

        template_type = type_record['Type']
        log.info(f'Using template {template_type}')
        run_snid_with_defaults(out_dir, spec, usetype=template_type, **kwargs)

    # The parameter file is overwritten for each target so is of little use
    # We remove it to avoid confusion.
    (out_dir / 'snid.param').unlink()


def construct_summary_table(results_dir):
    # Get subtyping results
    subtypes = compile_peak_subtypes(results_dir / 'subtype_rlap_5')
    subtypes.update(compile_peak_subtypes(results_dir / 'subtype_rlap_10'))

    # Add typing results to the dataframe
    types = compile_peak_types(results_dir / 'type_rlap_5')
    types.update(compile_peak_types(results_dir / 'type_rlap_10'))
    subtypes = subtypes.join(types[['nType', 'percType']])

    # Add SDSS classifications
    master = sako18spec.load_table('master')
    master['objId'] = np.array(master['CID'], dtype=float)
    master['SDSSClass'] = master['Classification']
    sdss_df = master.to_pandas('objId')
    subtypes = subtypes.join(sdss_df.SDSSClass, on='objId')

    return subtypes


if __name__ == '__main__':

    _results_dir = Path(__file__).resolve().parent / 'results' / f'snid_dage{dage}'
    _results_dir.mkdir(exist_ok=True, parents=True)

    formatter = logging.Formatter(f'%(levelname)8s (%(asctime)s): %(name)s - %(message)s')
    file_handler = logging.FileHandler(_results_dir / 'snid_typing.log')
    file_handler.setFormatter(formatter)
    log.addHandler(file_handler)

    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    log.addHandler(stream_handler)
    log.setLevel(logging.DEBUG)

    # Run typing
    log.info('Typing Spectra rlapmin=5')
    type_path = _results_dir / f'type_rlap_5'
    types_5 = run_snid_typing_on_sdss_with_defaults(type_path, rlapmin=5, inter=0, plot=0, verbose=0)

    log.info('Typing Spectra rlapmin=10')
    type_path = _results_dir / f'type_rlap_10'
    types_10 = run_snid_typing_on_sdss_with_defaults(type_path, rlapmin=10, inter=0, plot=0, verbose=0)

    # Combine typing results from rlap=5 and rlap=10.
    # Prefer rlap=10 results
    final_types = types_5.copy()
    final_types.update(types_10)

    for rlap in (5, 10):
        log.info(f'Subtyping Spectra rlapmin={rlap}')
        subtype_path = _results_dir / f'subtype_rlap_{rlap}'
        run_snid_subtyping_on_sdss(subtype_path, final_types, rlapmin=rlap, inter=0, plot=0, verbose=0)

    log.info('Writing summary table...')
    summary_table = construct_summary_table(_results_dir)
    summary_table.to_csv(_results_dir / 'summary.csv')
