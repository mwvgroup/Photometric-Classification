#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""This script runs SNID on SDSS Sako et. al 2018 spectra"""

import logging
import subprocess
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table, vstack
from sndata._utils import convert_to_jd
from sndata.sdss import sako18spec

log = logging.getLogger()
formatter = logging.Formatter(f'%(levelname)8s (%(asctime)s): %(name)s - %(message)s')

stream_handler = logging.StreamHandler(sys.stdout)
stream_handler.setFormatter(formatter)
log.addHandler(stream_handler)

file_handler = logging.FileHandler('./snid.log')
file_handler.setFormatter(formatter)
log.addHandler(file_handler)

log.setLevel(logging.DEBUG)

# Add custom project code to python path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Make sure data is downloaded to the local machine
sako18spec.download_module_data()

# Load some data tables from the Sako et al. 2018 publication
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    sdss_master_table = sako18spec.load_table('master').to_pandas(index='CID')

# Specify minimum and maximum phase of spectra to consider
min_phase = -15
max_phase = 15


###############################################################################
# Create an iterator over SDSS spectra formatted for use with SNID
###############################################################################

def get_sdss_t0(obj_id):
    """Get the t0 value for CSP targets

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

    # Remove galaxy spectra and spectra outside phase range
    table = table[table['type'] != 'Gal']
    table['time'] = convert_sdss_date_to_jd(table['date'])
    table['phase'] = table['time'] - t0
    table = table[(min_phase <= table['phase']) & (table['phase'] <= max_phase)]

    # Group data by individual spectra
    if table:  # Raises error for empty table
        table.group_by('phase')
        log.info(f'<{obj_id}> has {len(set(table["phase"]))} spectra')

    else:
	    log.info(f'Formatted table for <{obj_id}> is empty')

    return table


def sdss_data_iter():
    """Iterate over SDSS spectra

    Only includes objects that are spectroscopically confirmed Ia

    Yields:
        Individual spectra as an astropy Table
    """

    # Here we select object Id's for just SNe Ia
    spec_summary = sako18spec.load_table(9)
    obj_ids = spec_summary[spec_summary['Type'] == 'Ia']['CID']  # Todo: Include all SN types?
    obj_ids = sorted(obj_ids, key=int)

    for obj_id in obj_ids:
        log.info(f'Fetching data for <{obj_id}>')

        data = sako18spec.get_data_for_id(obj_id)
        processed_data = format_table(data)
        if processed_data:
            for individual_spectrum in processed_data.groups:
                yield individual_spectrum

        log.debug('Finished with target\n\n')


###############################################################################
# Python wrappers for running SNID
###############################################################################

# Todo: This was a quick, hacky solution with unnecessarily many type casts
def read_snid_types(path):
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
    )

    dataframe = data.to_pandas(index='type')
    types = ['Ia', 'Ib', 'Ic', 'II', 'NotSN']
    num_tempates = dataframe.loc[types].ntemp
    percent_templates = num_tempates / dataframe.loc[types].ntemp.sum()

    out_data = Table(
        [types, num_tempates, percent_templates],
        names=['type', 'ntemp', 'perctemp'])

    return out_data


def run_snid_on_spectrum(out_dir, spectrum, **kwargs):
    """Run SNID on a spectrum

    Uses the SDSS measured redshift.

    Args:
        out_dir   (Path): Directory to write results into
        spectrum (Table): Table with phase, wavelength, and flux columns

    Returns:
        A Table with SNID type results
    """

    obj_id = spectrum.meta['obj_id']
    phase = spectrum['phase'][0]

    # Create input file
    snid_input_path = out_dir / f'{obj_id}_{phase:.2f}.dat'
    snid_out_path = snid_input_path.parent / (snid_input_path.stem + '_snid.output')

    if not snid_out_path.exists():
        # Run SNID
        fortran_kwargs = ' '.join(f'{key}={val}' for key, val in kwargs.items())
        np.savetxt(snid_input_path, spectrum['wavelength', 'flux'])
        bash_command = f'snid {fortran_kwargs} {snid_input_path}'
        log.info(bash_command)
        subprocess.Popen(bash_command.split(), cwd=str(out_dir)).communicate()

        # Delete input file
        snid_input_path.unlink()

    else:
        log.info(f'SNID output already exists at: {snid_out_path}')

    # Parse SNID output file
    # If statement handles failed SNID run with no output file
    if snid_out_path.exists():
        data = read_snid_types(snid_out_path)
        data['obj_id'] = obj_id
        data['phase'] = phase
        data['minwave'] = min(spectrum['wavelength'])
        data['maxwave'] = max(spectrum['wavelength'])
        return data

    return Table()


def run_snid_on_sdss(out_dir, **kwargs):
    """Run SNID for all SDSS spectra classified by Sako18 as ``Ia``

    Args:
        out_dir (Path): Directory to write results into
    """

    out_dir.mkdir(exist_ok=True, parents=True)
    snid_outputs = []
    for spec in sdss_data_iter():
        z = spec.meta['z']
        snid_outputs.append(run_snid_on_spectrum(out_dir, spec, forcez=z, **kwargs))

    out_path = out_dir / 'all.csv'
    vstack(snid_outputs).write(out_path, format='ascii.csv', overwrite=True)

    # The parameter file is overwritten for each target so is of little use
    # We remove it to avoid confusion.
    (out_dir / 'snid.param').unlink()


###############################################################################
# Logic for Sub-typing spectra
###############################################################################


def read_peak_types(path, perc_cutoff=0.):
    """Read the top SNID classifications for each objects

    Args:
        path         (Path): File path to read from
        perc_cutoff (float): Drop values with ``perctemp`` < this value
    """

    types = pd.read_csv(path, engine='python')

    # Create a dummy value that allows us to sort measurements by how close
    # they are to peak phase in asc
    types['psuedo_phase'] = 100 - types.phase.abs()

    # Sort classifications by number of template matches for each class
    types.sort_values(['ntemp', 'psuedo_phase'], ascending=False, inplace=True)

    # Keep only the top classifications for each object
    types.drop_duplicates(keep='first', subset='obj_id', inplace=True)

    # Drop results that don't make the cut
    types = types[types.perctemp >= perc_cutoff]
    types['obj_id'] = types['obj_id'].astype('U100')

    return types.set_index('obj_id')


def combine_typing_results(*paths, perc_cutoff=0.):
    """Combining tables of peak types from SNID

    Values from the later tables take precidence over earlier tables.

    Args:
        *paths       (Path): Path of the types table
        perc_cutoff (float): Drop values with ``perctemp`` < this value

    Returns:
         A pandas data frame with types for each object
    """

    combined_data = read_peak_types(paths[0], perc_cutoff=perc_cutoff)
    for path in paths[1:]:
        new_data = read_peak_types(path, perc_cutoff=perc_cutoff)
        combined_data.update(new_data)

    return combined_data


def run_snid_subtyping_on_sdss(out_dir, object_types, **kwargs):

    out_dir.mkdir(exist_ok=True, parents=True)
    snid_outputs = []
    for spec in sdss_data_iter():
        obj_id = spec.meta['obj_id']

        try:
            type = object_types.loc[obj_id]

        except KeyError:
            log.info(f'No recorded type for <{obj_id}>')
            continue

        run_snid_on_spectrum(out_dir, spec, usetype=type['type'], **kwargs)

    out_path = out_dir / 'all.csv'
    vstack(snid_outputs).write(out_path, format='ascii.csv', overwrite=True)

    # The parameter file is overwritten for each target so is of little use
    # We remove it to avoid confusion.
    (out_dir / 'snid.param').unlink()


if __name__ == '__main__':
    results_dir = Path(__file__).resolve().parent.parent / 'results' / 'snid'

    no_interact_kwargs = dict(inter=0, plot=0, verbose=0)

    type_result_paths = []
    for rlap in (5, 10):
        log.info(f'Typing Spectra rlapmin={rlap}')
        type_path = results_dir / f'type_rlap_{rlap}'
        #run_snid_on_sdss(type_path, rlapmin=rlap, **no_interact_kwargs)
        type_result_paths.append(type_path / 'all.csv')

    log.info('Selecting best fit object types')
    types = combine_typing_results(*type_result_paths, perc_cutoff=.5)

    for rlap in (5, 10):
        log.info(f'Subtyping Spectra rlapmin={rlap}')
        subtype_path = results_dir / f'subtype_rlap_{rlap}'
        run_snid_subtyping_on_sdss(subtype_path, types, rlapmin=rlap, **no_interact_kwargs)
