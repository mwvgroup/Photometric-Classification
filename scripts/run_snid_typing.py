#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-


"""This script runs SNID on SDSS Sako et. al 2018 spectra"""

import subprocess
import sys
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack
from sndata._utils import convert_to_jd
from sndata.sdss import sako18spec
from tqdm import tqdm

# Add custom project code to python path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Make sure data is downloaded to the local machine
sako18spec.download_module_data()

# Load some data tables from the Sako et al. 2018 publication
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    sdss_master_table = sako18spec.load_table('master').to_pandas(index='CID')

# Specify minimum and maximum phase to include in returned data (inclusive)
min_phase = -15
max_phase = 15


def get_sdss_t0(obj_id):
    """Get the t0 value for CSP targets

    Args:
        obj_id (str): The object identifier

    Returns:
        The time of B-band maximum in units of
    """

    # Unknown object ID
    if obj_id not in sdss_master_table.index:
        raise ValueError(f't0 not available for {obj_id}')

    t0_mjd = sdss_master_table.loc[obj_id]['PeakMJDSALT2zspec']

    # Known object Id with unknown peak time
    if np.isnan(t0_mjd):
        raise ValueError(f't0 not available for {obj_id}')

    return convert_to_jd(t0_mjd)


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
    return (unix_time / day_in_seconds) + january_1_1970_in_julian


def format_table(table):
    """Formats data tables for use with SNID

    Args:
        table (Table): The data to format

    Returns:
        A new table with formatted data
    """

    obj_id = table.meta['obj_id']

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

    return table


def sdss_data_iter():
    """Iterate over SDSS spectra

    Only includes objects that are spectroscopically confirmed Ia

    Yields:
        Individual spectra as an astropy Table
    """

    # Here we select object Id's for just SNe Ia
    spec_summary = sako18spec.load_table(9)
    obj_ids = spec_summary[spec_summary['Type'] == 'Ia']['CID']  # Todo: Include all SN types
    obj_ids = sorted(obj_ids, key=int)

    for obj_id in tqdm(obj_ids, desc='Object Ids'):
        data = sako18spec.get_data_for_id(obj_id)
        processed_data = format_table(data)
        if processed_data:
            for individual_spectrum in processed_data.groups:
                yield individual_spectrum


# Todo: This was a quick / hacky solution with unnecessarily many type casts
def read_snid_types(path):
    """Return type summary from an SNID output file

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


def run_snid_on_spectrum(out_dir, spectrum, inter=0, plot=0, verbose=0, **kwargs):
    """Run SNID on a spectrum

    Use the SDSS measured redshift.

    Args:
        out_dir   (Path): Directory to write results into
        spectrum (Table): Table with phase, wavelength, and flux columns

    Returns:
        A Table with SNID type results
    """

    obj_id = spectrum.meta['obj_id']
    phase = spectrum['phase'][0]
    z = spectrum.meta['z']

    # Create input file
    snid_input_path = out_dir / f'{obj_id}_{phase:.2f}.dat'
    snid_out_path = snid_input_path.parent / (snid_input_path.stem + '_snid.output')

    if not snid_out_path.exists():
        # Run SNID
        fortran_kwargs = f'inter={inter} plot={plot} verbose={verbose}'
        for key, val in kwargs.items():
            fortran_kwargs += f' {key}={val}'

        np.savetxt(snid_input_path, spectrum['wavelength', 'flux'])
        bash_command = f'snid forcez={z} {fortran_kwargs} {snid_input_path}'
        subprocess.Popen(bash_command.split(), cwd=str(out_dir)).communicate()

        # Delete input file
        snid_input_path.unlink()

    # Parse SNID output file
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
    snid_outputs = [run_snid_on_spectrum(out_dir, spec, **kwargs) for spec in sdss_data_iter()]

    out_path = out_dir / 'all.csv'
    vstack(snid_outputs).write(out_path, format='ascii.csv', overwrite=True)

    # The parameter file is overwritten for each target so is of little use
    # We remove it to avoid confusion.
    (out_dir / 'snid.param').unlink()


if __name__ == '__main__':
    results_dir = Path(__file__).resolve().parent.parent / 'results' / 'snid'

    print('Typing Spectra rlapmin=10')
    run_snid_on_sdss(results_dir / 'type_rlap_10', rlapmin=10)

    print('Typing Spectra rlapmin=5')
    run_snid_on_sdss(results_dir / 'type_rlap_5', rlapmin=5)
