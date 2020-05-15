#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-


"""This script runs SNID on SDSS Sako et. al 2018 spectra"""

import subprocess
import sys
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


def pre_process(table):
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

    # Remove galaxy spectra
    table = table[table['type'] != 'Gal']

    time = []
    for row in table:
        observed_date = row['date']
        date_with_timezone = observed_date + '+0000'
        date = datetime.strptime(date_with_timezone, '%Y-%m-%d%z')

        unix_time = date.timestamp()
        january_1_1970_in_julian = 2440587.5
        day_in_seconds = 24 * 60 * 60
        time.append((unix_time / day_in_seconds) + january_1_1970_in_julian)

    # Remove spectra outside phase range
    table['time'] = time
    table['phase'] = table['time'] - t0
    table = table[(min_phase <= table['phase']) & (table['phase'] <= max_phase)]

    # Group data by individual spectra
    if table:
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
    obj_ids = spec_summary[spec_summary['Type'] == 'Ia']['CID']
    obj_ids = sorted(obj_ids, key=int)

    for obj_id in tqdm(obj_ids):
        data = sako18spec.get_data_for_id(obj_id)
        processed_data = pre_process(data)
        if processed_data:
            for individual_spectrum in processed_data.groups:
                yield individual_spectrum


def read_snid_output(path):
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

    del data.meta['comments']
    return data


def run_snid(out_dir, spectrum):
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

    snid_input_path = out_dir / f'{obj_id}_{phase:.2f}.dat'
    snid_out_path = snid_input_path.parent / (snid_input_path.stem + '_snid.output')

    # Run SNID
    np.savetxt(snid_input_path, spectrum['wavelength', 'flux'])
    bash_command = f'snid forcez={z} inter=0 plot=0 verbose=0 {snid_input_path}'
    subprocess.Popen(bash_command.split(), cwd=str(out_dir)).communicate()

    # Delete input file
    snid_input_path.unlink()
    if snid_out_path.exists():
        data = read_snid_output(snid_out_path)
        data['obj_id'] = obj_id
        data['phase'] = phase
        return data

    return Table()


def main(out_dir):
    """Run SNID for all SDSS spectra classified by Sako18 as ``Ia``

    Args:
        out_dir (Path): Directory to write results into
    """

    out_dir.mkdir(exist_ok=True, parents=True)

    snid_out = Table()
    for spec in sdss_data_iter():
        snid_out = vstack([snid_out, run_snid(out_dir, spec)])
        snid_out.write(out_dir / 'all.csv', format='ascii.csv', overwrite=True)


if __name__ == '__main__':
    snid_out_dir = Path(__file__).resolve().parent.parent / 'results' / 'snid'
    main(snid_out_dir)
