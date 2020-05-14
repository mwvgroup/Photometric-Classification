# Example call: snid forcez=0.528 inter=0 plot=0 sn2003jo.dat

import subprocess
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
from astropy.table import Table
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
    """Formats daa tables for use with the GUI

    Changes:
        - Remove galaxy spectra from data tables
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


def main(out_dir):
    """Run SNID for all SDSS spectra classified by Sako18 as ``Ia``

    Args:
        out_dir (Path): Directory to write results into
    """

    out_dir.mkdir(exist_ok=True, parents=True)

    for spectrum in sdss_data_iter():
        obj_id = spectrum.meta['obj_id']
        phase = spectrum['phase'][0]
        z = spectrum.meta['z']

        sndata_input_file = out_dir / f'{obj_id}_{phase:.2f}.dat'
        np.savetxt(sndata_input_file, spectrum['wavelength', 'flux'])

        bash_command = f"snid forcez={z} inter=0 plot=0 verbose=0 {sndata_input_file}"
        process = subprocess.Popen(bash_command.split(), cwd=str(out_dir))
        process.communicate()

        sndata_input_file.unlink()


if __name__ == '__main__':
    snid_out_dir = Path(__file__).resolve().parent.parent / 'results' / 'snid'
    main(snid_out_dir)
