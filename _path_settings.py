import os as _os

SDSS_URL = 'https://data.sdss.org/sas/dr10/boss/papers/supernova/'
FILE_DIR = _os.path.dirname(_os.path.realpath(__file__))
DATA_DIR = _os.path.join(FILE_DIR, 'data')

TEMP_DIR = _os.path.join(DATA_DIR, '.temp/')              # Temporary files
MASTER_PTH = _os.path.join(DATA_DIR, 'master_table.txt')  # Master table path
SMP_DIR = _os.path.join(DATA_DIR, 'smp_data/')            # SMP data files
SNANA_DIR = _os.path.join(DATA_DIR, 'snana_data/')        # SNANA data files
