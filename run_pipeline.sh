#!/usr/bin/env bash

# Spectroscopic classification
##############################

# First we fixed the bin size and vary the number of samples
python run_pipeline.py \
    -s sdss \
    -r sako18spec \
    spectroscopic \
    -n 0 \
    -o './results/';

python run_pipeline.py \
    -s sdss \
    -r sako18spec \
    spectroscopic \
    -n 2 \
    -o './results/';

python run_pipeline.py \
    -s sdss \
    -r sako18spec \
    spectroscopic \
    -n 5 \
    -o './results/';

# We also try an alternate bin size
python run_pipeline.py \
    -s sdss \
    -r sako18spec \
    spectroscopic \
    -n 5 \
    -b 3 \
    -o './results/';

# emcee fitting
###############

python run_pipeline.py \
    -s sdss \
    -r sako18 \
    photometric \
    -f mcmc_fit \
    -m band \
    -c 'config_files/sdss_config_ext.yml' \
    -o './results/band_fits/with_ext';

python run_pipeline.py \
    -s sdss \
    -r sako18 \
    photometric \
    -f mcmc_fit \
    -m collective \
    -c 'config_files/sdss_config_ext.yml' \
    -o './results/collective_fits/with_ext';

python run_pipeline.py \
    -s des \
    -r sn3yr \
    photometric \
    -f mcmc_fit \
    -m band \
    -c 'config_files/des_config_ext.yml' \
    -o './results/band_fits/with_ext';

python run_pipeline.py \
    -s des \
    -r sn3yr \
    photometric \
    -f mcmc_fit \
    -m collective \
    -c 'config_files/des_config_ext.yml' \
    -o './results/collective_fits/with_ext';

# Iminuit Fitting
#################

# Run DES SN3YR

python run_pipeline.py \
    -s des \
    -r sn3yr \
    photometric \
    -f simple_fit \
    -m band \
    -c 'config_files/des_config_ext.yml' \
    -o './results/band_fits/with_ext';

python run_pipeline.py \
    -s des \
    -r sn3yr \
    photometric \
    -f simple_fit \
    -m collective \
    -c 'config_files/des_config_ext.yml' \
    -o './results/collective_fits/with_ext';

# Run SDSS Sako 2018

python run_pipeline.py \
    -s sdss \
    -r sako18 \

    photometric \
    -f simple_fit \
    -m band \
    -c 'config_files/sdss_config_ext.yml' \
    -o './results/band_fits/with_ext';

python run_pipeline.py \
    -s sdss \
    -r sako18 \
    photometric \
    -f simple_fit \
    -m collective \
    -c 'config_files/sdss_config_ext.yml' \
    -o './results/collective_fits/with_ext';

# A repeat of the above, but without extinction
#################

python run_pipeline.py \
    -s des \
    -r sn3yr \
    photometric \
    -f simple_fit \
    -m band \
    -c 'config_files/des_config_noext.yml' \
    -o './results/band_fits/no_ext';

python run_pipeline.py \
    -s des \
    -r sn3yr \
    photometric \
    -f simple_fit \
    -m collective \
    -c 'config_files/des_config_noext.yml' \
    -o './results/collective_fits/no_ext';

# python run_pipeline.py \
#     -s sdss \
#     -r sako18 \
#     photometric \
#     -f simple_fit \
#     -m band \
#     -c 'config_files/sdss_config_noext.yml' \
#     -o './results/band_fits/no_ext';

# python run_pipeline.py \
#     -s sdss \
#     -r sako18 \
#     photometric \
#     -f collective \
#     -c 'config_files/sdss_config_noext.yml' \
#     -o './results/band_fits/no_ext';