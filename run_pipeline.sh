#!/usr/bin/env bash

# Run DES SN3YR
# python run_pipeline.py \
#     -s des \
#     -r sn3yr \
#     -f mcmc_fit \
#     -m band \
#     -c 'config_files/des_config_noext.yml' \
#     -o './results/band_fits/no_ext';

python run_pipeline.py \
    -s des \
    -r sn3yr \
    -f mcmc_fit \
    -m band \
    -c 'config_files/des_config_ext.yml' \
    -o './results/band_fits/with_ext';

# python run_pipeline.py \
#     -s des \
#     -r sn3yr \
#     -f mcmc_fit \
#     -m collective \
#     -c 'config_files/des_config_noext.yml' \
#     -o './results/collective_fits/no_ext';

python run_pipeline.py \
    -s des \
    -r sn3yr \
    -f mcmc_fit \
    -m collective \
    -c 'config_files/des_config_ext.yml' \
    -o './results/collective_fits/with_ext';

# Run SDSS Sako 2018
# python run_pipeline.py \
#     -s sdss \
#     -r sako18 \
#     -f mcmc_fit \
#     -m band \
#     -c 'config_files/sdss_config_noext.yml' \
#     -o './results/band_fits/no_ext';

python run_pipeline.py \
    -s sdss \
    -r sako18 \
    -f mcmc_fit \
    -m band \
    -c 'config_files/sdss_config_ext.yml' \
    -o './results/band_fits/with_ext';

python run_pipeline.py \
    -s sdss \
    -r sako18 \
    -f mcmc_fit \
    -m collective \
    -c 'config_files/sdss_config_ext.yml' \
    -o './results/collective_fits/with_ext';