#!/usr/bin/env bash

# # Run DES SN3YR
# python run_pipeline.py \
#     -s des \
#     -r sn3yr \
#     -f simple_fit \
#     -v t0 amplitude x1 c \
#     -c 'config_files/des_config.yml' \
#     -o './results';

# Run SDSS Sako 2018
python run_pipeline.py \
    -s sdss \
    -r sako18 \
    -f simple_fit \
    -c 'config_files/sdss_config.yml' \
    -o './results';