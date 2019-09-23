#!/usr/bin/env bash

# Run CSP DR3 and fix redshift to the prior
python run_pipeline.py \
    -s csp \
    -r dr3 \
    -f simple_fit \
    -v t0 x0 x1 c \
    -c 'config_files/csp_config.yml' \
    -o './results';

# Run DES SN3YR and fix redshift to the prior
python run_pipeline.py \
    -s des \
    -r sn3yr \
    -f simple_fit \
    -v t0 x0 x1 c \
    -t 120 \
    -c 'config_files/des_config.yml' \
    -o './results';

# Run DES SN3YR and fix redshift to the prior
python run_pipeline.py \
    -s sdss \
    -r sako18 \
    -f simple_fit \
    -v t0 x0 x1 c \
    -t 120 \
    -c 'config_files/sdss_config.yml' \
    -o './results';
