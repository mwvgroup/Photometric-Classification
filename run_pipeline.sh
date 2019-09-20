#!/usr/bin/env bash

# Run DES SN3YR and fix redshift to the prior
python run_pipeline.py \
    -s des \
    -r sn3yr \
    -f simple_fit \
    -v t0 x0 x1 c \
    -t 120 \
    -c 'priors_config.yml' \
    -o './results';

# Run CSP DR3 and fix redshift to the prior
python run_pipeline.py \
    -s csp \
    -r dr3 \
    -f simple_fit \
    -v t0 x0 x1 c \
    -c 'priors_config.yml' \
    -o './results';
