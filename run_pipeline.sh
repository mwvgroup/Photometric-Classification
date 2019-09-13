#!/usr/bin/env bash

# Run CSP DR3 and fix redshift to the published value
python run_pipeline.py \
    -s csp \
    -r dr3 \
    -f simple_fit \
    -v t0 x0 x1 c \
    -c 'config_files/csp.yml' \
    -o './results';
