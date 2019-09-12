#!/usr/bin/env bash

python run_pipeline.py \
    -s csp \
    -r dr3 \
    -f simple_fit \
    -v t0 x0 x1 c \
    -bga 'config_files/csp_sn91bg.yml' \
    -s2a 'config_files/csp_salt2.yml' \
    -o './results';
