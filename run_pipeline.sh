#!/usr/bin/env bash

python run_pipeline.py -s csp -r dr3 -f simple_fit -v t0 x0 x1 c -t 90 -o ./phot_class_results
python run_pipeline.py -s des -r sn3yr -f simple_fit -v t0 x0 x1 c -t 90 -o ./phot_class_results
