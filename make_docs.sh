#!/usr/bin/env bash

python setup.py install --user;
cd docs; rm -rf build/;
make html;
cd ../;
rm -rf $HOME/.local/lib/python3.7/site-packages/phot_class-0.0.0-py3.7.egg/
