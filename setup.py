#!/usr/bin/env python3
# -*- coding: UTF-8 -*-


from setuptools import find_packages, setup

setup(name='phot_class',
      packages=find_packages(),
      python_requires='>=3.7',
      setup_requires=['pytest-runner'],
      tests_require=['pytest'],
      package_data={'': ['*.npy']},
      include_pacakge_data=True)
