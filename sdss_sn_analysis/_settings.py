#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This document defines constants for use in the sdss_sn_analysis package."""

import os

SDSS_URL = 'https://data.sdss.org/sas/dr10/boss/papers/supernova/'
FILE_DIR = os.path.dirname(os.path.realpath(__file__))

TEMP_DIR = os.path.join(FILE_DIR, '.temp/')              # For temporary files
MASTER_PTH = os.path.join(FILE_DIR, 'master_table.txt')  # Path of master table
SMP_DIR = os.path.join(FILE_DIR, 'smp_data/')            # For SMP data files
SNANA_DIR = os.path.join(FILE_DIR, 'snana_data/')        # For SNANA data files
