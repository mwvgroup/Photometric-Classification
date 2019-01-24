#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""Test that our custom model for 91bg like supernovae was ported correctly
for use with SNCosmo.
"""

from unittest import TestCase

from astropy.table import Table

from data_access._utils import keep_restframe_bands


class SN91bgModel(TestCase):

    @classmethod
    def setUpClass(cls):
        test_data = Table()

    def test_data_selection(self):
        keep_restframe_bands()