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
        cls.band_names = list('ugriz')
        cls.redshift_vals = (0.001, 0.01, 0.1, 0.5)
        cls.lambda_effective = [3500, 4500, 6000, 7500, 8500]
        cls.test_data = Table(data=[range(5), cls.band_names],
                              names=['obs_id', 'band'])

    def test_ug_selection(self):

        self.test_data.meta['redshift'] = 0
        cut_data = keep_restframe_bands(
            self.test_data, ['u', 'g'], self.band_names, self.lambda_effective)

        self.assertItemsEqual(cut_data['band'], ['u', 'g'])

        self.test_data.meta['redshift'] = .3
        cut_data = keep_restframe_bands(
            self.test_data, ['u', 'g'], self.band_names, self.lambda_effective)

        self.assertItemsEqual(cut_data['band'], ['g', 'r'])

        self.test_data.meta['redshift'] = .7
        cut_data = keep_restframe_bands(
            self.test_data, ['u', 'g'], self.band_names, self.lambda_effective)

        self.assertItemsEqual(cut_data['band'], ['r', 'i', 'z'])
