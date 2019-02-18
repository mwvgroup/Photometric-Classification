#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""Test input data used for fitting light-curves with SNCosmo."""

import itertools
from unittest import TestCase

import numpy as np
from astropy.table import Table

from data_access import des
from data_access import sdss
from data_access._utils import keep_restframe_bands


class BandSelection(TestCase):
    """Test that `keep_restframe_bands` applies the correct data cuts"""

    @classmethod
    def setUpClass(cls):
        """Create a test table with dummy observation ids and band names"""

        cls.band_names = list('ugriz')
        cls.redshift_vals = (0.001, 0.01, 0.1, 0.5)
        cls.lambda_effective = [3500, 4500, 6000, 7500, 8500]
        cls.test_data = Table(data=[range(5), cls.band_names],
                              names=['obs_id', 'band'])

    def test_ug_selection(self):
        """Test that `keep_restframe_bands` applies the correct data cuts"""

        # At a redshift 0 observer frame ug is rest frame ug
        self.test_data.meta['redshift'] = 0
        cut_data = keep_restframe_bands(
            self.test_data, ['u', 'g'], self.band_names, self.lambda_effective)

        self.assertCountEqual(cut_data['band'], ['u', 'g'])

        # At a redshift .3 observer frame ug is rest frame gr
        self.test_data.meta['redshift'] = .3
        cut_data = keep_restframe_bands(
            self.test_data, ['u', 'g'], self.band_names, self.lambda_effective)

        self.assertCountEqual(cut_data['band'], ['g', 'r'])

        # At a redshift .7 observer frame ug is rest frame riz
        self.test_data.meta['redshift'] = .7
        cut_data = keep_restframe_bands(
            self.test_data, ['u', 'g'], self.band_names, self.lambda_effective)

        self.assertCountEqual(cut_data['band'], ['r', 'i', 'z'])


class EmptyInputTables(TestCase):
    """Test for any empty tables when fitting DES and SDSS data"""

    def check_no_empty_tables(self, input_iterable, band_cut):
        """Generic function to check for empty SNCosmo input tables

        Only checks the first 20 tables.

        Args:
            input_iterable (iter): An iterable of SNCosmo input tables
            band_cut  (iter[str]): The bands included in the input table
        """

        msg = 'Empty table for cid {} with bands {}'
        for input_table in itertools.islice(input_iterable, 20):
            cid = input_table.meta['cid']
            self.assertTrue(input_table, msg=msg.format(cid, band_cut))

    def test_empty_des_inputs(self):
        """Test the first 20 DES inputs aren't empty for various band cuts"""

        band_names = ('desg', 'desr', 'desi', 'desz', 'desy')
        band_cuts = (band_names[0: 2], band_names[2:], band_names)

        for band_cut in band_cuts:
            input_tables = des.iter_sncosmo_input(band_cut)
            self.check_no_empty_tables(input_tables, band_cut)

    def test_empty_sdss_inputs(self):
        """Test the first 20 SDSS inputs aren't empty for various band cuts"""

        band_names = ('sdssu', 'sdssg', 'sdssr', 'sdssi', 'sdssz')
        band_cuts = (band_names[0: 2], band_names[2:], band_names)

        for band_cut in band_cuts:
            input_tables = sdss.iter_sncosmo_input(band_cut)
            self.check_no_empty_tables(input_tables, band_cut)


class ZeroPoint(TestCase):
    """Test for correct zero points when fitting DES and SDSS data"""

    def check_iterable(self, input_iterable, expected_zero):
        """Generic function to check zero point of an SNCosmo input tables

        Only checks the first 20 tables.

        Args:
            input_iterable (iter): An iterable of SNCosmo input tables
            expected_zero  (float): The expected zero point
        """

        generic_msg = 'Incorrect zero point for cid {}. Found {}, expected {}'
        for table in itertools.islice(input_iterable, 20):
            cid = table.meta['cid']
            correct_zero = table['zp'] == expected_zero

            if not all(correct_zero):
                bad_indices = np.logical_not(correct_zero)
                example_val = table['zp'][bad_indices][0]

                err_msg = generic_msg.format(cid, example_val, expected_zero)
                self.assertTrue(False, msg=err_msg)

    def test_des_zero_point(self):
        """Test the first 20 DES inputs for a zero point of 27.5"""

        input_iterable = des.iter_sncosmo_input()
        self.check_iterable(input_iterable, 27.5)

    def test_sdss_zero_point(self):
        """Test the first 20 SDSS inputs for a zero point of 25"""

        input_iterable = sdss.iter_sncosmo_input()
        self.check_iterable(input_iterable, 25)
