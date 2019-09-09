#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``fitting`` module."""

from unittest import TestCase

import numpy as np
import sncosmo
from astropy.table import Table

from phot_class import classification


# Todo: Test the following functions
# - run_fits
# - tabulate_fit_results


class TestTableCreation(TestCase):
    """Tests for fitting.create_empty_table"""

    def test_is_empty(self):
        """Test the returned table is empty by default"""

        num_rows = len(classification.create_empty_table())
        self.assertEqual(0, num_rows, 'Table is not empty')

    def test_correct_columns(self):
        """Test the returned table has the correct columns"""

        returned = classification.create_empty_table().colnames
        expected = [
            'obj_id', 'band_set', 'source', 'pre_max', 'post_max',
            'z', 't0', 'x0', 'x1', 'c',
            'z_err', 't0_err', 'x0_err', 'x1_err', 'c_err',
            'chisq', 'ndof', 'b_max', 'delta_15', 'message', ]

        self.assertSequenceEqual(expected, returned)

    def test_is_masked(self):
        """Test the returned table is a masked Table by default"""

        is_masked = classification.create_empty_table().masked
        self.assertTrue(is_masked, 'Table has no mask')

    def test_accepts_kwargs(self):
        """Test that the table constructor uses the kwargs"""

        num_columns = len(classification.create_empty_table().colnames)
        dummy_row = np.ones(num_columns)
        table = classification.create_empty_table(rows=[dummy_row])
        self.assertEqual(1, len(table), 'Table did not add specified rows')


class TestFitsToTableRow(TestCase):
    """Tests for fitting.fit_results_to_table_row"""

    def runTest(self):
        """Run fits for data with a known result and check the returned
        row is correct.
        """

        # We use the example data from sncosmo since it is weel behaved
        data = sncosmo.load_example_data()
        model = sncosmo.Model('salt2')
        model.update(data.meta)
        data.meta['obj_id'] = 'dummy_id'

        fit_results, fitted_model = sncosmo.fit_lc(
            data=data,
            model=model,
            vparam_names=['z', 't0', 'x0', 'x1', 'c'],
            bounds={'z': (0.3, 0.7)})

        row = classification.fit_results_to_table_row(
            data, 'dummy_band_set', fit_results, fitted_model)

        expected_row = [
            'dummy_id',  # obj_id
            'dummy_band_set',  # band_set
            'salt2',  # source
            15,  # pre_max
            25,  # post_max
            0.5151630264692924,  # z
            55100.47790767062,  # t0
            1.1962747082468881e-05,  # x0
            0.46680671979205146,  # x1
            0.19392275740112822,  # c
            0.01541479638843532,  # z_err
            0.41542130042944336,  # t0_err
            3.985676691503247e-07,  # x0_err
            0.3464229401925888,  # x1_err
            0.03751626613909028,  # c_err
            3306.312724429728,  # chisq
            39,  # ndof
            -19.565602274514767,  # b_max
            0.9593115374152568,  # delta_15
            'Minimization exited successfully.'  # message
        ]

        self.assertCountEqual(expected_row, row)


class TestClassificationCoords(TestCase):
    """Tests for fitting.classify_targets"""

    expected_input_columns = ['obj_id', 'source', 'band_set', 'chisq', 'ndof']

    def test_correct_coordinates(self):
        """Test correct coordinates are returned for the given input data"""

        test_data = Table(names=self.expected_input_columns, rows=[
            ['dummy_id', 'salt2', 'blue', 10, 1],
            ['dummy_id', 'salt2', 'red', 20, 1],
            ['dummy_id', 'sn91bg', 'blue', 10, 1],
            ['dummy_id', 'sn91bg', 'red', 10, 1]
        ])

        expected_row = ['dummy_id', 0, 10]
        class_coordinates = classification.classify_targets(test_data)
        self.assertListEqual(list(class_coordinates[0]), expected_row)

    def test_masked_values(self):
        """Test handling of masked data"""

        test_data = Table(names=self.expected_input_columns, rows=[
            ['dummy_id', 'salt2', 'blue', 10, 1],
            ['dummy_id', 'salt2', 'red', 20, 1],  # Second row
            ['dummy_id', 'sn91bg', 'blue', 10, 1],
            ['dummy_id', 'sn91bg', 'red', 10, 1]
        ], masked=True)

        # Mask data in the second row
        test_data.mask = [
            [False, True, False, False],
            [False, True, False, False],
            [False, True, False, False],
            [False, True, False, False],
            [False, True, False, False]
        ]

        class_coordinates = classification.classify_targets(test_data)
        self.assertEqual(0, len(class_coordinates))

    def test_failed_fits(self):
        """Test handling of rows with NAN values"""

        # Instantiate a table representing a failed 91bg fit
        test_data = Table(names=self.expected_input_columns, rows=[
            ['dummy_id', 'salt2', 'blue', 10, 1],
            ['dummy_id', 'salt2', 'red', 20, 1],
            ['dummy_id', 'sn91bg', 'blue', np.NAN, np.NAN],
            ['dummy_id', 'sn91bg', 'red', np.NAN, np.NAN]
        ])

        class_coordinates = classification.classify_targets(test_data)
        self.assertEqual(0, len(class_coordinates))
