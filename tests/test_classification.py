#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``fitting`` module."""

from unittest import TestCase

import numpy as np
import sncosmo
from astropy.table import Table
from sncosmo.utils import Result
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


class TestFitResultsToTableRow(TestCase):
    """Tests for fitting.fit_results_to_table_row"""

    # Don't limit output messages on test failures
    maxDiff = None

    def runTest(self):
        """Run fits for data with a known result and check the results are
        correctly formatted as a table row.
        """

        # Define a mock "fitted model"
        params = ['z', 't0', 'x0', 'x1', 'c']
        model_params = {'x1': 0.5, 'c': 0.2, 'z': 0.5, 'x0': 1.0, 't0': 55100.0}
        model = sncosmo.Model('salt2')
        model.update(model_params)

        # Define mock fit results for a model fit for the last four parameters
        # I.e. for a fit where z is specified
        result = Result({
            'param_names': params,
            'vparam_names': params[1:],
            'parameters': [model_params[p] for p in params],
            'errors': {p: .1 * model_params[p] for p in params}
        })

        # Use demo data as mock data
        # We drop a known data point that is out of the phase range
        data = sncosmo.load_example_data()
        data = data[data['time'] > model_params['t0'] + 49]
        data.meta['obj_id'] = 'dummy_id'

        row = classification.fit_results_to_table_row(
            data, 'dummy_band_set', result, model)

        expected_row = [
            'dummy_id',  # obj_id
            'dummy_band_set',  # band_set
            'salt2',  # source
            15,  # pre_max
            25,  # post_max
            model_params['z'],
            model_params['t0'],
            model_params['x0'],
            model_params['x1'],
            model_params['c'],
            result.errors['z'],
            result.errors['t0'],
            result.errors['x0'],
            result.errors['x1'],
            result.errors['c'],
            11844.58,  # chisq
            len(data) - len(result.vparam_names),  # ndof
            -31.79,  # b_max
            0.953,  # delta_15
            'NONE'  # message
        ]

        self.assertListEqual(expected_row, row)


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
