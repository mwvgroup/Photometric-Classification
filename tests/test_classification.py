#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``classification`` module."""

from unittest import TestCase

import numpy as np
import sncosmo
from astropy.table import Table
from sncosmo.utils import Result

from phot_class import classification


# Todo: Test the following functions
# - run_classification_fits
# - tabulate_fit_results


class TableCreation(TestCase):
    """Tests for the ``create_empty_table`` function"""

    def test_is_empty(self):
        """Test the returned table is empty by default"""

        num_rows = len(classification.create_empty_table([]))
        self.assertEqual(0, num_rows, 'Table is not empty')

    def test_correct_columns(self):
        """Test the returned table has the correct columns"""

        parameters = ['z', 't0', 'x0', 'x1', 'c']
        returned = classification.create_empty_table(parameters).colnames
        expected = [
            'obj_id', 'band', 'source', 'pre_max', 'post_max',
            'z', 't0', 'x0', 'x1', 'c',
            'z_err', 't0_err', 'x0_err', 'x1_err', 'c_err',
            'chisq', 'ndof', 'b_max', 'delta_15', 'message', ]

        self.assertSequenceEqual(expected, returned)

    def test_is_masked(self):
        """Test the returned table is a masked Table by default"""

        is_masked = classification.create_empty_table([]).masked
        self.assertTrue(is_masked, 'Table has no mask')

    def test_accepts_kwargs(self):
        """Test that the table constructor uses the kwargs"""

        num_columns = len(classification.create_empty_table([]).colnames)
        dummy_row = np.ones(num_columns)
        table = classification.create_empty_table([], rows=[dummy_row])
        self.assertEqual(1, len(table), 'Table did not add specified rows')


class FitResultsToDict(TestCase):
    """Tests for the ``_fit_results_to_dict`` function"""

    # Don't limit output messages on test failures
    maxDiff = None

    def runTest(self):
        """Run fits for data with a known result and check the results are
        correctly formatted as a table row.
        """

        data = sncosmo.load_example_data()
        params = ['z', 't0', 'x0', 'x1', 'c']
        param_values = [data.meta[p] for p in params]

        # Define mock fit results for a model fit for the last four parameters
        # I.e. for a fit where z is specified
        result = Result({
            'param_names': params,
            'vparam_names': params[1:],
            'parameters': param_values,
            'errors': {p: .1 * v for p, v in data.meta.items()}
        })

        model = sncosmo.Model('salt2')
        model.update(data.meta)
        data.meta['obj_id'] = 'dummy_id'

        row = classification._fit_results_to_dict(
            data, 'dummy_id', 'dummy_band_set', result, model)

        expected_row = {
            'obj_id': 'dummy_id',
            'band': 'dummy_band_set',
            'source': 'salt2',
            'pre_max': 15,
            'post_max': 25,
            'z': result.parameters[0],
            't0': result.parameters[1],
            'x0': result.parameters[2],
            'x1': result.parameters[3],
            'c': result.parameters[4],
            'z_err': result.errors['z'],
            't0_err': result.errors['t0'],
            'x0_err': result.errors['x0'],
            'x1_err': result.errors['x1'],
            'c_err': result.errors['c'],
            'chisq': 36.44,
            'ndof': len(data) - len(result.vparam_names),
            'b_max': -19.5,
            'delta_15': 0.953,
            'message': 'NONE'
        }

        self.assertEqual(expected_row, row)


class RaiseUnspecifiedParams(TestCase):
    """Tests for ``_raise_unspecified_params``"""

    def test_unspecified_param(self):
        """Test a RuntimeError is raise for an unspecified parameter"""

        fixed_params = ['z']
        prior = {'t0': 1, 'x1': 1}
        args = (fixed_params, prior)
        func = classification._raise_unspecified_params
        self.assertRaises(RuntimeError, func, *args)

    def test_all_params_specified(self):
        """Test no error is raised when all fixed params are specified"""

        fixed_params = ['z']
        prior = {'z': .5, 't0': 1, 'x1': 1}
        classification._raise_unspecified_params(fixed_params, prior)


class ClassificationCoords(TestCase):
    """Tests for the ``classify_targets`` function"""

    expected_input_columns = ['obj_id', 'source', 'band', 'chisq', 'ndof']

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
