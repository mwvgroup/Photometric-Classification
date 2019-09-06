#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``classification`` module."""

from unittest import TestCase

import numpy as np

from analysis_pipeline import classification


# Todo: Test the following functions
# - classify_targets
# - create_empty_table
# - fit_results_to_table_row
# - run_fits
# - tabulate_fit_results


class TestTableCreation(TestCase):
    """Tests for classification.create_empty_table"""

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
