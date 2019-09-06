#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``fit_funcs`` module."""

from copy import deepcopy
from unittest import TestCase

import sncosmo

from analysis_pipeline import fit_funcs


class TestMutation(TestCase):
    """Test arguments are not mutated by fit functions"""

    def _generic_test(self, func):
        # Use sncosmo example data for testing
        test_data = sncosmo.load_example_data()
        test_model = sncosmo.Model('salt2')
        test_params = test_model.param_names
        test_bounds = {
            'z': (0.3, 0.7),
            't0': (test_data.meta['t0'] - 1, test_data.meta['t0'] + 1),
            'x0': (test_data.meta['x0'] - 1, test_data.meta['x0'] + 1),
            'x1': (test_data.meta['x1'] - 1, test_data.meta['x1'] + 1),
            'c': (test_data.meta['c'] - 1, test_data.meta['c'] + 1)}

        # Preserve original input data
        original_data = deepcopy(test_data)
        original_model = deepcopy(test_model)
        original_bounds = deepcopy(test_bounds)
        original_params = deepcopy(test_params)

        # Check for argument mutation
        func(test_data,
             test_model,
             vparam_names=test_model.param_names,
             bounds=test_bounds)

        self.assertTrue(
            all(original_data == test_data),
            '`data` argument was mutated')

        self.assertSequenceEqual(
            original_params, test_params,
            '`vparam_names` argument was mutated')

        self.assertEqual(
            original_bounds, test_bounds,
            '`bounds` argument was mutated')

        self.assertSequenceEqual(
            original_model.parameters.tolist(),
            test_model.parameters.tolist(),
            '`model` argument was mutated')

    def test_simple_fit(self):
        """Test fit_funcs.simple_fit"""

        self._generic_test(fit_funcs.simple_fit)

    def test_nest_fit(self):
        """Test fit_funcs.nest_fit"""

        self._generic_test(fit_funcs.nest_fit)

    def test_mcmc_fit(self):
        """Test fit_funcs.mcmc_fit"""

        self._generic_test(fit_funcs.mcmc_fit)

    def test_nested_simple_fit(self):
        """Test fit_funcs.nested_simple_fit"""

        self._generic_test(fit_funcs.nested_simple_fit)

    def test_nested_mcmc_fit(self):
        """Test fit_funcs.nested_mcmc_fit"""

        self._generic_test(fit_funcs.nested_mcmc_fit)
