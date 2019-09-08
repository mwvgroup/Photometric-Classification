#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``fit_funcs`` module."""

from copy import deepcopy
from unittest import TestCase

import sncosmo

from analysis_pipeline import fit_funcs


class TestMutation(TestCase):
    """Test arguments are not mutated by fit functions"""

    def _test_fit_function(self, func):
        """Test a given fitting function does not mutate arguments

        Args:
            func (Callable): A sncosmo style fitting function
        """

        # Use sncosmo example data for testing
        data = sncosmo.load_example_data()
        model = sncosmo.Model('salt2')
        params = model.param_names
        bounds = {
            'z': (0.3, 0.7),
            't0': (data.meta['t0'] - 1, data.meta['t0'] + 1),
            'x0': (data.meta['x0'] - 1, data.meta['x0'] + 1),
            'x1': (data.meta['x1'] - 1, data.meta['x1'] + 1),
            'c': (data.meta['c'] - 1, data.meta['c'] + 1)}

        # Preserve original input data
        original_data = deepcopy(data)
        original_model = deepcopy(model)
        original_bounds = deepcopy(bounds)
        original_params = deepcopy(params)

        # Check for argument mutation
        func(data, model, vparam_names=model.param_names, bounds=bounds)

        self.assertTrue(
            all(original_data == data),
            '`data` argument was mutated')

        self.assertSequenceEqual(
            original_params, params,
            '`vparam_names` argument was mutated')

        self.assertEqual(
            original_bounds, bounds,
            '`bounds` argument was mutated')

        self.assertSequenceEqual(
            original_model.parameters.tolist(),
            model.parameters.tolist(),
            '`model` argument was mutated')

    def test_simple_fit(self):
        """Test fit_funcs.simple_fit"""

        self._test_fit_function(fit_funcs.simple_fit)

    def test_nest_fit(self):
        """Test fit_funcs.nest_fit"""

        self._test_fit_function(fit_funcs.nest_fit)

    def test_mcmc_fit(self):
        """Test fit_funcs.mcmc_fit"""

        self._test_fit_function(fit_funcs.mcmc_fit)

    def test_nested_simple_fit(self):
        """Test fit_funcs.nested_simple_fit"""

        self._test_fit_function(fit_funcs.nested_simple_fit)

    def test_nested_mcmc_fit(self):
        """Test fit_funcs.nested_mcmc_fit"""

        self._test_fit_function(fit_funcs.nested_mcmc_fit)
