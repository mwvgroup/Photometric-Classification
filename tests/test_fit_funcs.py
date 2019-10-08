#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``fit_funcs`` module."""

from copy import deepcopy
from unittest import TestCase

import numpy as np
import sncosmo

from phot_class import fit_func_wraps


class BaseTestingClass(TestCase):
    """Tests for an arbitrary sncosmo model"""

    def _test_mutation(self):
        """Test a pipeline fitting function does not mutate arguments"""

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
        self.fit_func(data, model, vparam_names=model.param_names,
                      bounds=bounds)

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

    def _test_agrees_with_sncsomo(self, sncosmo_func):
        """Test a pipeline fitting function returns the same results as sncosmo

        Args:
            sncosmo_func (callable): The fit function to compare against
        """

        data = sncosmo.load_example_data()
        vparams = list(data.meta.keys())
        bounds = {p: (.9 * v, 1.1 * v) for p, v in data.meta.items()}
        model = sncosmo.Model('salt2')
        model.update(data.meta)

        results, _ = self.fit_func(data, model, vparams, bounds=bounds)
        sncosmo_results, _ = sncosmo_func(data, model, vparams, bounds=bounds)

        sncosmo_params = list(sncosmo_results.parameters)
        returned_params = list(results.parameters)
        agrees = np.isclose(sncosmo_params, returned_params, rtol=.1).all()
        if not agrees:
            self.assertListEqual(sncosmo_params, returned_params)


class SimpleFit(BaseTestingClass):
    """Tests for the ``simple_fit`` function"""

    @staticmethod
    def fit_func(*a, **kw):
        return fit_func_wraps.simple_fit(*a, **kw)

    def test_mutation(self):
        """Test arguments are not mutated"""

        self._test_mutation()

    def test_agrees_with_sncosmo(self):
        """Test fit results agree with the equivalent sncosmo function"""

        self._test_agrees_with_sncsomo(sncosmo.fit_lc)


class NestFit(BaseTestingClass):
    """Tests for the ``nest_fit`` function"""

    @staticmethod
    def fit_func(*a, **kw):
        return fit_func_wraps.nest_fit(*a, **kw)

    def test_mutation(self):
        """Test arguments are not mutated"""

        self._test_mutation()

    def test_agrees_with_sncosmo(self):
        """Test fit results agree with the equivalent sncosmo function"""

        self._test_agrees_with_sncsomo(sncosmo.nest_lc)


class MCMCFit(BaseTestingClass):
    """Tests for the ``mcmc_fit`` function"""

    @staticmethod
    def fit_func(*a, **kw):
        return fit_func_wraps.mcmc_fit(*a, **kw)

    def test_mutation(self):
        """Test arguments are not mutated"""

        self._test_mutation()

    def test_agrees_with_sncsomo(self):
        """Test fit results agree with the equivalent sncosmo function"""

        self._test_agrees_with_sncsomo(sncosmo.mcmc_lc)


class NestedSimpleFit(BaseTestingClass):
    """Tests for the ``nested_simple_fit`` function"""

    @staticmethod
    def fit_func(*a, **kw):
        return fit_func_wraps.nested_simple_fit(*a, **kw)

    def test_mutation(self):
        """Test arguments are not mutated"""

        self._test_mutation()


class NestedMCMCFit(BaseTestingClass):
    """Tests for the ``nested_mcmc_fit`` function"""

    @staticmethod
    def fit_func(*a, **kw):
        return fit_func_wraps.nested_mcmc_fit(*a, **kw)

    def test_mutation(self):
        """Test arguments are not mutated"""

        self._test_mutation()
