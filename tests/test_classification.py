#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``classification`` module."""

from copy import deepcopy
from unittest import TestCase

import numpy as np
import sncosmo
from astropy.table import Table
from sncosmo.utils import Result

from phot_class import classification, models


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


class PlotLc(TestCase):
    """Tests for the ``_plot_lc`` function"""

    def runTest(self):
        """Test the returned figure is not empty"""

        data = sncosmo.load_example_data()
        model = sncosmo.Model(source='salt2')
        result, fitted_model = sncosmo.fit_lc(
            data, model,
            ['z', 't0', 'x0', 'x1', 'c'],
            bounds={'z': (0.3, 0.7)})
        fig = classification._plot_lc(data, result, fitted_model, show=False)
        self.assertTrue(fig.get_axes(), 'Returned figure has no axes')


class BandFits(TestCase):
    """Tests for the ``run_band_fits`` function"""

    @classmethod
    def setUpClass(cls):
        """Define default arguments for running a successful set of fits"""

        models.register_sources(force=True)
        cls.data = sncosmo.load_example_data()
        cls.default_args = dict(
            obj_id='dummy_id',
            data=cls.data,
            vparams=['amplitude', 'x1', 'c'],
            fit_func=sncosmo.fit_lc,
            priors_hs={'z': cls.data.meta['z'], 't0': cls.data.meta['t0']},
            priors_bg={'z': cls.data.meta['z'], 't0': cls.data.meta['t0']},
            kwargs_hs={'bounds': {'x1': (-1, 1)}},
            kwargs_bg={'bounds': {'x1': (0.65, 1.25), 'c': (0, 1)}},
        )

        cls.returned = classification.run_band_fits(**cls.default_args)

    def test_returned_bands(self):
        """Test the returned table has two fits per band"""

        # Remember we fit each band with two models, but also fit all bands
        expected_bands = 2 * (list(set(self.data['band'])) + ['all'])
        self.assertCountEqual(expected_bands, self.returned['band'])

    def test_unconstrained_param(self):
        """Test an error is raised for missing values in priors"""

        # Missing param for hsiao_x1 model
        missing_hs_t0 = deepcopy(self.default_args)
        del missing_hs_t0['priors_hs']['t0']
        self.assertRaises(
            RuntimeError, classification.run_band_fits, **missing_hs_t0)

        # Missing param for sn91bg model
        missing_bg_t0 = deepcopy(self.default_args)
        del missing_bg_t0['priors_bg']['t0']
        self.assertRaises(
            RuntimeError, classification.run_band_fits, **missing_bg_t0)


# Todo test redshift dependence
class ClassificationCoords(TestCase):
    """Tests for the ``classify_targets`` function"""

    expected_input_columns = ['obj_id', 'source', 'band', 'chisq', 'ndof']

    def test_coordinate_calculation(self):
        """Test correct coordinates are returned for the given input data"""

        h_chisq_u, h_dof_u = 10, 1
        h_chisq_g, h_dof_g = 20, 1
        h_chisq_i, h_dof_i = 30, 2
        h_chisq_z, h_dof_z = 40, 2

        b_chisq_u, b_dof_u = 1, 1
        b_chisq_g, b_dof_g = 2, 1
        b_chisq_i, b_dof_i = 3, 2
        b_chisq_z, b_dof_z = 4, 2

        expected_x = (
                ((h_chisq_u + h_chisq_g) / (h_dof_u + h_dof_g))
                - ((b_chisq_u + b_chisq_g) / (b_dof_u + b_dof_g))
        )

        expected_y = (
                ((h_chisq_i + h_chisq_z) / (h_dof_i + h_dof_z))
                - ((b_chisq_i + b_chisq_z) / (b_dof_i + b_dof_z))
        )

        test_data = Table(names=self.expected_input_columns, rows=[
            ['dummy_id', 'hsiao_x1', 'sdssu', h_chisq_u, h_dof_u],
            ['dummy_id', 'hsiao_x1', 'sdssg', h_chisq_g, h_dof_g],
            ['dummy_id', 'hsiao_x1', 'sdssi', h_chisq_i, h_dof_i],
            ['dummy_id', 'hsiao_x1', 'sdssz', h_chisq_z, h_dof_z],

            ['dummy_id', 'sn91bg', 'sdssu', b_chisq_u, b_dof_u],
            ['dummy_id', 'sn91bg', 'sdssg', b_chisq_g, b_dof_g],
            ['dummy_id', 'sn91bg', 'sdssi', b_chisq_i, b_dof_i],
            ['dummy_id', 'sn91bg', 'sdssz', b_chisq_z, b_dof_z],
        ])

        test_data.meta['band_names'] = ['sdssu', 'sdssg', 'sdssi', 'sdssz']
        test_data.meta['lambda_eff'] = \
            [sncosmo.get_bandpass(b).wave_eff for b in
             test_data.meta['band_names']]

        expected_row = ['dummy_id', expected_x, expected_y]
        class_coordinates = classification.classify_targets(test_data)
        self.assertListEqual(list(class_coordinates[0]), expected_row)
