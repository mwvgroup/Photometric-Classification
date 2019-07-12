#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Test core functionality of light curve fitting from
analysis_pipeline.lc_fitting._fit_funcs
"""

from copy import copy
from copy import deepcopy
from unittest import TestCase

import numpy as np
import sncosmo
from astropy.table import Column, Table
from sndata.csp import dr3

from analysis_pipeline.lc_fitting import calc_chisq
from analysis_pipeline.lc_fitting import create_empty_summary_table
from analysis_pipeline.lc_fitting import fit_lc
from analysis_pipeline.lc_fitting import nest_lc

TEST_ID = '2005kc'  # Use the same CSP object for all applicable tests
dr3.download_module_data()
dr3.register_filters(force=True)


def copy_data(*args):
    """Return a deep copy of passed objects"""

    if len(args) == 1:
        return deepcopy(args[0])

    return (deepcopy(a) for a in args)


class TestSummaryTable(TestCase):
    """Test the fit result summary table is created with the correct columns"""

    def runTest(self):
        test_param_names = ['var1', 'var2']
        expected_names = set(create_empty_summary_table([]).colnames)
        expected_names.update(test_param_names)
        expected_names.update((n + '_err' for n in test_param_names))

        returned_table = create_empty_summary_table(test_param_names)
        returned_names = set(returned_table.colnames)

        self.assertEqual(expected_names, returned_names,
                         'Incorrect column names.')


class TestChisqCalculation(TestCase):
    """Tests for analysis_pipeline.lc_fitting.calc_chisq"""

    @staticmethod
    def create_test_table(model):
        test_table = Table()
        test_table['time'] = np.arange(0, 10, 1)
        test_table['band'] = Column(np.full(len(test_table), 'sdssu'),
                                    dtype='U100')

        modeled_flux = model.bandflux(test_table['band'], test_table['time'])
        test_table['flux'] = 10 * modeled_flux
        test_table['fluxerr'] = .1 * test_table['flux']

        return test_table

    def test_mutable_args(self):
        """Test that input data tables are not mutated"""

        model = sncosmo.Model('salt2')
        test_table_new = self.create_test_table(model)
        test_table_old = copy(test_table_new)

        calc_chisq(test_table_new, model)
        self.assertTrue(all(test_table_new == test_table_old),
                        'Data table was mutated')

    def test_drop_bands(self):
        """Test that out of range data and unknown bands are dropped"""

        model = sncosmo.Model('salt2')
        test_table = self.create_test_table(model)
        expected_chisq = calc_chisq(test_table, model)

        # Add values that are out of the model range
        test_table.add_row([100000000, 'sdssu', 0, 1])  # Out of time range
        test_table.add_row([1, 'csp_dr3_H', 1, 1])  # Out of wavelength range
        actual_chisq = calc_chisq(test_table, model)

        self.assertEqual(expected_chisq[0], actual_chisq[0],
                         'Incorrect chisq value')

        self.assertEqual(expected_chisq[1], actual_chisq[1],
                         'Incorrect number of data points')


class TestLCFitting(TestCase):
    """Test arguments are mutated as expected during nested sampling"""

    @staticmethod
    def create_test_data():
        test_data = dr3.get_data_for_id(TEST_ID, format_sncosmo=True)
        model = sncosmo.Model('salt2')
        bounds = {
            'z': [0.0035, 0.084],
            'x0': [0, 0.05],
            'x1': [-5, 5],
            'c': [-1., 1.]
        }

        return test_data, model, bounds

    def test_nest_lc_mutation(self):
        """Test arguments are mutated as expected during nested sampling"""

        test_data, model, bounds = self.create_test_data()

        # Preserve original input data
        original_data, original_model, original_bounds = \
            copy_data(test_data, model, bounds)

        # Check for argument mutation
        nest_lc(test_data, model,
                vparam_names=model.param_names,
                bounds=bounds)

        self.assertTrue(all(original_data == test_data), 'Data was mutated')
        self.assertEqual(original_bounds, bounds, 'Bounds were mutated')
        self.assertSequenceEqual(
            original_model.parameters.tolist(),
            model.parameters.tolist(),
            'Model was mutated')

    def test_fit_lc_mutation(self):
        test_data, model, bounds = self.create_test_data()

        # Preserve original input data
        original_data, original_model, original_bounds = \
            copy_data(test_data, model, bounds)

        # Check for argument mutation
        nest_lc(test_data, model,
                vparam_names=model.param_names,
                bounds=bounds)

        self.assertTrue(all(original_data == test_data), 'Data was mutated')
        self.assertEqual(original_bounds, bounds, 'Bounds were mutated')
        self.assertSequenceEqual(
            original_model.parameters.tolist(),
            model.parameters.tolist(),
            'Model was mutated')

    def test_fit_5_params(self):
        """
        Test that fit results agree between the analysis pipeline and
        directly calling SNCosmo
        """

        test_data, model, bounds = self.create_test_data()
        fitting_kwargs = dict(
            data=test_data,
            model=model,
            warn=False,
            vparam_names=['z', 't0', 'x0', 'x1', 'c'],
            bounds=bounds)

        # Run fits with analysis pipeline and SNCosmo seperatly
        pipeline_result = fit_lc(**deepcopy(fitting_kwargs))
        sncosmo_result, _ = sncosmo.fit_lc(**deepcopy(fitting_kwargs))

        for i, value in enumerate(sncosmo_result.parameters):
            param = fitting_kwargs['vparam_names'][i]
            err_msg = f'Parameter {param} disagrees'
            self.assertEqual(pipeline_result[i + 2], value, err_msg)
