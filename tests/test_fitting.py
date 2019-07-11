#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Test functions used in light curve fitting.
"""

from copy import copy
from copy import deepcopy
from unittest import TestCase

import numpy as np
import sncosmo
from astropy.table import Column, Table
from sndata.csp import dr3

from analysis_pipeline.lc_fitting import calc_chisq
from analysis_pipeline.lc_fitting import fit_lc

dr3.download_module_data()
dr3.register_filters(force=True)


class TestSNCosmoAgreement(TestCase):
    """Test that fit results agree between the analysis pipeline and
    directly calling SNCosmo"""

    @classmethod
    def setUpClass(cls):
        cls.test_cid = '2005kc'
        cls.input_table = dr3.get_data_for_id(cls.test_cid,
                                              format_sncosmo=True)

    def test_5_params(self):
        """Run test for a 5 parameter fit"""

        # Create model
        model_source = sncosmo.get_source('salt2', version='2.4')
        salt_2_4 = sncosmo.Model(source=model_source)
        fitting_kwargs = dict(
            data=self.input_table,
            model=salt_2_4,
            warn=False,
            vparam_names=['z', 't0', 'x0', 'x1', 'c'],
            bounds={'z': [0.002, 0.085]})

        # Run fits with analysis pipeline and SNCosmo seperatly
        pipeline_result = fit_lc(**deepcopy(fitting_kwargs))
        sncosmo_result, _ = sncosmo.fit_lc(**deepcopy(fitting_kwargs))

        for i, value in enumerate(sncosmo_result.parameters):
            param = fitting_kwargs['vparam_names'][i]
            err_msg = f'Parameter {param} disagrees'
            self.assertEqual(pipeline_result[i + 2], value, err_msg)


class TestChisqCalculation(TestCase):
    """Tests for analysis_pipeline.lc_fitting.calc_chisq"""

    @staticmethod
    def create_test_table(model):
        test_table = Table()
        test_table['time'] = np.arange(0, 10, 1)
        test_table['band'] = Column(np.full(len(test_table), 'sdssu'),
                                    dtype='U100')

        modeled_flux = np.array([model.bandflux(b, t) for b, t in
                                 zip(test_table['band'], test_table['time'])])
        test_table['flux'] = 10 * modeled_flux
        test_table['fluxerr'] = .1 * test_table['flux']

        return test_table

    def test_mutable_args(self):
        """Test that input data tables are not mutated"""

        model = sncosmo.Model('salt2')
        test_table_new = self.create_test_table(model)
        test_table_old = copy(test_table_new)

        calc_chisq(test_table_new, model)
        self.assertTrue(all(test_table_new == test_table_old))

    def test_drop_bands(self):
        """Test that out of range data and unknown bands are dropped"""

        model = sncosmo.Model('salt2')
        test_table = self.create_test_table(model)
        expected_chisq = calc_chisq(test_table, model)

        # Add values that are out of the model range
        test_table.add_row([100000000, 'sdssu', 0, 1])  # Out of time range
        test_table.add_row([1, 'csp_dr3_H', 1, 1])  # out of wavelength range
        actual_chisq = calc_chisq(test_table, model)

        self.assertEqual(expected_chisq[0], actual_chisq[0])
        self.assertEqual(expected_chisq[1], actual_chisq[1])
