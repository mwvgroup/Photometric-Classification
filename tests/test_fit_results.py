#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Test that fit results generated iterably match results obtained by an
independent run.
"""

import os
from copy import deepcopy
from unittest import TestCase, skipIf

import sncosmo
from astropy.table import Table

from analysis_pipeline.data_access import csp, des, sdss
from analysis_pipeline.lc_fitting import fit_lc

pipeline_output_dir = '../pipeline_outputs'
csp_dir = os.path.join(pipeline_output_dir, 'csp')
des_dir = os.path.join(pipeline_output_dir, 'des')
sdss_dir = os.path.join(pipeline_output_dir, 'sdss')
salt_2_4 = sncosmo.Model(source=sncosmo.get_source('salt2', version='2.4'))


class Test4Param(TestCase):
    """Run independent fits and compare with tabulated results"""

    def check_pipeline_entry(self, pipeline_entry, input_table, fit_kwargs):
        """Test fit results for an arbitrary survey"""

        test_model = deepcopy(salt_2_4)
        test_model.set(z=pipeline_entry['z'])
        test_results = fit_lc(
            input_table,
            test_model,
            ['t0', 'x0', 'x1', 'c'],
            **fit_kwargs)

        for i, key in enumerate(('z', 't0', 'x0', 'x1', 'c')):
            self.assertEqual(test_results[i], pipeline_entry[key])

    @skipIf(not os.path.exists(csp_dir), 'No CSP pipeline outputs found.')
    def test_csp(self):
        """Test fit results for CSP"""

        path = os.path.join(csp_dir, 'salt2_2.4_4param_all.ecsv')
        pipeline_results = Table.read(path)
        fit_kwargs = pipeline_results.meta

        for i in (0, 1, 5):
            test_cid = pipeline_results[i]['cid']
            test_input_table = csp.get_input_for_id(test_cid)
            self.check_pipeline_entry(
                pipeline_results[i], test_input_table, fit_kwargs)

    @skipIf(not os.path.exists(des_dir), 'No DES pipeline outputs found.')
    def test_des(self):
        """Test fit results for DES"""

        path = os.path.join(des_dir, 'salt2_2.4_4param_all.ecsv')
        pipeline_results = Table.read(path)
        fit_kwargs = pipeline_results.meta

        for i in (0, 1, 5):
            test_cid = pipeline_results[i]['cid']
            test_input_table = des.get_input_for_id(test_cid)
            self.check_pipeline_entry(
                pipeline_results[i], test_input_table, fit_kwargs)

    @skipIf(not os.path.exists(sdss_dir), 'No SDSS pipeline outputs found.')
    def test_sdss(self):
        """Test fit results for SDSS"""

        path = os.path.join(sdss_dir, 'salt2_2.4_4param_all.ecsv')
        pipeline_results = Table.read(path)
        fit_kwargs = pipeline_results.meta

        for i in (0, 1, 5):
            test_cid = pipeline_results[i]['cid']
            test_input_table = sdss.get_input_for_id(test_cid)
            self.check_pipeline_entry(
                pipeline_results[i], test_input_table, fit_kwargs)


class TestSNCosmoAgreement(TestCase):

    def runTest(self):
        test_cid = '2005kc'
        input_table = csp.get_input_for_id(test_cid)

        fitting_kwargs = dict(modelcov=True, warn=False,
                              bounds={'z': [0.002, 0.085]})

        pipeline_result = fit_lc(
            data=input_table,
            model=salt_2_4,
            vparam_names=['z', 't0', 'x0', 'x1', 'c'],
            **fitting_kwargs)

        sncosmo_result, _ = sncosmo.fit_lc(
            data=input_table,
            model=salt_2_4,
            vparam_names=['z', 't0', 'x0', 'x1', 'c'],
            **fitting_kwargs)

        params = ['z', 't0', 'x0', 'x1', 'c']
        for i, value in enumerate(sncosmo_result.parameters):
            err_msg = f'Parameter {params[i]} disagrees'
            self.assertEqual(pipeline_result[i + 2],
                             value,
                             err_msg)
