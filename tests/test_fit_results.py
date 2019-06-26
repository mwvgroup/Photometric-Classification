#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Test that fit results generated iterably match results obtained by an
independent run.
"""

from copy import deepcopy
from unittest import TestCase

import sncosmo
from SNData.csp import dr3

from analysis_pipeline import fit_lc

dr3.download_module_data()
dr3.register_filters(force=True)


class TestSNCosmoAgreement(TestCase):
    """Test that fit results agree between the analysis pipeline and
    directly calling SNCosmo"""

    @classmethod
    def setUpClass(cls):
        cls.test_cid = '2005kc'
        cls.input_table = dr3.get_sncosmo_input(cls.test_cid)

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
