#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""Test that our custom model for 91bg like supernovae was ported correctly
for use with SNCosmo.
"""

from unittest import TestCase

import sncosmo
from astropy.table import Table

from sn91bg_model import SN91bgSource


class SN91bgModel(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.model = sncosmo.Model(source=SN91bgSource())

    def test_model_runs(self):
        data = Table.read('sn91bg_50003.csv')
        result, fitted_model = sncosmo.fit_lc(
            data, self.model, ['z','t0','amplitude', 'stretch', 'color'], bounds={'z': (0.01, .1)})

        self.assertLessEqual(result.chisq / result.ndof, 1.01)
