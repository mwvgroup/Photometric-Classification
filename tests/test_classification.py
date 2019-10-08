#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``classification`` module."""

from unittest import TestCase

import sncosmo
from astropy.table import Table

from phot_class import classification


# Todo test redshift dependence
class ClassificationCoords(TestCase):
    """Tests for the ``classify_targets`` function"""

    @classmethod
    def setUpClass(cls):
        """Define some moke data for testing"""

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

        cls.test_data = Table(
            names=['obj_id', 'source', 'band', 'chisq', 'ndof', 'z',
                   'message'],
            dtype=['U10', 'U10', 'U10', float, float, float, 'U10'],
            rows=[
                ['dummy_id', 'hsiao_x1', 'all', 0, 0, 0, ''],
                ['dummy_id', 'hsiao_x1', 'sdssu', h_chisq_u, h_dof_u, 0, ''],
                ['dummy_id', 'hsiao_x1', 'sdssg', h_chisq_g, h_dof_g, 0, ''],
                ['dummy_id', 'hsiao_x1', 'sdssi', h_chisq_i, h_dof_i, 0, ''],
                ['dummy_id', 'hsiao_x1', 'sdssz', h_chisq_z, h_dof_z, 0, ''],

                ['dummy_id', 'sn91bg', 'sdssu', b_chisq_u, b_dof_u, 0, ''],
                ['dummy_id', 'sn91bg', 'sdssg', b_chisq_g, b_dof_g, 0, ''],
                ['dummy_id', 'sn91bg', 'sdssi', b_chisq_i, b_dof_i, 0, ''],
                ['dummy_id', 'sn91bg', 'sdssz', b_chisq_z, b_dof_z, 0, ''],
            ]
        )

        cls.test_data.meta['band_names'] = ['sdssu', 'sdssg', 'sdssi', 'sdssz']
        cls.test_data.meta['lambda_eff'] = \
            [sncosmo.get_bandpass(b).wave_eff for b in
             cls.test_data.meta['band_names']]

        cls.expected_row = ['dummy_id', expected_x, expected_y]

    def test_coordinate_calculation(self):
        """Test correct coordinates are returned for the given input data"""

        class_coordinates = classification.classify_targets(self.test_data)
        self.assertListEqual(list(class_coordinates[0]), self.expected_row)

    def test_failed_fits(self):
        """Test failed fits are ignored"""

        test_data = self.test_data.copy()
        test_data[0]['message'] = 'Failed'
        class_coordinates = classification.classify_targets(test_data)
        self.assertEqual(0, len(class_coordinates))
