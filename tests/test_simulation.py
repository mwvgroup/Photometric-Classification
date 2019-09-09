#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``models`` module."""

from unittest import TestCase

import numpy as np

from analysis_pipeline import models
from analysis_pipeline.simulation import sncosmo_sims

models.register_sources(force=True)


# Todo: Test ``sim_bg_params`` and ``generate_lc``

class StretchColorSimulation(TestCase):

    @classmethod
    def setUpClass(cls):
        """Simulate a set of parameters to test against"""

        np.random.seed(12345)  # So we have consistent test results
        cls.sim_stretch, cls.sim_color = sncosmo_sims.bg_stretch_color(1e5)

        # Since we are dealing with a sampling from a distribution, we define
        # an absolute tolerance for determining if returned values match
        # expected values
        cls.tolerance = .0005

    def test_average_stretch(self):
        """Test simulated params have correct average stretch"""

        average_stretch = np.average(self.sim_stretch)
        is_correct_stretch = np.isclose(
            sncosmo_sims.AVG_STRETCH,
            average_stretch,
            rtol=0, atol=self.tolerance)

        err_msg = f'Expected {sncosmo_sims.AVG_STRETCH}, found {average_stretch}'
        self.assertTrue(is_correct_stretch, err_msg)

    def test_average_color(self):
        """Test simulated params have correct average color"""

        average_color = np.average(self.sim_color)
        is_correct_stretch = np.isclose(
            sncosmo_sims.AVG_COLOR,
            average_color,
            rtol=0, atol=self.tolerance)

        err_msg = f'Expected {sncosmo_sims.AVG_COLOR}, found {average_color}'
        self.assertTrue(is_correct_stretch, err_msg)

    def test_covariance(self):
        """Test simulated params have correct covariance matrix"""

        covariance = np.cov([self.sim_stretch, self.sim_color])
        is_correct_cov = np.isclose(sncosmo_sims.COVARIANCE, covariance).all()

        err_msg = f'Expected {sncosmo_sims.COVARIANCE}, found {covariance}'
        self.assertTrue(is_correct_cov, err_msg)

    def test_sim_bounds(self):
        """test simulated parameters are within specified boundaries"""

        # We use arbitrarily chosen boundaries
        min_stretch = 0.5
        max_stretch = .8
        min_color = .3
        max_color = .7

        np.random.seed(12345)
        stretch, color = sncosmo_sims.bg_stretch_color(
            size=1000,
            min_stretch=min_stretch,
            max_stretch=max_stretch,
            min_color=min_color,
            max_color=max_color
        )

        self.assertLessEqual(
            min_stretch, min(stretch), 'Minimum stretch out of bounds')

        self.assertGreaterEqual(
            max_stretch, max(stretch), 'Maximum stretch out of bounds')

        self.assertLessEqual(
            min_color, min(color), 'Minimum stretch out of bounds')

        self.assertGreaterEqual(
            max_color, max(color), 'Maximum stretch out of bounds')
