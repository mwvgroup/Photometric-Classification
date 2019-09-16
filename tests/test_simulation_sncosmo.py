#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``simulation.sncosmo_sims`` module."""

from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase

import numpy as np
import sncosmo
from astropy.table import Table

from phot_class import models
from phot_class.simulation import sncosmo_sims

models.register_sources(force=True)


class StretchColorSimulation(TestCase):
    """Tests for the bg_stretch_color function"""

    @classmethod
    def setUpClass(cls):
        """Simulate a set of parameters to test against"""

        np.random.seed(1)  # Seed so we have consistent test results
        cls.sim_stretch, cls.sim_color = sncosmo_sims.bg_stretch_color(1e5)

    def test_average_stretch(self):
        """Test simulated params have correct average stretch +/- .001"""

        average_stretch = np.average(self.sim_stretch)
        is_correct_stretch = np.isclose(
            sncosmo_sims.AVG_STRETCH,
            average_stretch,
            rtol=0, atol=.001)

        err_msg = f'Expected {sncosmo_sims.AVG_STRETCH}, found {average_stretch}'
        self.assertTrue(is_correct_stretch, err_msg)

    def test_average_color(self):
        """Test simulated params have correct average color +/- .001"""

        average_color = np.average(self.sim_color)
        is_correct_stretch = np.isclose(
            sncosmo_sims.AVG_COLOR,
            average_color,
            rtol=0, atol=.001)

        err_msg = f'Expected {sncosmo_sims.AVG_COLOR}, found {average_color}'
        self.assertTrue(is_correct_stretch, err_msg)

    def test_covariance(self):
        """Test simulated params have correct covariance matrix

        Comparison requires an relative tolerance of 10%
        """

        covariance = np.cov([self.sim_stretch, self.sim_color])
        is_correct_cov = np.isclose(
            sncosmo_sims.COVARIANCE, covariance, atol=0, rtol=.1).all()

        err_msg = f'Expected {sncosmo_sims.COVARIANCE}, found {covariance}'
        self.assertTrue(is_correct_cov, err_msg)

    def test_simulation_bounds(self):
        """Test simulated parameters are within specified boundaries"""

        # We use arbitrarily chosen boundaries
        min_stretch = 0.5
        max_stretch = .8
        min_color = .3
        max_color = .7

        np.random.seed(1)
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


class LCParameterSimulation(TestCase):
    """Tests for the ``sim_bg_params`` function."""

    @classmethod
    def setUpClass(cls):
        cls.zmin = 0
        cls.zmax = 1
        cls.tmin = 100
        cls.tmax = 1000

        np.random.seed(1)  # Seed so we have consistent test results
        parameter_dicts = sncosmo_sims.sim_bg_params(
            zmin=cls.zmin, zmax=cls.zmax, tmin=cls.tmin, tmax=cls.tmax)

        cls.parameters = Table(rows=parameter_dicts)

    def test_redshift_range(self):
        """Test simulated redshifts are all within the specified range"""

        is_in_range = np.logical_and(
            self.parameters['z'] < self.zmax,
            self.parameters['z'] > self.zmin).all()

        self.assertTrue(is_in_range)

    def test_time_range(self):
        """Test simulated times are all within the specified range"""

        is_in_range = np.logical_and(
            self.parameters['t0'] < self.tmax,
            self.parameters['t0'] > self.tmin).all()

        self.assertTrue(is_in_range)


class GenerateLC(TestCase):
    """Tests for the ``generate_lc`` function."""

    @classmethod
    def setUpClass(cls):
        cls.temp_dir = TemporaryDirectory()
        cls.model = sncosmo.Model('sn91bg')
        cls.phase_range = (-10, 40)
        cls.model_params = sncosmo_sims.sim_bg_params(
            zmin=0, zmax=1, tmin=100, tmax=150
        )

        sncosmo_sims.generate_lc(
            cls.model, cls.phase_range, cls.model_params, cls.temp_dir.name)

    @classmethod
    def tearDownClass(cls):
        cls.temp_dir.cleanup()

    def test_number_simulations(self):
        """Test one light-curve file exists for every set of parameters"""

        files = list(Path(self.temp_dir.name).glob('*.ecsv'))
        self.assertGreaterEqual(len(files), 0)
        self.assertEqual(len(self.model_params), len(files))
