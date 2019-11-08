#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``simulation.spectra`` module."""

from unittest import TestCase

import numpy as np

from phot_class import spectra


class Area(TestCase):
    """Tests for the ``feature_area`` function"""

    def test_tophat(self):
        """Test the correct area is returned for an inverse top-hat feature"""

        wave = np.arange(1000, 3000)
        flux = np.ones_like(wave)
        flux[100:-100] = 0

        # We use a simulated flux that will remain unchanged when normalized
        # This means the feature area is the same as the width of the feature
        expected_area = len(wave) - 200
        returned_area = spectra.feature_area(wave, flux)
        self.assertEqual(expected_area, returned_area)

    def test_no_feature(self):
        """Test zero is returned for a spectrum without a feature
        (i.e. for y=x)
        """

        wave = np.arange(1000, 3000)
        self.assertEqual(0, spectra.feature_area(wave, wave))


class PEW(TestCase):
    """Tests for the ``feature_pew`` function"""

    def test_tophat(self):
        """Test the correct pew is returned for an inverse top-hat"""

        wave = np.arange(1000, 3000)
        flux = np.ones_like(wave)
        flux[100:-100] = 0

        # We use a simulated flux that will remain unchanged when normalized
        # This means the feature pew is the same as the width of the feature
        expected_area = len(wave) - 200
        normed_flux, returned_area = spectra.feature_pew(wave, flux)

        self.assertEqual(expected_area, returned_area)

    def test_no_feature(self):
        """Pass a dummy spectra that is a straight line (f = 2 * lambda)
        and check that the pew is zero.
        """

        wave = np.arange(1000, 3000)
        flux = 2 * wave

        self.assertEqual(0, spectra.feature_pew(wave, flux)[1])

    def test_normalization(self):
        """Pass a dummy spectra that is a straight line (f = 2 * lambda)
        and check that the normalized flux is an array of ones.
        """

        wave = np.arange(1000, 3000)
        flux = 2 * wave

        norm_flux, _ = spectra.feature_pew(wave, flux)
        expected_norm_flux = np.ones_like(flux).tolist()

        self.assertListEqual(expected_norm_flux, norm_flux.tolist())


class Velocity(TestCase):
    # Todo: Write tests for ``feature_velocity``

    def runTest(self):
        gaussian = lambda x, amplitude, avg, stddev, offset: \
            amplitude * np.exp(-((x - avg) ** 2) / (2 * stddev ** 2)) + offset

        amp = 10
        avg = 2000
        stddev = 15
        offset = 500
        wave = np.arange(1000, 3000)
        flux = gaussian(wave, amp, avg, stddev, offset)
        eflux = np.zeros_like(flux)

        vel = spectra.feature_velocity(wave, flux, eflux, avg - 100)
        print(vel)
        self.fail()


class Properties(TestCase):
    # Todo: Write tests for ``calc_feature_properties``

    def runTest(self):
        self.fail()
