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

        # We use a simulated flux that will remain unchanged when normalized
        # This means the feature area is the same as the width of the feature
        wave = np.arange(1000, 3000)
        flux = np.ones_like(wave)
        flux[100:-100] = 0

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


class FeatureIdentification(TestCase):
    """Test the identification of feature boundaries"""

    @classmethod
    def setUpClass(cls):
        # Create dummy spectrum
        wave_range = (7000, 8001, 100)
        cls.wavelength = np.arange(*wave_range)
        cls.flux = np.ones_like(cls.wavelength)  # Define continuum

        cls.peak_wavelengths = (7100, 7500)
        for peak in cls.peak_wavelengths:
            cls.flux[cls.wavelength == peak] = 10  # Add delta function peak

    def test_peak_coordinates(self):
        """Test the correct peak wavelength is found for a single flux spike"""

        expected_peak = self.peak_wavelengths[0]
        recovered_peak = spectra.get_peak_wavelength(
            self.wavelength,
            self.flux,
            expected_peak - 10,
            expected_peak + 10
        )

        self.assertEqual(expected_peak, recovered_peak)

    def test_unobserved_feature(self):
        """Test an error is raise if the feature is out of bounds"""

        max_wavelength = max(self.wavelength)
        with self.assertRaises(ValueError):
            spectra.get_peak_wavelength(
                self.wavelength,
                self.flux,
                max_wavelength + 10,
                max_wavelength + 20
            )

    def test_double_peak(self):
        """Test the correct feature wavelengths are found"""

        lower_peak_wavelength = min(self.peak_wavelengths)
        upper_peak_wavelength = max(self.peak_wavelengths)
        recovered_lower_peak = spectra.get_peak_wavelength(
            self.wavelength,
            self.flux,
            lower_peak_wavelength - 10,
            upper_peak_wavelength + 10,
            'min'
        )

        recovered_upper_peak = spectra.get_peak_wavelength(
            self.wavelength,
            self.flux,
            lower_peak_wavelength - 10,
            upper_peak_wavelength + 10,
            'max'
        )

        self.assertEqual(
            lower_peak_wavelength, recovered_lower_peak, 'Incorrect min peak')

        self.assertEqual(
            upper_peak_wavelength, recovered_upper_peak, 'Incorrect max peak')

    def test_feature_bounds(self):
        lower_peak_wavelength = min(self.peak_wavelengths)
        upper_peak_wavelength = max(self.peak_wavelengths)
        feature_dict = {
            'lower_blue': lower_peak_wavelength - 10,
            'upper_blue': lower_peak_wavelength + 10,
            'lower_red': upper_peak_wavelength - 10,
            'upper_red': upper_peak_wavelength + 10
        }

        feat_start, feat_end = spectra.get_feature_bounds(
            self.wavelength, self.flux, feature_dict)

        self.assertEqual(
            lower_peak_wavelength, feat_start, 'Incorrect min peak')

        self.assertEqual(
            upper_peak_wavelength, feat_end, 'Incorrect max peak')
