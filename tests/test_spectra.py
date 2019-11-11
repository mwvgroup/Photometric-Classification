#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``simulation.spectra`` module."""

from unittest import TestCase

import numpy as np
from uncertainties.unumpy import uarray

from phot_class import spectra


class SimulatedSpectrum:
    """Functions for simulating dummy spectra"""

    @staticmethod
    def tophat(wave, m, b, start, end, height, seed=0):
        """Simulate a top-hat absorption feature with normal errors

        Setting ``height=None`` will simulate just the continuum

        Args:
            wave (ndarray): Array of wavelengths to simulate flux for
            m      (float): Slope of the continuum
            b      (float): Y-intercept of the continuum
            start    (int): Starting index for the top-hat
            end      (int): Ending index for the top-hat
            height (float): Height of the top-hat
            seed   (float): Seed for random number generator (default: 0)

        Returns:
            - An array of flux values
            - An array of error values
        """

        flux = m * wave + b
        if height is not None:
            flux[start: end] = height

        np.random.seed(seed)
        return flux, np.random.random(flux.size)

    @staticmethod
    def gaussian(wave, amplitude=-1, mean=None, stddev=0, offset=0, seed=0):
        """Simulate gaussian flux with normal errors

        Args:
            wave    (ndarray): Array of wavelengths to simulate flux for
            amplitude (float): Amplitude of the Gaussian (default: -1)
            mean      (float): Average of the Gaussian (default: mean of wave)
            stddev    (float): Standard deviation of the Gaussian (default: 0)
            offset    (float): Vertical offset of the Gaussian (default: 0)
            seed      (float): Seed for random number generator (default: 0)

        Returns:
            - An array of flux values
            - An array of error values
        """

        mean = np.mean(wave) if mean is None else mean
        flux = amplitude * np.exp(
            -((wave - mean) ** 2) / (2 * stddev ** 2)
        ) + offset

        np.random.seed(seed)
        return flux, np.random.random(flux.size)


class Area(TestCase):
    """Tests for the ``feature_area`` function"""

    def test_tophat_area(self):
        """Test the correct area is returned for an inverse top-hat feature"""

        # We use a simulated flux that will remain unchanged when normalized
        # This means the feature area is the same as the width of the feature
        wave = np.arange(1000, 3000)
        flux, eflux = SimulatedSpectrum.tophat(wave, 0, 1, 100, -100, 0)

        expected_area = len(wave) - 200
        returned_area = spectra.feature_area(wave, flux)
        self.assertEqual(expected_area, returned_area)

    def test_no_feature(self):
        """Test zero is returned for a spectrum without a feature
        (i.e. for y=x)
        """

        wave = np.arange(1000, 3000)
        self.assertEqual(0, spectra.feature_area(wave, wave))

    def test_uarray_support(self):
        """Test the function supports input arrays with ufloat objects"""

        wave = np.arange(1000, 2000)
        uflux = uarray(*SimulatedSpectrum.gaussian(wave))
        returned_area = spectra.feature_area(wave, uflux)
        self.assertLess(0, returned_area.std_dev)


class PEW(TestCase):
    """Tests for the ``feature_pew`` function"""

    def test_tophat(self):
        """Test the correct pew is returned for an inverse top-hat"""

        wave = np.arange(1000, 3000)
        flux, eflux = SimulatedSpectrum.tophat(wave, 0, 1, 100, -100, 0)

        expected_area = len(wave) - 200
        normed_flux, returned_area = spectra.feature_pew(wave, flux)

        self.assertEqual(expected_area, returned_area)

    def test_no_feature(self):
        """Pass a dummy spectra that is a straight line (f = 2 * lambda)
        and check that the pew is zero.
        """

        wave = np.arange(1000, 3000)
        flux = 2 * wave

        norm_flux, pew = spectra.feature_pew(wave, flux)
        self.assertEqual(0, pew)

    def test_normalization(self):
        """Pass a dummy spectra that is a straight line (f = 2 * lambda)
        and check that the normalized flux is an array of ones.
        """

        wave = np.arange(1000, 3000)
        flux = 2 * wave

        norm_flux, pew = spectra.feature_pew(wave, flux)
        expected_norm_flux = np.ones_like(flux).tolist()

        self.assertListEqual(expected_norm_flux, norm_flux.tolist())

    def test_uarray_support(self):
        """Test the function supports input arrays with ufloat objects"""

        # Simulate gaussian flux with normal errors
        wave = np.arange(1000, 2000)
        flux = - np.exp(-((wave - np.mean(wave)) ** 2) / 2)
        np.random.seed(0)
        eflux = np.random.random(flux.size)

        uflux = uarray(flux, eflux)
        norm_flux, pew = spectra.feature_pew(wave, uflux)
        print(norm_flux, pew)


class Velocity(TestCase):
    """Tests for the ``feature_area`` function"""

    def runTest(self):
        self.fail()


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
        recovered_peak = spectra._calc_properties._get_peak_wavelength(
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
            spectra._calc_properties._get_peak_wavelength(
                self.wavelength,
                self.flux,
                max_wavelength + 10,
                max_wavelength + 20
            )

    def test_double_peak(self):
        """Test the correct feature wavelengths are found"""

        lower_peak_wavelength = min(self.peak_wavelengths)
        upper_peak_wavelength = max(self.peak_wavelengths)
        recovered_lower_peak = spectra._calc_properties._get_peak_wavelength(
            self.wavelength,
            self.flux,
            lower_peak_wavelength - 10,
            upper_peak_wavelength + 10,
            'min'
        )

        recovered_upper_peak = spectra._calc_properties._get_peak_wavelength(
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
