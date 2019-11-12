#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``simulation.spectra`` module."""

from unittest import TestCase

import numpy as np
from astropy.constants import c
from uncertainties.unumpy import uarray

from phot_class import spectra


class SimulatedSpectrum:
    """Functions for simulating dummy spectra"""

    @staticmethod
    def tophat(wave, m=0, b=1, start=100, end=-100, height=0, seed=0):
        """Simulate a top-hat absorption feature with normal errors

        Setting ``height=None`` will simulate just the continuum

        Args:
            wave (ndarray): Array of wavelengths to simulate flux for
            m      (float): Slope of the continuum (default: 0)
            b      (float): Y-intercept of the continuum (default: 1)
            start    (int): Starting index for the top-hat (default: 100)
            end      (int): Ending index for the top-hat (default: -100)
            height (float): Height of the top-hat  (default: 0)
            seed   (float): Seed for random number generator (default: 0)

        Returns:
            - An array of flux values
            - An array of error values
        """

        flux = m * wave + b
        if height is not None:
            flux[start: end] = height

        np.random.seed(seed)
        inverse_snr = np.random.randint(1, high=10, size=flux.size) / 1000
        return flux, flux * inverse_snr

    @staticmethod
    def gaussian(wave, amplitude=-1, mean=None, stddev=1, offset=100, seed=0):
        """Simulate gaussian flux with normal errors

        Args:
            wave    (ndarray): Array of wavelengths to simulate flux for
            amplitude (float): Amplitude of the Gaussian (default: -1)
            mean      (float): Average of the Gaussian (default: mean of wave)
            stddev    (float): Standard deviation of the Gaussian (default: 1)
            offset    (float): Vertical offset of the Gaussian (default: 100)
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
        inverse_snr = np.random.randint(1, high=10, size=flux.size) / 1000
        return flux, flux * inverse_snr

    @staticmethod
    def delta_func(wave, m=0, b=0, peak_wave=(), amplitude=1, seed=0):
        """Simulate linear flux with interspersed delta functions

        Args:
            wave    (ndarray): Array of wavelengths to simulate flux for
            m         (float): Slope of the continuum (default: 0)
            b         (float): Y-intercept of the continuum (default: 1)
            peak_wave  (iter): Starting index for the top-hat (default: 100)
            amplitude (float): Height of the delta functions
            seed      (float): Seed for random number generator (default: 0)

        Returns:
            - An array of flux values
            - An array of error values
        """

        flux = m * wave + b
        flux[np.isin(wave, peak_wave)] = amplitude

        np.random.seed(seed)
        inverse_snr = np.random.randint(1, high=10, size=flux.size) / 1000
        return flux, flux * inverse_snr


class Area(TestCase):
    """Tests for the ``feature_area`` function"""

    def test_tophat_area(self):
        """Test the correct area is returned for an inverse top-hat feature"""

        # We use a simulated flux that will remain unchanged when normalized
        # This means the feature area is the same as the width of the feature
        wave = np.arange(1000, 3000)
        flux, eflux = SimulatedSpectrum.tophat(wave)

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
        uflux = uarray(*SimulatedSpectrum.gaussian(wave, stddev=100))
        returned_area = spectra.feature_area(wave, uflux)
        self.assertLess(0, returned_area.std_dev)


class PEW(TestCase):
    """Tests for the ``feature_pew`` function"""

    def test_tophat(self):
        """Test the correct pew is returned for an inverse top-hat"""

        wave = np.arange(1000, 3000)
        flux, eflux = SimulatedSpectrum.tophat(wave)

        expected_area = len(wave) - 200
        normed_flux, returned_area = spectra.feature_pew(wave, flux)

        self.assertEqual(expected_area, returned_area)

    def test_no_feature(self):
        """Pass a dummy spectra that is a straight line (f = 2 * lambda)
        and check that the pew is zero.
        """

        wave = np.arange(1000, 3000)
        flux, eflux = SimulatedSpectrum.tophat(wave, height=None)

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

        wave = np.arange(1000, 2000)
        uflux = uarray(*SimulatedSpectrum.gaussian(wave, stddev=100))
        norm_flux, pew = spectra.feature_pew(wave, uflux)
        self.assertLess(0, pew.std_dev)


class Velocity(TestCase):
    """Tests for the ``feature_area`` function"""

    def test_velocity_estimation(self):
        wave = np.arange(1000, 2000)
        lambda_rest = np.mean(wave)
        lambda_observed = lambda_rest - 100
        flux, eflux = SimulatedSpectrum.gaussian(
            wave, mean=lambda_observed, stddev=100)

        # Doppler equation: λ_observed = λ_source (c − v_source) / c
        # v_expected = (c * (1 - (lambda_observed / lambda_rest))).value
        lambda_ratio = ((lambda_rest - lambda_observed) / lambda_rest) + 1
        v_expected = c.value * (
                (lambda_ratio ** 2 - 1) / (lambda_ratio ** 2 + 1)
        )

        v_returned = spectra.feature_velocity(
            lambda_rest, wave, flux, unit=c.unit)

        self.assertEqual(v_expected, v_returned)

    def test_uarray_support(self):
        """Test the function supports input arrays with ufloat objects"""

        wave = np.arange(1000, 2000)
        lambda_rest = np.mean(wave)

        flux, eflux = SimulatedSpectrum.gaussian(wave, stddev=100)
        uflux = uarray(flux, eflux)

        velocity_no_err = spectra.feature_velocity(lambda_rest, wave, flux)
        velocity_w_err = spectra.feature_velocity(lambda_rest, wave, uflux)
        self.assertAlmostEqual(velocity_no_err, velocity_w_err.nominal_value)


class FindPeakWavelength(TestCase):
    """Tests for the ``find_peak_wavelength`` function"""

    @classmethod
    def setUpClass(cls):
        # Note that we are simulating delta function emission features
        cls.wave = np.arange(100, 300)
        cls.peak_wavelengths = (210, 250)
        cls.flux, cls.eflux = SimulatedSpectrum.delta_func(
            cls.wave, m=1, b=10, amplitude=500, peak_wave=cls.peak_wavelengths)

    def test_peak_coordinates(self):
        """Test the correct peak wavelength is found for a single flux spike"""

        expected_peak = self.peak_wavelengths[0]
        recovered_peak = spectra.find_peak_wavelength(
            wave=self.wave,
            flux=self.flux,
            lower_bound=expected_peak - 10,
            upper_bound=expected_peak + 10
        )

        self.assertEqual(expected_peak, recovered_peak)

    def test_unobserved_feature(self):
        """Test an error is raise if the feature is out of bounds"""

        max_wavelength = max(self.wave)
        with self.assertRaises(ValueError):
            spectra._calc_properties.find_peak_wavelength(
                wave=self.wave,
                flux=self.flux,
                lower_bound=max_wavelength + 10,
                upper_bound=max_wavelength + 20
            )

    def test_double_peak(self):
        """Test the correct feature wavelengths are found corresponding to
        ``behavior = 'min'`` and ``behavior = 'max'``
        """

        lower_peak_wavelength = min(self.peak_wavelengths)
        upper_peak_wavelength = max(self.peak_wavelengths)
        returned_lower_peak = spectra.find_peak_wavelength(
            self.wave,
            self.flux,
            lower_peak_wavelength - 10,
            upper_peak_wavelength + 10,
            'min'
        )

        self.assertEqual(
            lower_peak_wavelength, returned_lower_peak, 'Incorrect min peak')

        returned_upper_peak = spectra.find_peak_wavelength(
            self.wave,
            self.flux,
            lower_peak_wavelength - 10,
            upper_peak_wavelength + 10,
            'max'
        )

        self.assertEqual(
            upper_peak_wavelength, returned_upper_peak, 'Incorrect max peak')


class FindFeatureBounds(TestCase):
    """Tests for the ``find_feature_bounds`` function"""

    def test_bounds_for_simulated_feature(self):
        """Test correct boundaries are returned for a simulated feature"""

        # Note that we are simulating delta function absorption features
        wave = np.arange(7000, 8001)
        peak_wavelengths = (7100, 7500)
        flux, eflux = SimulatedSpectrum.delta_func(
            wave, peak_wave=peak_wavelengths)

        lower_peak_wavelength = min(peak_wavelengths)
        upper_peak_wavelength = max(peak_wavelengths)
        feature_dict = {
            'lower_blue': lower_peak_wavelength - 10,
            'upper_blue': lower_peak_wavelength + 10,
            'lower_red': upper_peak_wavelength - 10,
            'upper_red': upper_peak_wavelength + 10
        }

        feat_start, feat_end = spectra.find_feature_bounds(
            wave, flux, feature_dict)

        self.assertEqual(
            lower_peak_wavelength, feat_start, 'Incorrect min peak')

        self.assertEqual(
            upper_peak_wavelength, feat_end, 'Incorrect max peak')


class SampleFeatureProperties(TestCase):
    """Tests for the ``sample_feature_properties`` function"""

    @classmethod
    def setUpClass(cls):
        # Define test spectrum
        cls.wave = np.arange(1000, 2000)
        cls.observed_wave = np.mean(cls.wave)
        cls.rest_wave = cls.observed_wave - 100
        cls.flux, cls.error = SimulatedSpectrum.gaussian(
            cls.wave, mean=cls.observed_wave, stddev=100)

        # Select just the feature
        cls.feat_start = cls.wave[100]
        cls.feat_end = cls.wave[-100]
        feat_wave = cls.wave[100:-100]
        feat_flux = cls.flux[100:-100]

        # Calculate feature properties
        cls.area = spectra.feature_area(feat_wave, feat_flux)
        norm_flux, cls.pequiv_width = spectra.feature_pew(feat_wave, feat_flux)
        cls.velocity = spectra.feature_velocity(
            cls.rest_wave, feat_wave, feat_flux)

        cls.feat_name = 'test_CalcFeatureProperties'
        line_properties = {
            'restframe': cls.rest_wave,
        }

        spectra.line_locations[cls.feat_name] = line_properties

    @classmethod
    def tearDownClass(cls):
        del spectra.line_locations[cls.feat_name]

    def assertNumberSamples(self, nstep, nsamp=None):
        velocity, pequiv_width, area = spectra.sample_feature_properties(
            self.feat_name, self.feat_start, self.feat_end, self.wave,
            self.flux, nstep=nstep, debug=True)

        msg = 'Wrong number of samples for n={}'
        nsamp = ((2 * nstep) + 1) ** 2 if nsamp is None else nsamp
        self.assertEqual(nsamp, len(velocity), msg.format(nstep))
        self.assertEqual(nsamp, len(pequiv_width), msg.format(nstep))
        self.assertEqual(nsamp, len(area), msg.format(nstep))

    def test_number_of_samples(self):
        """Test the correct number of samples are performed"""

        for n in [0, 5]:
            self.assertNumberSamples(n)

        self.assertNumberSamples(0, 1)
        self.assertNumberSamples(-1, 0)

    def test_return_order(self):
        """Test values are returned in the correct order"""

        velocity, pequiv_width, area = spectra.sample_feature_properties(
            self.feat_name, self.feat_start, self.feat_end, self.wave,
            self.flux, nstep=0)

        msg = '{} values do not match.'
        self.assertAlmostEqual(
            self.velocity, velocity[0], msg=msg.format('Velocity'))

        self.assertAlmostEqual(
            self.pequiv_width, pequiv_width[0], msg=msg.format('pEW'))

        self.assertAlmostEqual(self.area, area[0], msg=msg.format('Area'))

    def test_uarray_support(self):
        """Test the function supports input arrays with ufloat objects"""

        uflux = uarray(self.flux, self.error)
        velocity, pequiv_width, area = spectra.sample_feature_properties(
            self.feat_name, self.feat_start, self.feat_end, self.wave, uflux,
            nstep=0)

        msg = '{} values do not match.'
        self.assertAlmostEqual(
            self.velocity, velocity[0], msg=msg.format('Velocity'))

        self.assertAlmostEqual(
            self.pequiv_width, pequiv_width[0], msg=msg.format('pEW'))

        self.assertAlmostEqual(self.area, area[0], msg=msg.format('Area'))


class SpectrumProperties(TestCase):
    """Tests for the ``_spectrum_properties`` function"""

    @classmethod
    def setUpClass(cls):
        # Define test spectrum
        cls.wave = np.arange(7000, 8000)
        cls.observed_wave = np.mean(cls.wave)
        cls.rest_wave = cls.observed_wave - 100
        cls.flux, error = SimulatedSpectrum.gaussian(
            cls.wave, mean=cls.observed_wave, stddev=100)

        line_properties = {
            'restframe': cls.rest_wave,
            'lower_blue': cls.wave[100],
            'upper_blue': cls.wave[100],
            'lower_red': cls.wave[-100],
            'upper_red': cls.wave[-100]
        }

        cls.old_lines = spectra._calc_properties.line_locations
        spectra._calc_properties.line_locations = \
            {'test_SpectrumProperties': line_properties}

    @classmethod
    def tearDownClass(cls):
        spectra._calc_properties.line_locations = cls.old_lines

    def test_spectrum_restframing(self):
        """Test properties of rest-framed spectra are independent of z"""

        z0_return = spectra._calc_properties._spectrum_properties(
            self.wave, self.flux, z=0, ra=1, dec=1)

        test_z = 1
        wave_at_z = self.wave * (1 + test_z)
        z1_return = spectra._calc_properties._spectrum_properties(
            wave_at_z, self.flux, z=test_z, ra=1, dec=1)

        self.assertListEqual(z0_return, z1_return)

    def test_extinction_dependence(self):
        """Test properties scale correctly with extinction"""

        self.fail()


class TabulateSpectrumProperties(TestCase):
    """Tests for the ``tabulate_spectral_properties`` function"""

    @classmethod
    def setUpClass(cls):
        # Define test spectrum
        cls.wave = np.arange(7000, 8000)
        cls.observed_wave = np.mean(cls.wave)
        cls.rest_wave = cls.observed_wave - 100
        cls.flux, cls.error = SimulatedSpectrum.gaussian(
            cls.wave, mean=cls.observed_wave, stddev=100)

        line_properties = {
            'restframe': cls.rest_wave,
            'lower_blue': cls.wave[100],
            'upper_blue': cls.wave[100],
            'lower_red': cls.wave[-100],
            'upper_red': cls.wave[-100]
        }

        cls.old_lines = spectra._calc_properties.line_locations
        spectra._calc_properties.line_locations = \
            {'test_TabulateSpectrumProperties': line_properties}

    @classmethod
    def tearDownClass(cls):
        spectra.line_locations = cls.old_lines

    def test_column_names(self):
        """Test the returned table has the correct column names"""

        expected_names = [
            'date',
            'test_TabulateSpectrumProperties_vel',
            'test_TabulateSpectrumProperties_vel_err',
            'test_TabulateSpectrumProperties_vel_samperr',
            'test_TabulateSpectrumProperties_pew',
            'test_TabulateSpectrumProperties_pew_err',
            'test_TabulateSpectrumProperties_pew_samperr',
            'test_TabulateSpectrumProperties_area',
            'test_TabulateSpectrumProperties_area_err',
            'test_TabulateSpectrumProperties_area_samperr'
        ]

        returned_table = spectra.tabulate_spectral_properties(
            [], [], [], [], [], [])

        self.assertListEqual(expected_names, returned_table.colnames)
