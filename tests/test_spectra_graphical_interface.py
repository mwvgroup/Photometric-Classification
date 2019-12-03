#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``spectra.graphical_interface`` module."""

from unittest import TestCase

import matplotlib as mpl
import numpy as np
from astropy.table import Table

from phot_class import spectra
from phot_class.spectra import SpectrumInspector
from .test_spectra_calc_properties import SimulatedSpectrum


class TestPrepareSpectrum(TestCase):

    @classmethod
    def setUp(cls):
        """Simulate a spectrum and create a ``SpectrumInspector`` object"""

        # Define test spectrum
        wave = np.arange(7000, 8000)
        flux, error = SimulatedSpectrum.gaussian(wave, stddev=100)
        spectrum = Table([wave, flux], names=['wavelength', 'flux'])
        spectrum['date'] = ['some date']
        spectrum.meta['z'] = 0
        spectrum.meta['ra'] = 30
        spectrum.meta['dec'] = -29
        spectrum.meta['obj_id'] = 'dummy_id'

        cls.inspector = SpectrumInspector(spectrum)

    def test_prepare_spectrum_sets_attributes(self):
        """Test ``prepare_spectrum`` save results to instance attributes

        Checked attributes include ``bin_wave``, ``bin_flux``, ``rest_wave``,
        and ``rest_flux``.
        """

        self.inspector.prepare_spectrum(5, 'avg')
        self.assertIsNotNone(self.inspector.bin_wave)
        self.assertIsNotNone(self.inspector.bin_flux)
        self.assertIsNotNone(self.inspector.rest_wave)
        self.assertIsNotNone(self.inspector.rest_flux)

    def test_successive_preparation_error(self):
        """Test successive ``prepare_spectrum`` calls raise an error"""

        args = (5, 'avg')
        self.inspector.prepare_spectrum(*args)
        self.assertRaises(ValueError, self.inspector.prepare_spectrum, *args)

    def test_interactive_plotting_on(self):
        """Check interactive plotting is turned on"""

        self.assertTrue(mpl.is_interactive())


# noinspection PyTypeChecker
class SampleFeatureProperties(TestCase):
    """Tests for the ``sample_feature_properties`` function"""

    @classmethod
    def setUpClass(cls):
        # Define test spectrum
        wave = np.arange(7000, 9000)
        flux, error = SimulatedSpectrum.gaussian(wave, stddev=100)
        spectrum = Table([wave, flux], names=['wavelength', 'flux'])
        spectrum['date'] = ['some date']
        spectrum.meta['z'] = 0
        spectrum.meta['ra'] = 30
        spectrum.meta['dec'] = -29
        spectrum.meta['obj_id'] = 'dummy_id'

        cls.inspector = SpectrumInspector(spectrum)
        cls.bin_size = 3
        cls.inspector.prepare_spectrum(cls.bin_size, 'avg')

        # Define test feature
        cls.feat_name = 'test_CalcFeatureProperties'
        line_properties = dict(
            restframe=np.average(wave),
            feature_id='dummy_feature_id'
        )

        spectra.line_locations[cls.feat_name] = line_properties

    @classmethod
    def tearDownClass(cls):
        del spectra.line_locations[cls.feat_name]

    # Todo: test range changes with binning
    def test_range_too_narrow(self):
        """Test an error is raised if too narrow a feature is specified"""

        args = dict(
            feat_name=self.feat_name,
            feat_start=self.inspector.rest_wave[100],
            feat_end=self.inspector.rest_wave[109],
            debug=True
        )

        self.assertRaisesRegex(
            ValueError,
            'Range too small.*',
            self.inspector._sample_feature_properties,
            **args
        )

    def assertNumberSamples(self, nstep, nsamp):
        """Assert the correct number of samples were performed for a given
         number of steps.

         If ``nsamp`` is not specified use:
             nsamp = ((2 * nstep) + 1) ** 2

        Args:
            nstep (int): Number of steps taken in each direction
            nsamp (int): Number of expected samples
        """

        velocity, pew, area = self.inspector._sample_feature_properties(
            feat_name=self.feat_name,
            feat_start=self.inspector.rest_wave[100],
            feat_end=self.inspector.rest_wave[-100],
            nstep=nstep, return_samples=True, debug=True)

        msg = 'Wrong number of samples for n={}'
        self.assertEqual(nsamp, len(velocity), msg.format(nstep))
        self.assertEqual(nsamp, len(pew), msg.format(nstep))
        self.assertEqual(nsamp, len(area), msg.format(nstep))

    def test_number_of_samples(self):
        """Test the correct number of samples are performed"""

        for n in [0, 5]:
            nsamp = ((2 * n) + 1) ** 2
            self.assertNumberSamples(n, nsamp)

        self.assertNumberSamples(-1, 0)


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

        cls.old_lines = spectra.calc_properties.line_locations
        spectra.calc_properties.line_locations = \
            {'test_TabulateSpectrumProperties': line_properties}

    @classmethod
    def tearDownClass(cls):
        spectra.line_locations = cls.old_lines

    def test_column_names(self):
        """Test the returned table has the correct column names"""

        expected_names = [
            'obj_id',
            'sid',
            'date',
            'type',
            'feat_name',
            'feat_start',
            'feat_end',
            'vel',
            'vel_err',
            'vel_samperr',
            'pew',
            'pew_err',
            'pew_samperr',
            'area',
            'area_err',
            'area_samperr',
            'msg'
        ]

        returned_table = spectra.tabulate_spectral_properties(iter([]))
        self.assertListEqual(expected_names, returned_table.colnames)
