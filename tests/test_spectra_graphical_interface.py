#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``simulation.spectra`` module."""

from unittest import TestCase

import numpy as np
from uncertainties.unumpy import uarray

from phot_class import spectra
from phot_class.spectra import SpectrumInspector
from .test_spectra_calc_properties import SimulatedSpectrum


# noinspection PyTypeChecker
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
        continuum, norm_flux, cls.pequiv_width = \
            spectra.feature_pew(feat_wave, feat_flux)

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

    def test_range_too_narrow(self):
        """Test an error is raised if too narrow a feature is specified"""

        args = dict(
            feat_name=self.feat_name,
            feat_start=self.wave[100],
            feat_end=self.wave[109],
            wave=self.wave,
            flux=self.flux
        )

        self.assertRaisesRegex(
            ValueError,
            'Range too small.*',
            SpectrumInspector._sample_feature_properties,
            **args
        )

    def assertNumberSamples(self, nstep, nsamp=None):
        """Assert the correct number of samples were performed for a given
         number of steps.

         If ``nsamp`` is not specified use:
             nsamp = ((2 * nstep) + 1) ** 2

        Args:
            nstep (int): Number of steps taken in each direction
            nsamp (int): Number of expected samples
        """

        velocity, pequiv_width, area = SpectrumInspector._sample_feature_properties(
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

        velocity, pequiv_width, area = SpectrumInspector._sample_feature_properties(
            self.feat_name, self.feat_start, self.feat_end, self.wave,
            self.flux, nstep=0)

        msg = '{} values do not match.'
        # Todo: Add back in the velocity calculation
        # self.assertAlmostEqual(
        #     self.velocity, velocity[0], msg=msg.format('Velocity'))

        self.assertAlmostEqual(
            self.pequiv_width, pequiv_width[0], msg=msg.format('pEW'))

        self.assertAlmostEqual(self.area, area[0], msg=msg.format('Area'))

    def test_uarray_support(self):
        """Test the function supports input arrays with ufloat objects"""

        uflux = uarray(self.flux, self.error)
        velocity, pequiv_width, area = SpectrumInspector._sample_feature_properties(
            self.feat_name, self.feat_start, self.feat_end, self.wave, uflux,
            nstep=0)

        msg = '{} values do not match.'
        # Todo: Add back in the velocity calculation
        # self.assertAlmostEqual(
        #     self.velocity, velocity[0], msg=msg.format('Velocity'))

        self.assertAlmostEqual(
            self.pequiv_width, pequiv_width[0], msg=msg.format('pEW'))

        self.assertAlmostEqual(self.area, area[0], msg=msg.format('Area'))


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
