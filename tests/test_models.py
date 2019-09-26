#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``models`` module."""

from os import environ
from unittest import TestCase

import numpy as np
import sncosmo

from phot_class import models

models.register_sources(force=True)
running_in_travis = 'TRAVIS' in environ


class TemplateLoading(TestCase):
    """Tests for the ``load_template`` function"""

    coord_names = ('stretch', 'color', 'phase', 'wave')

    # noinspection PyTypeChecker
    def test_correct_template_coordinates(self):
        """Test the returned template coordinates have the correct values"""

        coords, template = models.load_template()
        self.assertEqual(len(coords), len(self.coord_names),
                         'Wrong number of coordinate arrays')

        # Define expected coordinates
        stretch = np.array([0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25])
        color = np.array([0.0, 0.25, 0.5, 0.75, 1])
        phase = np.arange(-18, 101, 1)
        wave = np.arange(1000, 12001, 10)
        expected_coords = [stretch, color, phase, wave]

        test_data = zip(self.coord_names, expected_coords, coords)
        for coord_name, expected, returned in test_data:
            self.assertListEqual(expected.tolist(), returned.tolist(),
                                 f'Incorrect {coord_name} coordinates')

    def test_template_dimensions_match_coords(self):
        """Test the flux template dimensions match the number of coordinates"""

        coords, template = models.load_template()
        test_data = zip(self.coord_names, coords, template.shape)
        for coord_name, coord_array, template_len in test_data:
            self.assertEqual(
                len(coord_array), template_len,
                f'Shape mismatch for coordinate {coord_name}')

    def test_phase_limiting(self):
        """Test the returned phase matches the function arguments"""

        min_phase = -10
        max_phase = 50

        (_, _, phase, _), _ = models.load_template(min_phase, max_phase)
        self.assertEqual(min_phase, min(phase))
        self.assertEqual(max_phase, max(phase))

    def test_default_phase(self):
        """Test the default phase range extends over the full model"""

        min_model_phase = -18
        max_model_phase = 100

        (_, _, phase, _), _ = models.load_template()
        self.assertEqual(min_model_phase, min(phase))
        self.assertEqual(max_model_phase, max(phase))


class BisectSearch(TestCase):
    """Tests for the ``_bi_search`` function"""

    def test_element_in_array(self):
        """Test the correct index is returned for an element in an array"""

        test_array = np.array([1, 2, 3])
        test_elt_index = 1
        test_elt = test_array[test_elt_index]

        self.assertEqual(1, models._sources._bi_search(test_array, test_elt))

    def test_element_in_range(self):
        """Test the correct indices are returned for an element in the
        range of an array
        """

        test_array = np.array([1, 2, 3])
        test_elt = 1.5
        expected_indices = [0, 1]
        returned_indices = models._sources._bi_search(test_array, test_elt)
        self.assertSequenceEqual(expected_indices, returned_indices)

    def test_element_outside_range(self):
        """Test an error is raised for a parameter outside the array's range"""

        args = (np.array([1, 2, 3]), 5)
        self.assertRaises(RuntimeError, models._sources._bi_search, *args)

    def test_unsorted_array(self):
        """Test an error is raised when an array is not sorted"""

        args = (np.array([1, 3, 2]), 5)
        self.assertRaises(RuntimeError, models._sources._bi_search, *args)


class Registration(TestCase):
    """Models are correctly registered with sncosmo"""

    def test_source_loading(self):
        """Test the correct source object is returned for different args"""

        bg_full_phase = models._load_source('sn91bg', version='full_phase')
        self.assertEqual('sn91bg', bg_full_phase.name)
        self.assertEqual('full_phase', bg_full_phase.version)

        bg_salt2_phase = models._load_source('sn91bg', version='salt2_phase')
        self.assertEqual('sn91bg', bg_salt2_phase.name)
        self.assertEqual('salt2_phase', bg_salt2_phase.version)

        bg_hsiao_phase = models._load_source('sn91bg', version='hsiao_phase')
        self.assertEqual('sn91bg', bg_hsiao_phase.name)
        self.assertEqual('hsiao_phase', bg_hsiao_phase.version)

        hsiao_x1 = models._load_source('hsiao_x1', version='shouldn\'t matter')
        self.assertEqual('hsiao_x1', hsiao_x1.name)
        self.assertEqual('3.0.x1', hsiao_x1.version)

        # Test an error is raised for an unknown source
        args = 'dummy_name', 'dummy_version'
        self.assertRaises(ValueError, models._load_source, *args)

    def test_registration(self):
        """Test the correct default version is returned."""

        bg_default = sncosmo.get_source('sn91bg')
        self.assertEqual('salt2_phase', bg_default.version)

        bg_full_phase = sncosmo.get_source('sn91bg', version='full_phase')
        self.assertEqual('full_phase', bg_full_phase.version)

        bg_salt2_phase = sncosmo.get_source('sn91bg', version='salt2_phase')
        self.assertEqual('salt2_phase', bg_salt2_phase.version)

        bg_hsiao_phase = sncosmo.get_source('sn91bg', version='hsiao_phase')
        self.assertEqual('hsiao_phase', bg_hsiao_phase.version)

        hsiao_default = sncosmo.get_source('hsiao')
        self.assertEqual('3.0', hsiao_default.version)

        hsiao_x1 = sncosmo.get_source('hsiao_x1')
        self.assertEqual('3.0.x1', hsiao_x1.version)


class BaseSourceTestingClass(TestCase):
    """Tests for an arbitrary sncosmo Source object"""

    @classmethod
    def setUpClass(cls):
        cls.source = sncosmo.get_source(cls.source_name, cls.source_version)

    def _test_phase_range(self, min_phase, max_phase):
        """Test the model has the expected phase range"""

        self.assertEqual(min_phase, self.source.minphase())
        self.assertEqual(max_phase, self.source.maxphase())

    def _test_zero_flux_outside_phase_range(self):
        """Test the modeled flux outside the model's phase range is zero"""

        min_phase = self.source.minphase()
        early_phase_flux = self.source.bandflux('sdssg', min_phase - 1)
        late_phase_zero = np.isclose(0, early_phase_flux, atol=1e-6)
        self.assertTrue(late_phase_zero, f'Non-zero flux {early_phase_flux}')

        max_phase = self.source.minphase()
        late_phase_flux = self.source.bandflux('sdssg', max_phase + 1)
        early_phase_zero = np.isclose(0, late_phase_flux, atol=1e-6)
        self.assertTrue(early_phase_zero, f'Non-zero flux {late_phase_flux}')

    def _test_flux_at_coords_matches_template(self):
        """Test return of model.flux agrees with the source template"""

        # Get template spanning the phase range of the model
        (stretch, color, phase, wave), template = models.load_template(
            self.source.minphase(), self.source.maxphase())

        for i, x1 in enumerate(stretch):
            for j, c in enumerate(color):
                self.source.set(x1=x1, c=c)
                model_flux = self.source.flux(phase, wave)
                template_flux = template[i, j]
                is_close = np.isclose(model_flux, template_flux)
                self.assertTrue(is_close.all())


class SN91bgSalt2Phase(BaseSourceTestingClass):
    """Tests the 'phase_limited' version of the sn91bg model"""

    source_name = 'sn91bg'
    source_version = 'salt2_phase'

    def test_flux_matches_template(self):
        """Test return of model.flux agrees with the source template at
        the template coordinates.
        """

        self._test_flux_at_coords_matches_template()

    def test_correct_phase_range(self):
        """Test the model has the correct phase range"""

        self._test_phase_range(-18, 50)

    def test_zero_flux_outside_phase_range(self):
        """Test the modeled flux outside the modeled phase range is zero"""

        self._test_zero_flux_outside_phase_range()


class SN91bgHsiaoPhase(BaseSourceTestingClass):
    """Tests the 'phase_limited' version of the sn91bg model"""

    source_name = 'sn91bg'
    source_version = 'hsiao_phase'

    def test_flux_matches_template(self):
        """Test return of model.flux agrees with the source template at
        the template coordinates.
        """

        self._test_flux_at_coords_matches_template()

    def test_correct_phase_range(self):
        """Test the model has the correct phase range"""

        self._test_phase_range(-18, 85)

    def test_zero_flux_outside_phase_range(self):
        """Test the modeled flux outside the modeled phase range is zero"""

        self._test_zero_flux_outside_phase_range()


class SN91bgFullPhase(BaseSourceTestingClass):
    """Tests the 'full_phase' version of the sn91bg model"""

    source_name = 'sn91bg'
    source_version = 'full_phase'

    def test_flux_matches_template(self):
        """Test return of model.flux agrees with the source template"""

        self._test_flux_at_coords_matches_template()

    def test_correct_phase_range(self):
        """Test the model has the correct phase range"""

        self._test_phase_range(-18, 100)

    def test_zero_flux_outside_phase_range(self):
        """Test the modeled flux outside the modeled phase range is zero"""

        self._test_zero_flux_outside_phase_range()


class HsiaoX1(BaseSourceTestingClass):
    """Tests the 'phase_limited' version of the sn91bg model"""

    source_name = 'hsiao_x1'
    source_version = '3.0.x1'

    def test_correct_phase_range(self):
        """Test the model has a phase range matching the original model"""

        hsiao = sncosmo.get_source('hsiao')
        self._test_phase_range(hsiao.minphase(), hsiao.maxphase())

    def test_matches_original_model(self):
        """Test the flux of our custom model at x1=default matches hsiao 3.0"""

        custom = sncosmo.get_source(self.source_name, self.source_version)
        original = sncosmo.get_source('hsiao', '3.0')

        phase_arr = np.arange(original.minphase(), original.maxphase())
        wave_arr = np.arange(original.minwave(), original.maxwave())

        custom_flux = custom.flux(phase_arr, wave_arr)
        original_flux = original.flux(phase_arr, wave_arr)
        self.assertListEqual(original_flux.tolist(), custom_flux.tolist())
