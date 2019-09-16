#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``models`` module."""

from os import environ
from unittest import TestCase, skipIf

import numpy as np
import sncosmo

from phot_class import models

models.register_sources(force=True)
running_in_travis = 'TRAVIS' in environ


class TemplateLoading(TestCase):
    """Tests for models.load_template"""

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


class BaseSourceTestingClass(TestCase):
    """Tests for an arbitrary sncosmo Source object"""

    @classmethod
    def setUpClass(cls):
        cls.source = sncosmo.get_source(cls.source_name, cls.source_version)

    def _test_flux_matches_template(self):
        """Test return of model.flux agrees with the source template"""

        # Get template spanning the phase range of the model
        (stretch, color, phase, wave), template = models.load_template(
            self.source.minphase(), self.source.maxphase())

        for i, x1 in enumerate(stretch):
            for j, c in enumerate(color):
                self.source.set(x1=x1, c=c)
                model_flux = self.source.flux(phase, wave)
                template_flux = template[i, j]
                self.assertTrue(np.all(np.isclose(model_flux, template_flux)))

    def _test_phase_range(self, min_phase, max_phase):
        self.assertEqual(min_phase, self.source.minphase())
        self.assertEqual(max_phase, self.source.maxphase())

    def _test_correct_version(self):
        """Test the source was registered with the correct version name

        This ensures that ``sncosmo.get_source`` returns the correct version.
        """

        self.assertEqual(self.source.version, self.source_version)

    def _test_zero_flux_outside_phase_range(self):
        """Test the modeled flux outside the model's phase range is zero"""

        min_phase = self.source.minphase()
        late_phase_flux = self.source.bandflux('sdssg', min_phase - 1)
        late_phase_zero = np.isclose(0, late_phase_flux, atol=1e-6)
        self.assertTrue(late_phase_zero, f'Non-zero flux {late_phase_flux}')

        max_phase = self.source.minphase()
        early_phase_flux = self.source.bandflux('sdssg', max_phase + 1)
        early_phase_zero = np.isclose(0, early_phase_flux, atol=1e-6)
        self.assertTrue(early_phase_zero, f'Non-zero flux {early_phase_flux}')


class PhaseLimited(BaseSourceTestingClass):
    """Tests the 'phase_limited' version of the sn91bg model"""

    source_name = 'sn91bg'
    source_version = 'phase_limited'

    def test_flux_matches_template(self):
        """Test return of model.flux agrees with the source template"""

        self._test_flux_matches_template()

    def test_correct_version(self):
        """Test the source was registered with the correct version name"""

        self._test_correct_version()

    def test_phase_matches_salt2(self):
        """Test the model has the correct phase range"""

        self._test_phase_range(-18, 50)

    def test_zero_flux_outside_phase_range(self):
        """Test the modeled flux outside the modeled phase range is zero"""

        self._test_zero_flux_outside_phase_range()


class FullPhase(BaseSourceTestingClass):
    """Tests the 'phase_limited' version of the sn91bg model"""

    source_name = 'sn91bg'
    source_version = 'full_phase'

    def test_flux_matches_template(self):
        """Test return of model.flux agrees with the source template"""

        self._test_flux_matches_template()

    def test_correct_version(self):
        """Test the source was registered with the correct version name"""

        self._test_correct_version()

    def test_phase_matches_salt2(self):
        """Test the model has the correct phase range"""

        self._test_phase_range(-18, 100)

    def test_zero_flux_outside_phase_range(self):
        """Test the modeled flux outside the modeled phase range is zero"""

        self._test_zero_flux_outside_phase_range()
