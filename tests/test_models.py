#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``models`` module."""

from unittest import TestCase

import numpy as np
import sncosmo

from analysis_pipeline import models

models.register_sources(force=True)


class TestSourceFluxTemplate(TestCase):
    """Test sources return flux values in agreement with their templates"""

    def _test_source(self, source):
        """Test return of model.flux agrees with the source template

        Args:
            source (Source): The sncosmo source object to test
        """

        stretch = source._stretch
        color = source._color
        phase = source._phase
        wave = source._wave
        template = source._template

        for i, x1 in enumerate(stretch):
            for j, c in enumerate(color):
                source.set(x1=x1, c=c)
                model_flux = source.flux(phase, wave)
                template_flux = template[i, j]
                self.assertTrue(np.all(np.isclose(model_flux, template_flux)))

    def test_salt2_phase(self):
        """Test return of model.flux agrees with the source template"""

        source = sncosmo.get_source('sn91bg', 'salt2_phase')
        self._test_source(source)

    def test_color_interpolation(self):
        """Test return of model.flux agrees with the source template"""
        source = sncosmo.get_source('sn91bg', 'color_interpolation')
        self._test_source(source)


class TestVersionAgreement(TestCase):
    """Test the various versions of the ported models agree with each other"""

    def runTest(self):
        """Compare modeled fluxes for each stretch, color, phase, and wavelength
        in the flux templates.
        """

        # Load models for different source versions
        salt2_phase = sncosmo.Model(sncosmo.get_source('sn91bg', 'salt2_phase'))
        color_interp = sncosmo.Model(sncosmo.get_source('sn91bg', 'color_interpolation'))

        # Read in the coordinates of the flux template from file
        (stretch, color, phase, wave), _ = color_interp.source.get_template()

        for x1 in stretch:
            for c in color:
                salt2_phase.set(x1=x1, c=c)
                color_interp.set(x1=x1, c=c)

                salt2_phase_flux = salt2_phase.flux(phase, wave)
                color_interp_phase_flux = color_interp.flux(phase, wave)
                models_agree = np.all(salt2_phase_flux == color_interp_phase_flux)
                self.assertTrue(models_agree, f'Models disagree for x1={x1}, c={c}')


# Todo: This might not be a necessary test.
class TestInterpolation(TestCase):
    """Test models are correctly interpolating from their flux templates"""

    def test_salt2_phase(self):
        self.fail()

    def test_color_interpolation(self):
        self.fail()
