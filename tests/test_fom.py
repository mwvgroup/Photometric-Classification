#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``fom`` module."""

from unittest import TestCase

import numpy as np

from phot_class import fom


class FOM(TestCase):
    """Tests for the ``fom`` function"""

    def test_correct_classifications(self):
        """Test FOM is 1 for a set of completely correct classifications"""

        truth = ['type1', 'type1', 'type1', 'type2', 'type2']
        class_ = ['type1', 'type1', 'type1', 'type2', 'type2']

        type2_fom = fom.fom(truth, class_, check_type='type2')
        self.assertEqual(1, type2_fom)

    def test_incorrect_classifications(self):
        """Test FOM is 1 for a set of completely correct classifications"""

        truth = ['type1', 'type1', 'type1', 'type2', 'type2']
        class_ = ['type2', 'type2', 'type2', 'type1', 'type1']

        type2_fom = fom.fom(truth, class_, check_type='type2')
        self.assertEqual(0, type2_fom)

    def test_independent_of_extraneous_types(self):
        """Test the FOM is independent of types other than test_type"""

        truth = ['type1', 'type2', 'type3', 'test_type', 'test_type']
        class_ = ['type4', 'type5', 'type6', 'test_type', 'test_type']

        type2_fom = fom.fom(truth, class_, check_type='test_type')
        self.assertEqual(1, type2_fom)


def create_test_grid(xrange, yrange, fill_class):

    xx, yy = np.meshgrid(np.arange(*xrange), np.arange(*yrange))
    x = xx.flatten()
    y = yy.flatten()
    truth = np.full(len(x), fill_class)

    return x, y, truth


class RectangularFOM(TestCase):
    """Tests for the ``rectangular`` function"""

    def runTest(self):
        """Test the FOM is 1 when evaluated at a known, simulated boundary"""

        x, y, truth = create_test_grid((-10, 10), (-10, 10), 'type1')
        truth[(x > 0) & (y > 0)] = 'type2'

        type2_fom = fom.rectangular(truth, x, y, 0, 0, 'type2')
        type1_fom = fom.rectangular(truth, x, y, 0, 0, 'type1')

        self.assertEqual(type2_fom, 1)
        self.assertEqual(type1_fom, 0)


class VerticalFOM(TestCase):
    """Tests for the ``vertical`` function"""

    def runTest(self):
        """Test the FOM is 1 when evaluated at a known, simulated boundary"""
        
        x, y, truth = create_test_grid((-10, 10), (-10, 10), 'type1')
        truth[x > 0] = 'type2'

        type2_fom = fom.vertical(truth, x, 0, 'type2')
        type1_fom = fom.vertical(truth, x, 0, 'type1')

        self.assertEqual(type2_fom, 1)
        self.assertEqual(type1_fom, 0)


class HorizontalFOM(TestCase):
    """Tests for the ``horizontal`` function"""

    def runTest(self):
        """Test the FOM is 1 when evaluated at a known, simulated boundary"""
        
        x, y, truth = create_test_grid((-10, 10), (-10, 10), 'type1')
        truth[y > 0] = 'type2'

        type2_fom = fom.horizontal(truth, y, 0, 'type2')
        type1_fom = fom.horizontal(truth, y, 0, 'type1')

        self.assertEqual(type2_fom, 1)
        self.assertEqual(type1_fom, 0)


class LinearFOM(TestCase):
    """Tests for the ``linear`` function"""

    def runTest(self):
        """Test the FOM is 1 when evaluated at a known, simulated boundary"""
        
        m, b = .5, 2
        x, y, truth = create_test_grid((-10, 10), (-10, 10), 'type1')
        truth[y > m * x + b] = 'type2'

        type2_fom = fom.linear(truth, x, y, m, b, 'type2')
        type1_fom = fom.linear(truth, x, y, m, b, 'type1')

        self.assertEqual(type2_fom, 1)
        self.assertEqual(type1_fom, 0)


class DiagonalFOM(TestCase):
    """Tests for the ``diagonal`` function"""

    def runTest(self):
        """Test the FOM is 1 when evaluated at a known, simulated boundary"""
        
        b = 2
        x, y, truth = create_test_grid((-10, 10), (-10, 10), 'type1')
        truth[y > x + b] = 'type2'

        type2_fom = fom.diagonal(truth, x, y, b, 'type2')
        type1_fom = fom.diagonal(truth, x, y, b, 'type1')

        self.assertEqual(type2_fom, 1)
        self.assertEqual(type1_fom, 0)
