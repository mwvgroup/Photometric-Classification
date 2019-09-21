#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``fom`` module calculates the Figure of Merit (FOM) parameter
defined as:

.. math::

    FOM = \frac{N_{true}}{N_{total}} \times \frac{N_{true}}{N_{true} + N_{false}}

where :math:`N_{true}` is the number of objects correctly classified as a given
type, :math:`N_{false}` is the number of objects falsely classified as the
given type, and :math:`N_{total}` is the total number of objects with that
type.

This module specifically focuses on classification techniques that classify
targets based on their (x, y) position in some two dimensional phase space.
Functions are provided for calculating the FOM according to multiple
classification schemes:

+----------------------+------------------------------------------------------+
| Function             | Classification Scheme                                |
+======================+======================================================+
| ``rectangular``      | Uses a pair of rectangular axes to classify targets  |
|                      | in the upper right quadrant as a given type.         |
+----------------------+------------------------------------------------------+
| ``vertical``         | Classifies targets with an x coordinate greater than |
|                      | some value as a given type.                          |
+----------------------+------------------------------------------------------+
| ``horizontal``       | Classifies targets with a y coordinate greater than  |
|                      | some value as a given type.                          |
+----------------------+------------------------------------------------------+
| ``linear``           | Classifies targets above a linear boundary           |
|                      | y = mx + b as a given type.                          |
+----------------------+------------------------------------------------------+
| ``diagonal``         | Classifies targets above a diagonal boundary         |
|                      | y = x + b as a given type. Note the difference in    |
|                      | slope between this option and the ``linear`` option. |
+----------------------+------------------------------------------------------+

Alternatively, the ``fom`` method can be used to calculate the FOM for a set
of predetermined classifications, independent of the chosen classification
technique (i.e. using the results of an already run classification).

Usage Example
-------------

>>> from phot_class.fom import fom
>>>
>>> # The true "spectroscopic" classification of an object
>>> truth = ['normal', 'normal', 'normal', '91bg', '91bg']
>>>
>>> # The result of some classification piepline
>>> classifications = ['normal', 'normal', '91bg', 'normal', '91bg']
>>>
>>> # Classify the FOM for the "91bg" type
>>> bg_fom = fom(truth, classifications, '91bg')
>>>
>>> # Classify the FOM for the "normal" type
>>> normal_fom = fom(truth, classifications, 'normal')

Function Documentation
----------------------
"""

import numpy as np


def fom(truth, classification, check_type):
    """Calculate the figure of merit for a given set of classifications

    Args:
        truth          (ndarray): The true classifications
        classification (ndarray): The assigned classifications
        check_type         (str): The type to calculate the FOM for

    Returns:
        (n_true / n_total_type) * (n_true / (n_true + n_false))
    """

    truth = np.asarray(truth)
    classification = np.asarray(classification)

    # The total number of objects of the given type
    ntotal = sum(truth == check_type)

    # The number of objects correctly classified as the given type
    ntrue = sum((truth == classification) & (classification == check_type))

    # The number of objects incorrectly classified as the given type
    nfalse = sum((truth != classification) & (classification == check_type))

    return (ntrue / ntotal) * (ntrue / (ntrue + nfalse))


def rectangular(truth, x, y, x_cutoff, y_cutoff, check_type):
    """Calculate the figure of merit using a pair of rectangular lower bounds

    Expected classifications are "normal" or "91bg" (case insensitive).
    All other classifications are ignored.
    
    Args:
        truth   (ndarray): Array of classifications (str) to take as truth
        x       (ndarray): The classification x coordinate
        y       (ndarray): The classification y coordinate
        x_cutoff  (float): The x boundary to use for separating 91bg / normal
        y_cutoff  (float): The y boundary to use for separating 91bg / normal
        check_type  (str): The classification to calculate the FOM for

    Returns:
        The figure of merit value
    """

    classification = np.full_like(truth, f'not_{check_type}')
    classification[(x > x_cutoff) & (y > y_cutoff)] = check_type
    return fom(truth, classification, check_type)


def horizontal(truth, x, x_cutoff, check_type):
    """Calculate the figure of merit using a vertical lower boundary
    (i.e. using only the x coordinate)

    Expected classifications are "normal" or "91bg" (case insensitive).
    All other classifications are ignored.

    Args:
        truth  (ndarray): Array of classifications to take as truth
        x       (ndarray): The classification x coordinate
        x_cutoff  (float): The x boundary to use for separating 91bg / normal
        check_type  (str): The classification to calculate the FOM for

    Returns:
        The figure of merit value
    """

    classification = np.full_like(truth, f'not_{check_type}')
    classification[x > x_cutoff] = check_type
    return fom(truth, classification, check_type)


def vertical(truth, y, y_cutoff, check_type):
    """Calculate the figure of merit using a horizontal lower boundary
    (i.e. using only the y coordinate)

    Expected classifications are "normal" or "91bg" (case insensitive).
    All other classifications are ignored.

    Args:
        truth  (ndarray): Array of classifications to take as truth
        y      (ndarray): The classification y coordinate
        y_cutoff (float): The y boundary to use for separating 91bg / normal
        check_type (str): The classification to calculate the FOM for

    Returns:
        The figure of merit value
    """

    classification = np.full_like(truth, f'not_{check_type}')
    classification[y > y_cutoff] = check_type
    return fom(truth, classification, check_type)


def linear(truth, x, y, m, b, check_type):
    """Calculate the figure of merit using a diagonal lower bound: y = mx + b

    Args:
        truth  (ndarray): Array of classifications to take as truth
        x      (ndarray): The classification x coordinate
        y      (ndarray): The classification y coordinate
        m        (float): The slope of the boundary
        b        (float): The y-intercept of the boundary
        check_type (str): The classification to calculate the FOM for

    Returns:
        The figure of merit value
    """

    classification = np.full_like(truth, f'not_{check_type}')
    classification[y > (m * x + b)] = check_type
    return fom(truth, classification, check_type)


def diagonal(truth, x, y, b, check_type):
    """Calculate the figure of merit using a diagonal lower bound: y = x + b

    Expected classifications are "normal" or "91bg" (case insensitive).
    All other classifications are ignored.

    Args:
        truth  (ndarray): Array of classifications to take as truth
        x      (ndarray): The classification x coordinate
        y      (ndarray): The classification y coordinate
        b        (float): The y-intercept of the boundary
        check_type (str): The classification to calculate the FOM for

    Returns:
        The figure of merit value
    """

    return linear(truth, x, y, m=1, b=b, check_type=check_type)
