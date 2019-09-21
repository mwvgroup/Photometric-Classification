#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``fom`` module calculates the Figure of Merit (FOM) parameter which is
used to optimize the classification of targets as either 91bg or normal type
Ia.

Functions are provided for calculating the FOM according to multiple
classification schemes. This includes the following functions:

+----------------------+--------------------------------------------------------------------------------+
| Function             | Description                                                                    |
+======================+================================================================================+
| ``rectangular``      | Uses a pair of rectangular axes to classify targets in the                     |
|                      | upper right / lower left quadrant as 91bg / normal type Ia                     |
+----------------------+--------------------------------------------------------------------------------+
| ``vertical``         | Classifies targets using only the x coordinate (difference in red band chisq)  |
+----------------------+--------------------------------------------------------------------------------+
| ``horizontal``       | Classifies targets using only the y coordinate (difference in blue band chisq) |
+----------------------+--------------------------------------------------------------------------------+
| ``linear``           | Classifies targets using a diagonal boundary: y = mx + b                       |
+----------------------+--------------------------------------------------------------------------------+
| ``diagonal``         | Classifies targets using a diagonal boundary with a slope of 1: y = x + b      |
+----------------------+--------------------------------------------------------------------------------+
"""


def _fom(equals_truth, equals_class):
    """Calculate the figure of merit for a given set of classifications

    Args:
        equals_truth (ndarray): Whether the true classification of each
                                 object matches a given type
        equals_class (ndarray): Whether the derived classification of each
                                 object matches a given type

    Returns:
        (ntrue / ntotal) * (ntrue / (ntrue + nfalse))
    """

    ntotal = sum(equals_truth)
    ntrue = sum(equals_truth & equals_class)
    nfalse = sum(~equals_truth & equals_class)
    return (ntrue / ntotal) * (ntrue / (ntrue + nfalse))


def rectangular(truth, x, y, x_cutoff, y_cutoff, check_type='91bg'):
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

    type_matches_truth = (truth == check_type)
    type_matches_classification = ((x > x_cutoff) & (y > y_cutoff))
    return _fom(type_matches_truth, type_matches_classification)


def vertical(truth, x, x_cutoff, check_type='91bg'):
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

    type_matches_truth = (truth == check_type)
    type_matches_classification = (x > x_cutoff)
    return _fom(type_matches_truth, type_matches_classification)


def horizontal(truth, y, y_cutoff, check_type='91bg'):
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

    type_matches_truth = (truth == check_type)
    type_matches_classification = (y > y_cutoff)
    return _fom(type_matches_truth, type_matches_classification)


def linear(truth, x, y, m, b, check_type='91bg'):
    """Calculate the figure of merit using a diagonal lower bound: y = mx + b

    Args:
        truth  (ndarray): Array of classifications to take as truth
        x      (ndarray): The classification x coordinate
        y      (ndarray): The classification y coordinate
        m        (float): The slope of the boundary
        b        (float): The y-intercept of the boundary
        check_type (str): The classification to calculate the FOM for

    Returns:
    :math:`N_{total}` is the total
number of objects with that type    The figure of merit value
    """

    type_matches_truth = (truth == check_type)
    type_matches_classification = y > (m * x + b)
    return _fom(type_matches_truth, type_matches_classification)


def diagonal(truth, x, y, b, check_type='91bg'):
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
