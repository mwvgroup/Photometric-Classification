#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Utilities used by different versions of SN 1991bg source classes"""

from bisect import bisect


def bi_search(a, x):
    """Binary search for value ``x`` in array ``a``

    Args:
        a (ndarray): The sorted list in which the number x will be searched
        x     (num): The number to be searched

    Returns:
        The position of nearest left neighbor of ``x``
        The position of nearest right neighbor of ``x``
    """

    index = bisect(a, x)
    return index - 1, index


def linear_interp(x0, x1, f, x):
    """Linear interpolation

    Args:
        x0, x1  (num): Grid points
        x       (num): Coordinates of interpolation point
        f      (list): The list contains values at x0 and x1

    Returns:
        Interpolated value
    """

    return f[0] + ((x - x0) * f[1] - (x - x0) * f[0]) / (x1 - x0)
