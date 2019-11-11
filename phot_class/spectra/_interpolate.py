# !/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module handles the interpolation of equivalent width values to the
time of B-band maximum using the prescription  of Folatelli et al. 2013.
"""

import numpy as np


# Todo: Write the logic for interpolating pew
# noinspection PyUnusedLocal
def interpolate_bmax_pew(time, ew_vals, tmax):
    """Interpolate time dependent equivalent widths

    Interpolate equivalent width values to determine the value at the time
    of B band maximum.

    From Folatelli et al. 2013:
        When several measurements are available within one week before and
        after the time of maximum, a smooth polynomial fit is used. If only
        two data points were obtained encompassing maximum light within −4
        and +4 days, then an interpolation was performed. If only data before
        or after maximum is available within [−7, +7] days, an extrapolation
        is performed if the closest point to maximum light is not farther than
        one day. For observations which only have one measurement in the
        range [−4, +4] days or when measurements in that range do not
        encompass maximum light, average slopes determined from the best
        observed SNe are used for correcting the pW value to the time of
        maximum light.

    Args:
        time (ndarray): A list of observation times for each equivalent width
        ew_vals (list): A list of equivalent width values
        tmax   (float): Time of maximum to interpolate for

    Returns:
        The equivalent width interpolated at ``tmax``
    """

    delta_t = np.abs(time - tmax)
    delta_t_less_than_7 = delta_t[time - tmax < 2]

    # When several measurements are available within one week before and after
    # the time of maximum, a smooth polynomial fit is used.
    if sum(delta_t <= 7) > 2:
        pass

    # If only two data points were obtained encompassing maximum light within
    # −4 and +4 days, then an interpolation was performed.
    elif sum(delta_t <= 4) == 2:
        pass

    # If only data before or after maximum is available within [−7, +7] days,
    # an extrapolation is performed if the closest point to maximum light is
    # not farther than one day.
    elif (all(delta_t_less_than_7 < 0) or all(delta_t_less_than_7 > 0)) and (
            min(delta_t) <= 1):
        pass

    # For observations which only have one measurement in the range [−4, +4]
    # days or when measurements in that range do not encompass maximum light,
    # average slopes determined from the best observed SNe are used for
    # correcting the pW value to the time of maximum light.
    elif sum(delta_t <= 4) == 1:
        pass
