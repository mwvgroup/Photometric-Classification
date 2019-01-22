#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-

"""This module combines a collection of individual SED template files into a
single .npy file.
"""

import os
from glob import glob

import numpy as np


def read_template_file(path):
    """Read a SED template text file and return data as a 2D array

    The first index of the returned array selects the day until maximum.
    The second index selects the wavelength.

    Args:
        path (str): Path of the file to read

    Returns:
        An 114 x 1101 array of flux values.
    """

    number_dates = 114
    number_wavelengths = 1101
    data = np.genfromtxt(path)
    return data[:, 2].reshape(number_dates, number_wavelengths)


def combine_sed_templates(out_path, in_dir):
    """Combine SED template files into a single array and save it
    to a .npy file.

    Template spans 7 stretch values, 5 color values, 114 dates,
    and 1101 wavelengths. Output shape is (7, 5, 114, 1101).

    Args:
        out_path (str): Path of the desired output file
        in_dir   (str): Directory of SED template files
    """

    path_list = glob(os.path.join(in_dir, '*.SED'))

    # Create array to store all SED templates
    number_of_stretch = 7
    number_of_color = 5
    number_of_dates = 114
    number_of_wavelengths = 1101
    out_data = np.empty(
        (number_of_stretch,
         number_of_color,
         number_of_dates,
         number_of_wavelengths)
    )

    # Add SED templates to out_data
    for path in sorted(path_list):
        model_params = os.path.basename(path).split('_')
        stretch = int(model_params[1].lstrip('ST'))
        color = int(model_params[2].lstrip('C').rstrip('.SED'))
        out_data[stretch, color] = read_template_file(path)

    np.save(out_path, out_data)
