#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module defines a custom 91bg model for use with SNCosmo"""

import os
from glob import glob

import numpy as np
import sncosmo
from scipy.interpolate import interpn

__file__ = './'
FILE_DIR = os.path.abspath(os.path.dirname(__file__))
MODEL_DIR = os.path.join(FILE_DIR, 'sed_templates')


def read_template_file(path):
    """Read a SED template text file and return data as a 2D array

    The first index of the returned array selects the day untill maximum.
    The second index selects the wavlength.

    Args:
        path (str): Path of the file to read

    Returns:
        An 114 x 1101 array of flux values.
    """

    number_dates = 114
    number_wavelengths = 1101
    data = np.genfromtxt(path)
    return data[:, 2].reshape(number_dates, number_wavelengths)


def combine_sed_templates(out_path, in_dir=MODEL_DIR):
    """Combine SED template files into a single array and save it
    to a .npy file.

    Args:
        out_path (str): Path of the desired output file
        in_dir   (str): Directory of SED template files
    """

    path_list = glob(os.path.join(in_dir, '*.SED'))
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

    for path in sorted(path_list):
        model_params = os.path.basename(path).split('_')
        stretch = int(model_params[1].lstrip('ST'))
        color = int(model_params[2].lstrip('C').rstrip('.SED'))
        out_data[stretch, color] = read_template_file(path)

    np.save(out_path, out_data)


class SN91bgSource(sncosmo.Source):
    _param_names = ['stretch', 'color', 'amplitude']
    param_names_latex = ['x1', 'c', 'a']

    def __init__(self, model_path, name=None, version=None):
        self.name = '91bg model'
        self.version = 'None'

        self.stretch_vals = np.arange(0.65, 1.26, .1)
        self.color_vals = np.array([0.0, 0.25, 0.5, 0.75, 1])
        self.day_values = np.arange(-13, 101, 1)
        self.wavelengths = np.arange(1000, 12001, 10)

        self._flux_points = \
            [self.stretch_vals, self.color_vals, self.day_values,
             self.wavelengths]
        self._flux_values = np.load(model_path)

    def _flux(self, stretch, color, day, wavelength, amplitude):
        requested_xi = [stretch, color, day, wavelength]
        flux = interpn(points=self._flux_points,
                       values=self._flux_values,
                       xi=requested_xi)

        return amplitude * flux


model_path = 'complete_template.npy'
if not os.path.exists(model_path):
    print('91bg template not found. Compiling a new template.')
    combine_sed_templates(model_path)
