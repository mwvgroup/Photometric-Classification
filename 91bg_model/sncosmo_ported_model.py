#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

"""This module defines a custom 91bg model for use with SNCosmo"""

import os
from glob import glob

import numpy as np
import sncosmo
from scipy.interpolate import interpn

FILE_DIR = os.path.abspath(os.path.dirname(__file__))
MODEL_SOURCE_DIR = os.path.join(FILE_DIR, 'sed_templates')
COMPILED_MODEL_PATH = os.path.join(FILE_DIR, 'complete_template.npy')


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


def combine_sed_templates(
        out_path=COMPILED_MODEL_PATH, in_dir=MODEL_SOURCE_DIR):
    """Combine SED template files into a single array and save it
    to a .npy file.

    Template spans 7 stretch values, 5 color values, 114 dates,
    and 1101 wavelengths. Output shape is (7, 5, 114, 1101).

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

    def __init__(self):
        super(SN91bgSource, self).__init__()

        self.name = '91bg model'
        self.version = 'None'

        # Model spectra and initial parameters
        self._flux_values = self._get_91bg_model()
        self._parameters = np.array([0, 0, 1])

        # Create an array of xi points for self._flux_values
        self.stretch_vals = np.arange(0.65, 1.26, .1)
        self.color_vals = np.array([0.0, 0.25, 0.5, 0.75, 1])
        self._phase = np.arange(-13, 101, 1)
        self._wave = np.arange(1000, 12001, 10)
        self._flux_points = [self.stretch_vals,
                             self.color_vals,
                             self._phase,
                             self._wave]

    @staticmethod
    def _get_91bg_model():
        """Load template spectra for SN 1991bg

        Template spans 7 stretch values, 5 color values, 114 dates,
        and 1101 wavelengths.

        Returns:
            An array of shape (7, 5, 114, 1101)
        """

        if not os.path.exists(COMPILED_MODEL_PATH):
            print('Compiled 91bg template not found. Creating a new template.')
            combine_sed_templates()
            print('Done')

        print(np.load(COMPILED_MODEL_PATH).shape)
        return np.load(COMPILED_MODEL_PATH)

    def _flux(self, phase, wave):
        stretch, color, amplitude = self._parameters
        requested_xi = [stretch, color, phase, wave]
        flux = interpn(points=self._flux_points,
                       values=self._flux_values,
                       xi=requested_xi)

        return amplitude * flux


if __name__ == '__main__':
    # Test our custom source on the example data
    # The example data aren't 91bg light curves - we are only looking for
    # the code to execute

    model = sncosmo.Model(source=SN91bgSource())
    data = sncosmo.load_example_data()
    result, fitted_model = sncosmo.fit_lc(data, model,
                                          ['stretch', 'color', 'amplitude'])

    sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)
