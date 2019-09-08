#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``utils`` module."""

import time
from unittest import TestCase

import numpy as np
import sncosmo
from astropy.table import Column, Table
from sndata.csp import dr3

from analysis_pipeline import utils


class TestTimeout(TestCase):
    """Tests for utils.timeout"""

    @staticmethod
    def sleep(seconds):
        """Run time.sleep wrapped with utils.timeout for a given seconds

        Args:
            seconds (int): Number of seconds to sleep for
        """

        with utils.timeout(seconds):
            time.sleep(seconds + 1)

    def test_raises_timeout_error(self):
        """Test a TimeoutError is raised"""

        self.assertRaises(TimeoutError, self.sleep, 2)

    def test_context_duration(self):
        """Test an error is raised after the correct amount of time"""

        seconds = 1
        start = time.time()
        try:
            self.sleep(seconds)

        except TimeoutError:
            duration = np.round(time.time() - start, 1)
            self.assertEqual(seconds, duration)

    def test_error_for_float_arg(self):
        """Test an error is raised by the context manager when passed a float

        The signal manager expects integers, so if a TypeError is not raised
        for a float then make sure you know what is actually going on.
        """

        self.assertRaises(TypeError, self.sleep, 2.1)


# Todo: Write a test that checks we get the correct value for chisq
class TestCalcModelChisq(TestCase):
    """Tests for analysis_pipeline.lc_fitting.calc_chisq"""

    @property
    def model(self):
        """The model to use for testing"""

        return sncosmo.Model('salt2')

    @property
    def data(self):
        """The data to use for testing"""

        # Type cast the 'band' column to U15 so we can use longer band names
        data = sncosmo.load_example_data()
        data['band'] = Column(data['band'], dtype='U15')
        return data

    # noinspection PyTypeChecker
    def test_arg_mutation(self):
        """Test that the input model and data table are not mutated"""

        model = self.model
        data = self.data
        utils.calc_model_chisq(data, model)

        self.assertTrue(all(data == self.data), 'Data table was mutated')
        self.assertListEqual(
            list(self.model.parameters),
            list(model.parameters),
            'Model parameters were mutated')

    def test_out_of_range_data(self):
        """Test that out of range phase values and bands are dropped"""

        # The sncosmo docs assure us that the example data is within
        # range of the salt2 model
        data = self.data
        expected_chisq, expected_dof = utils.calc_model_chisq(data, self.model)

        # Add values that are out of the model's phase range
        # Columns: 'time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'
        data.add_row([100000000, 'sdssu', 0, 0, 25, 'ab'])

        # Add values that are out of the model's wavelength range
        # Here we use the H band from CSP
        dr3.register_filters(force=True)
        data.add_row([1, 'csp_dr3_H', 0, 0, 25, 'ab'])

        returned_chisq, returned_dof = utils.calc_model_chisq(data,
                                                              self.model)

        self.assertEqual(expected_chisq, returned_chisq,
                         'Incorrect chisq value')
        self.assertEqual(expected_dof, returned_dof,
                         'Incorrect number of data points')

    def test_unregistered_bands(self):
        """Test unregistered band names in the data table cause an error
        instead of being dropped.
        """

        # Add values with an unregistered filter
        data = self.data
        data['band'][0] = 'made up band'
        self.assertRaises(
            Exception, utils.calc_model_chisq, data, self.model)

    def test_unsorted_times(self):
        """Test the function return is independent of the time order"""

        data = self.data
        data.sort('time')
        initial_chisq = utils.calc_model_chisq(data, self.model)

        data[0], data[-1] = data[-1], data[0]
        new_chisq = utils.calc_model_chisq(data, self.model)
        self.assertEqual(initial_chisq, new_chisq)


class TestSplitBands(TestCase):
    """Tests for utils.split_bands"""

    def runTest(self):
        """Test dummy bands are correctly seperated into blue and red"""

        # Define dummy bands and their effective wavelengths
        band_names = ['blue1', 'red1', 'red2']
        lambda_eff = [4000, 5500, 6000]
        expected_blue_bands = band_names[:1]
        expected_red_bands = band_names[1:]

        blue_bands, red_bands = utils.split_bands(band_names, lambda_eff)
        self.assertCountEqual(expected_blue_bands, blue_bands)
        self.assertCountEqual(expected_red_bands, red_bands)


class TestSplitData(TestCase):
    """Tests for utils.split_data"""

    def test_redshift_dependency(self):
        """Assert whether correct bands were returned for a given redshift

        Tested for redshifts 0, .18, .55, and .78. Runs test with the cutoff
        wavelength set to infinity.
        """

        band_names = np.array(['u', 'g', 'r', 'i', 'z'])
        lambda_eff = np.array([3550, 4680, 6160, 7480, 8930])

        for z in (0, .18, .55, .78):
            rest_frame_cutoff = 5500 * (1 + redshift)
            expected_blue = band_names[lambda_eff < rest_frame_cutoff]
            expected_red = band_names[lambda_eff > rest_frame_cutoff]

            data = Table([band_names], names=['band'])
            blue_table, red_table = utils.split_data(
                data, band_names, lambda_eff, redshift, cutoff=float('inf'))

            err_msg = f'Wrongs bands for z={redshift}'
            print(list(blue_table['band']), list(red_table['band']))
            self.assertCountEqual(expected_blue, blue_table['band'], err_msg)
            self.assertCountEqual(expected_red, red_table['band'], err_msg)

    def test_maintains_metadata(self):
        """Test whether passed and returned tables have same metadata"""

        test_data = Table([['u', 'g']], names=['band'])
        test_data.meta['dummy_key'] = 12345
        blue_data, red_data = utils.split_data(
            test_data, ['u', 'g'], [3550, 4680], 0)

        self.assertIn('dummy_key', blue_data.meta,
                      'Blue table missing metadata')
        self.assertIn('dummy_key', red_data.meta, 'Red table missing metadata')

    def test_missing_effective_wavelength(self):
        """Test a value error is raised for missing effective wavelengths"""

        test_data = Table([['u', 'g']], names=['band'])
        args = test_data, ['u'], [3550], 0
        self.assertRaises(ValueError, utils.split_data, *args)


class TestFilterFactory(TestCase):
    """Tests for utils.classification_filter_factory"""

    def test_is_factory_func(self):
        """Test the returned value is a function"""

        dummy_arg = []
        returned_obj = utils.classification_filter_factory(dummy_arg)
        self.assertTrue(callable(returned_obj), 'Returned object not callable')

    def test_no_classification(self):
        """Test handling of tables without a classification in the metadata"""

        no_classification_table = Table()
        filter_func = utils.classification_filter_factory(['class1'])
        self.assertTrue(
            filter_func(no_classification_table),
            "Did not return true for table with no classification")

    def test_filtering(self):
        """Test the returned filter function correctly filters data table"""

        class1_table, class2_table = Table(), Table()
        class1_table.meta['classification'] = 'class1'
        class2_table.meta['classification'] = 'class2'

        filter_func = utils.classification_filter_factory(['class1'])
        self.assertTrue(
            filter_func(class1_table),
            "Returned False for desired classification")

        self.assertFalse(
            filter_func(class2_table),
            "Returned True for un-desired classification")
