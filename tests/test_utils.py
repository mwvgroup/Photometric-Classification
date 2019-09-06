#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``utils`` module."""

import time
from unittest import TestCase

import numpy as np
import sncosmo
import sndata
from astropy.table import Column, Table
from sndata.csp import dr3

from analysis_pipeline import utils


class TestTimeout(TestCase):
    """Tests for utils.timeout"""

    @staticmethod
    def sleep(seconds):
        """time.sleep wrapped with utils.timeout for a given number of seconds

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

        except:
            duration = np.round(time.time() - start, 1)
            self.assertEqual(seconds, duration)

    def test_error_for_float_arg(self):
        """Test an error is raised by the context manager when passed a float

        The signal manager expects integers, so if an error is not raised for
        a float then make sure you know what is actually going on.
        """

        self.assertRaises(TypeError, self.sleep, 2.1)


# Todo: Write a test that checks we get the expected value for chisq
class TestCalcModelChisq(TestCase):
    """Tests for analysis_pipeline.lc_fitting.calc_chisq"""

    @property
    def test_model(self):
        return sncosmo.Model('salt2')

    @property
    def test_data(self):
        # Type cast the 'band' column to U15 so we can use longer band names
        data = sncosmo.load_example_data()
        data['band'] = Column(data['band'], dtype='U15')
        return data

    # noinspection PyTypeChecker
    def test_arg_mutation(self):
        """Test that the input model and data table are not mutated"""

        model = self.test_model
        data = self.test_data
        utils.calc_model_chisq(data, model)

        self.assertTrue(all(data == self.test_data), 'Data table was mutated')
        self.assertSequenceEqual(
            list(self.test_model.parameters),
            list(model.parameters),
            'Model parameters were mutated')

    def test_out_of_range(self):
        """Test that out of range phase values and bands are dropped"""

        # The sncosmo docs assure us that the example data is within
        # range of the salt2 model
        data = self.test_data
        expected_chisq, expected_dof = utils.calc_model_chisq(data,
                                                              self.test_model)

        # Add values that are out of the model's phase range
        # Columns: 'time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'
        data.add_row([100000000, 'sdssu', 0, 0, 25, 'ab'])

        # Add values that are out of the model's wavelength range
        # Here we use the H band from CSP
        dr3.register_filters(force=True)
        data.add_row([1, 'csp_dr3_H', 0, 0, 25, 'ab'])

        returned_chisq, returned_dof = utils.calc_model_chisq(data,
                                                              self.test_model)

        self.assertEqual(
            expected_chisq, returned_chisq, 'Incorrect chisq value')
        self.assertEqual(
            expected_dof, returned_dof, 'Incorrect number of data points')

    def test_unregistered_band(self):
        """Test unregistered band names in the data table cause an error
        instead of being dropped.
        """

        # Add values with an unregistered filter
        data = self.test_data
        data['band'][0] = 'made up band'
        self.assertRaises(
            Exception, utils.calc_model_chisq, data, self.test_model)

    def test_unsorted_times(self):
        """Test the function return is independent of the time order"""

        data = self.test_data
        data.sort('time')
        initial_chisq = utils.calc_model_chisq(data, self.test_model)

        data[0], data[-1] = data[-1], data[0]
        new_chisq = utils.calc_model_chisq(data, self.test_model)
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

    @property
    def test_data(self):
        return sncosmo.load_example_data()

    @classmethod
    def setUpClass(cls):
        """Define test data and expected return values"""

        cls.band_names = sndata.sdss.sako18.band_names
        cls.lambda_effective = sndata.sdss.sako18.lambda_effective

        # Define what bands should be rest frame blue / rest frame red at
        # different redshift values
        cls.test_conditions = {
            0: {'blue': ['sdssu', 'sdssg'], 'red': ['sdssr', 'sdssi']},
            .5: {'blue': ['sdssu'], 'red': ['sdssg', 'sdssr', 'sdssi']},
            1: {'blue': [], 'red': ['sdssu', 'sdssg', 'sdssr', 'sdssi']}
        }

    def _test_bands_for_redshift(self, redshift, blue_data, red_data):
        """Assert whether correct tables were returned for a given redshift

        Args:
            redshift  (float): The redshift value
            blue_data (Table): The table of returned blue data
            red_data  (Table): The table of returned red data
        """

        blue_bands = set(blue_data['bands'])
        red_bands = set(red_data['bands'])

        expected_blue = self.test_conditions[redshift]['blue']
        expected_red = self.test_conditions[redshift]['red']

        err_msg = f'Wrongs bands for z={redshift}'
        self.assertCountEqual(expected_blue, blue_bands, err_msg)
        self.assertCountEqual(expected_red, red_bands, err_msg)

    def test_maintaines_metadata(self):
        """Test whether passed and returned tables have same metadata"""

        test_data = self.test_data
        test_data.meta['redshift'] = 0
        test_data.meta['test_key'] = 12345
        blue_data, red_data = utils.split_data(
            test_data, self.band_names, self.lambda_effective)

        self.assertEqual(
            blue_data.meta['test_key'], test_data.meta['test_key'])

        self.assertEqual(
            red_data.meta['test_key'], test_data.meta['test_key'])

    def test_redshift_0(self):
        """Test whether correct tables were returned for z=0"""

        test_data = self.test_data
        test_data.meta['redshift'] = 0
        blue_data, red_data = utils.split_data(
            test_data, self.band_names, self.lambda_effective)

        self._test_bands_for_redshift(1, blue_data, red_data)

    def test_redshift_point_5(self):
        """Test whether correct tables were returned for z=.1"""

        test_data = self.test_data
        test_data.meta['redshift'] = .5
        blue_data, red_data = utils.split_data(
            test_data, self.band_names, self.lambda_effective)

        self._test_bands_for_redshift(.5, blue_data, red_data)

    def test_redshift_1(self):
        """Test whether correct tables were returned for z=1"""

        test_data = self.test_data
        test_data.meta['redshift'] = 1
        blue_data, red_data = utils.split_data(
            test_data, self.band_names, self.lambda_effective)

        self._test_bands_for_redshift(1, blue_data, red_data)


class TestFilterFactory(TestCase):
    """Tests for utils.classification_filter_factory"""

    def test_is_factory_func(self):
        """Test the returned value is a function"""

        dummy_arg = []
        returned_obj = utils.classification_filter_factory(dummy_arg)
        self.assertTrue(callable(returned_obj), 'Returned object not callable')

    def test_filtering(self):
        """Test the returned filter function correctly filters data table"""

        no_classification_table = Table()
        class1_table = Table()
        class1_table.meta['classification'] = 'class1'

        class2_table = Table()
        class2_table.meta['classification'] = 'class2'

        filter_func = utils.classification_filter_factory(['class1'])

        self.assertTrue(
            filter_func(no_classification_table),
            "Did not return true for table with no classification")

        self.assertTrue(
            filter_func(class1_table),
            "Returned False for desired classification")

        self.assertFalse(
            filter_func(class2_table),
            "Returned True for un-desired classification")
