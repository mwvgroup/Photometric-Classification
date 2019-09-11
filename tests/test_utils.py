#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Tests for the ``utils`` module."""

import time
from unittest import TestCase

import numpy as np
import sncosmo
from astropy.table import Column, Table
from sndata.csp import dr3

from phot_class import utils


class TestTimeout(TestCase):
    """Tests for the utils.timeout context manager"""

    @staticmethod
    def run_timeout(seconds):
        """Run time.run_timeout wrapped with utils.timeout for given seconds

        Sleep for one second longer than the timeout context manager.

        Args:
            seconds (int): Number of seconds to run_timeout for
        """

        with utils.timeout(seconds):
            time.sleep(seconds + 1)

    def test_raises_timeout_error(self):
        """Test a TimeoutError is raised"""

        self.assertRaises(TimeoutError, self.run_timeout, 2)

    def test_context_duration(self):
        """Test an error is raised after the correct amount of time"""

        seconds = 1
        start = time.time()
        try:
            self.run_timeout(seconds)

        except TimeoutError:
            duration = np.round(time.time() - start, 1)
            self.assertEqual(seconds, duration)

    def test_error_for_float_arg(self):
        """Test an error is raised when the context manager is passed a float

        The signal manager expects integers, so if a TypeError is not raised
        for a float then make sure you know what is actually going on.
        """

        self.assertRaises(TypeError, self.run_timeout, 2.1)


class TestCalcModelChisq(TestCase):
    """Tests for utils.calc_model_chisq"""

    # We define data, model, and result properties to ensure a consistent,
    # un-mutated set of arguments when testing

    @property
    def model(self):
        """The model to use for testing"""

        model = sncosmo.Model('salt2')
        model.update(self.data.meta)
        return model

    @property
    def data(self):
        """The data to use for testing"""

        # Type cast the 'band' column to U15 so we can use longer band names
        data = sncosmo.load_example_data()
        data['band'] = Column(data['band'], dtype='U15')
        return data

    @property
    def results(self):
        """The results of fitting self.data with self.model"""

        vparams = ['t0', 'x0', 'x1', 'c']
        results, fitted_model = sncosmo.fit_lc(self.data, self.model, vparams)
        return results

    # noinspection PyTypeChecker
    def test_arg_mutation(self):
        """Test the model, data table, and results object are not mutated"""

        model = self.model
        data = self.data
        results = self.results

        utils.calc_model_chisq(data, results, model)
        self.assertTrue(all(self.data == data), 'Data table was mutated')
        self.assertListEqual(
            list(self.results.parameters),
            list(results.parameters),
            'Results were mutated')

        self.assertListEqual(
            list(self.model.parameters),
            list(model.parameters),
            'Model parameters were mutated')

    def test_out_of_range_data(self):
        """Test that out of range phase values and bands are dropped"""

        data = self.data
        expected_chisq, expected_dof = utils.calc_model_chisq(
            data, self.results, self.model)

        # Add values that are out of the model's phase range
        # Columns: 'time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys'
        data.add_row([100000000, 'sdssu', 0, 0, 25, 'ab'])

        # Add values that are out of the model's wavelength range
        # Here we use the H band from CSP
        dr3.download_module_data()
        dr3.register_filters(force=True)
        data.add_row([1, 'csp_dr3_H', 0, 0, 25, 'ab'])

        chisq, dof = utils.calc_model_chisq(data, self.results, self.model)
        self.assertEqual(expected_chisq, chisq, 'Chisq values are not equal')
        self.assertEqual(expected_dof, dof, 'Degrees of freedom are not equal')

    def test_unregistered_bands(self):
        """Test unregistered band names in the data table cause an error
        instead of being dropped.
        """

        # Add values with an unregistered filter
        data = self.data
        data['band'][0] = 'made up band'

        func = utils.calc_model_chisq
        args = (data, self.model, self.results)
        self.assertRaises(Exception, func, *args)

    def test_empty_table(self):
        """Test an error is raised the data table is empty"""

        col_names = self.data.colnames
        empty_table = Table(names=col_names)
        args = empty_table, self.results, self.model
        self.assertRaises(ValueError, utils.calc_model_chisq, *args)

    def test_correct_chisq(self):
        """Test the correct chisq and dof are returned for simulated data"""

        # Define necessary information to simulate a table of flux data
        model = self.model
        t0 = model.parameters[1]
        zp_system = 'ab'
        zero_point = 25
        band_name = 'sdssg'
        flux_offset = 100  # Difference between "observed" and model flux
        flux_err_coeff = .1  # flux error / flux

        # Create table of simulated flux
        phase = np.arange(-10, 30) + t0
        model_flux = model.bandflux(band_name, phase, zero_point, zp_system)
        flux = model_flux + flux_offset
        flux_err = flux * flux_err_coeff
        band = np.full(len(phase), band_name)
        zpsys = np.full(len(phase), zp_system)
        zp = np.full(len(phase), zero_point)

        data = Table(
            [phase, band, flux, flux_err, zp, zpsys],
            names=['time', 'band', 'flux', 'fluxerr', 'zp', 'zpsys']
        )

        result = sncosmo.utils.Result(
            {'vparam_names': ['t0', 'x0', 'x1', 'c']})

        expected_chisq = np.sum(((flux - model_flux) / flux_err) ** 2)
        chisq, dof = utils.calc_model_chisq(data, result, model)
        self.assertEqual(expected_chisq, chisq)


class TestSplitBands(TestCase):
    """Tests for utils.split_bands"""

    def runTest(self):
        """Test dummy bands are correctly separated into blue and red"""

        # Define dummy bands and their effective wavelengths
        band_names = ['blue1', 'red1', 'red2']
        lambda_eff = [4000, 5500, 6000]
        expected_blue_bands = band_names[:1]
        expected_red_bands = band_names[1:]

        blue_bands, red_bands = utils.split_bands(band_names, lambda_eff)
        self.assertListEqual(expected_blue_bands, blue_bands.tolist())
        self.assertListEqual(expected_red_bands, red_bands.tolist())


# Todo: Test cutoff wavelength
class TestSplitData(TestCase):
    """Tests for utils.split_data"""

    def test_redshift_dependency(self):
        """Assert whether correct bands were returned for a given redshift

        Tested for redshifts 0, .18, .55, and .78. Runs test with the cutoff
        wavelength set to infinity.
        """

        band_names = np.array(['u', 'g', 'r', 'i', 'z'])
        lambda_eff = np.array([3550, 4680, 6160, 7480, 8930])

        for redshift in (0, .18, .55, .78):
            # Determine expected blue and red wavelengths
            rest_frame_cutoff = 5500 * (1 + redshift)
            expected_blue = band_names[lambda_eff < rest_frame_cutoff].tolist()
            expected_red = band_names[lambda_eff > rest_frame_cutoff].tolist()

            # Split data according to the given redshift
            data = Table([band_names], names=['band'])
            blue_table, red_table = utils.split_data(
                data, band_names, lambda_eff, redshift, cutoff=float('inf'))

            err_msg = f'Wrongs bands for z={redshift}'
            self.assertListEqual(expected_blue, list(blue_table['band']),
                                 err_msg)
            self.assertListEqual(expected_red, list(red_table['band']),
                                 err_msg)

    def test_maintains_metadata(self):
        """Test whether passed and returned tables have same metadata"""

        test_data = Table([['u', 'g']], names=['band'])
        test_data.meta['dummy_key'] = 12345
        blue_data, red_data = utils.split_data(
            data_table=test_data,
            band_names=['u', 'g'],
            lambda_eff=[3550, 4680],
            z=0)

        self.assertIn('dummy_key', blue_data.meta, 'Blue table missing metadata')
        self.assertIn('dummy_key', red_data.meta, 'Red table missing metadata')

    def test_missing_effective_wavelength(self):
        """Test a value error is raised for missing effective wavelengths"""

        test_data = Table([['u', 'g']], names=['band'])
        args = (test_data, ['u'], [3550], 0)
        self.assertRaises(ValueError, utils.split_data, *args)


class TestFilterFactory(TestCase):
    """Tests for utils.classification_filter_factory"""

    def test_return_is_callable(self):
        """Test the returned value is a function"""

        dummy_arg = []
        returned_obj = utils.classification_filter_factory(dummy_arg)
        self.assertTrue(callable(returned_obj), 'Returned object not callable')

    def test_no_metadata_classification(self):
        """Test handling of tables without a classification in the metadata

        The returned filter function relies on the "classification"
        value from the meta data to filter out targets of different classes.
        Test the returned filter function returns ``True`` for a table
        without a classification in the meta data.
        """

        no_classification_table = Table()
        filter_func = utils.classification_filter_factory(['class1'])
        self.assertTrue(
            filter_func(no_classification_table),
            "Did not return true for table with no fitting")

    def test_filtering(self):
        """Test the returned filter function correctly filters data tables"""

        class1_table, class2_table = Table(), Table()
        class1_table.meta['fitting'] = 'class1'
        class2_table.meta['fitting'] = 'class2'

        filter_func = utils.classification_filter_factory(['class1'])
        self.assertTrue(
            filter_func(class1_table),
            "Returned False for desired fitting")

        self.assertFalse(
            filter_func(class2_table),
            "Returned True for un-desired fitting")
