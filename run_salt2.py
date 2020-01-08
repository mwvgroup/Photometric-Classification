#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""Fit SDSS data with Salt2"""

from copy import deepcopy
from pathlib import Path

import sncosmo
import yaml
from astropy.table import Table, vstack
from sndata.sdss import sako18

from phot_class import utils
from phot_class.fit_func_wraps import mcmc_fit, simple_fit
from phot_class.fitting import create_empty_table, fit_results_to_dict


def run_salt2_fits(obj_id, data, model, fit_func, fit_prior, fit_kwargs):
    """Run light curve fits on a given target using the salt2 model

    Args:
        obj_id        (str): Id of the object being fitted
        data        (Table): Table of photometric data
        model       (Model): The model to fit to the data
        fit_func (callable): The minimization function to use
        fit_prior    (dict): Initial parameter values
        fit_kwargs   (dict): Arguments for sncosmo.mcmc_lc to use when fitting

    Returns:
       Fit results and the fitted model as a dictionary
    """

    # Protect against mutation
    fit_kwargs = deepcopy(fit_kwargs) if fit_kwargs else None
    model = deepcopy(model)

    # Set initial values fo mcmc
    model.update(fit_prior)
    vparams = set(model.param_names) - {'mwebv'}
    if 'z' in fit_prior:
        vparams -= {'z'}

    # Remove bounds for unvaried parameters to prevent mcmc error
    fit_kwargs['bounds'] = \
        {p: v for p, v in fit_kwargs.get('bounds', {}).items() if p in vparams}

    # Fit data
    fit_kwargs.setdefault('warn', False)
    result_all, fit_all = fit_func(data, model, vparams, **fit_kwargs)
    return fit_results_to_dict(data, obj_id, 'all', result_all, fit_all)


def tabulate_fit_results(data_iter, model, fit_func, config, out_path=None):
    """Tabulate fit results for a collection of data tables

    Args:
        data_iter    (iter): Iterable of photometric data for different SN
        model       (Model): The model to fit to the data
        fit_func (callable): The minimization function to use
        config       (dict): Specifies priors / kwargs for fitting each model
        out_path      (str): Optionally cache progressive results to file

    Returns:
       An astropy table with fit results
    """

    if out_path and Path(out_path).exists():
        out_table = Table.read(out_path)

    else:
        # Includes meta data
        out_table = create_empty_table(model.param_names)

    for data in data_iter:

        # Set model parameters for the current object
        obj_id = data.meta['obj_id']
        prior = config[obj_id]['priors']
        kwargs = config[obj_id]['kwargs']

        if obj_id in out_table['obj_id']:
            continue

        try:
            fit_results_table = run_salt2_fits(
                obj_id=obj_id,
                data=data,
                model=model,
                fit_func=fit_func,
                fit_kwargs=kwargs,
                fit_prior=prior
            )

            out_table = vstack([out_table, fit_results_table])

        except KeyboardInterrupt:
            raise

        except Exception as e:
            out_table.add_row({
                'obj_id': obj_id,
                'message': str(e).replace('\n', '')
            })

        if out_path:
            out_table.write(out_path, overwrite=True)

    return out_table


if __name__ == '__main__':

    # Define the model to fit
    dust = sncosmo.F99Dust()
    salt2_model = sncosmo.Model(
        source='salt2',
        effects=[dust],
        effect_names=['mw'],
        effect_frames=['obs'])

    # Create iterable over light-curve data
    filter_func = utils.classification_filter_factory(
        ['AGN', 'Variable']
    )
    sdss_data = sako18.iter_data(verbose=True, filter_func=filter_func)

    # Define priors for salt2by modifying hsiao_x1 parameters from
    # an existing config file. This allows us to use pre-tabulated values like
    # redshift, E(B-V), etc.
    config_path = Path('./config_files/sdss_config_ext.yml')
    with open(config_path) as config_file:
        sdss_config = yaml.load(config_file, Loader=yaml.Loader)['hsiao_x1']

    for obj_config in sdss_config.values():
        obj_config['kwargs']['bounds']['x1'] = [-5, 5]
        obj_config['kwargs']['bounds']['c'] = [-.5, .5]
        obj_config['kwargs']['bounds']['z'] = [0, 1]

    # Run fits
    sako18.register_filters()
    iminuit_results = tabulate_fit_results(
        data_iter=sdss_data,
        model=salt2_model,
        fit_func=simple_fit,
        out_path='./results/sdss_salt2_fits.ecsv',
        config=sdss_config
    )

    mcmc_results = tabulate_fit_results(
        data_iter=sdss_data,
        model=salt2_model,
        fit_func=mcmc_fit,
        out_path='./results/sdss_salt2_mcmc_fits.ecsv',
        config=sdss_config
    )
