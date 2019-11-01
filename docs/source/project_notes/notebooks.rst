Notebooks
=========

Although the ``phot_class`` package houses the core logic of our analysis, the
visualization and interaction with results from the package is performed using
Jupyter Notebooks. Descriptions of each notebook and what they inspect are
provided below. Online interactive versions are provided for each notebook via
`BinderHub`_. To launch a BinderHub server click `here`_.

.. note:: Please note that the BinderHub server might take a while to
   initialize and may require patience. Refreshing the page may help somewhat.

+------------------------------------+------------------------------------------------------------------------------+
| Notebook                           | Description                                                                  |
+====================================+==============================================================================+
|  apply_classification.ipynb        | Applies photometric classifications to all supernovae in our data sample.    |
+------------------------------------+------------------------------------------------------------------------------+
|  classifying_single_target.ipynb   | Demonstrates the classification technique of Gonzalez-Gaitan+ 14 on          |
|                                    | a single supernova.                                                          |
+------------------------------------+------------------------------------------------------------------------------+
| creating_config_files.ipynb        | Creates config files for CSP, DES, and SDSS.                                 |
+------------------------------------+------------------------------------------------------------------------------+
| fit_inspection.ipynb               | Inspects fit results for individual light curves.                            |
+------------------------------------+------------------------------------------------------------------------------+
| fitting_method_comparison.ipynb    | Comparison of classifcation results when using band-by-band vs. collective   |
|                                    | fitting.                                                                     |
+------------------------------------+------------------------------------------------------------------------------+
| iminuit_vs_emcee.ipynb             | A simple comparison of the fit_lc and mcmc_lc minimization routines.         |
+------------------------------------+------------------------------------------------------------------------------+
| inspecting_91bg_model.ipynb        | Demonstrates the properties of the 91bg model we use for classification.     |
+------------------------------------+------------------------------------------------------------------------------+
| salt2_fit_results.ipynb            | Minimal investigation of fit results from the Salt2 model.                   |
+------------------------------------+------------------------------------------------------------------------------+
| sdss_redshift_distribution.ipynb   | Plots of redshift distributions for the SNe data set.                        |
+------------------------------------+------------------------------------------------------------------------------+
| sncosmo_chisq_bug.ipynb            | Outlines a bug in the calculation of chi-squared in SNCosmo and demonstrates |
|                                    | that our results do not suffer from this bug.                                |
+------------------------------------+------------------------------------------------------------------------------+

.. _BinderHub: https://binderhub.readthedocs.io/en/latest/
.. _here: https://mybinder.org/v2/gh/mwvgroup/Photometric-Classification/master?filepath=notebooks%2F
