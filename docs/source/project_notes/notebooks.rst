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
|  compare_fits_to_published.ipynb   | Compares fit results from our pipeline against published values.             |
+------------------------------------+------------------------------------------------------------------------------+
|  fit_inspection.ipynb              | Inspects fit results for individual light curves.                            |
+------------------------------------+------------------------------------------------------------------------------+
|  inspecting_91bg_model.ipynb       | Demonstrates the properties of the 91bg model we use for classification.     |
+------------------------------------+------------------------------------------------------------------------------+
|  redshift_distributions.ipynb      | Plots of redshift distributions for various data sets.                       |
+------------------------------------+------------------------------------------------------------------------------+
|  sncosmo_chisq_bug.ipynb           | Outlines a bug in the calculation of chi-squared in SNCosmo and demonstrates |
|                                    | that our results do not suffer from this bug.                                |
+------------------------------------+------------------------------------------------------------------------------+

.. _BinderHub: https://binderhub.readthedocs.io/en/latest/
.. _here: https://mybinder.org/v2/gh/mwvgroup/Photometric-Classification/master?filepath=notebooks%2F
