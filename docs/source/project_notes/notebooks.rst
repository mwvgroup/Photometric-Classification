Notebooks
=========

Although the ``phot_class`` package houses the core logic of our analysis, the
visualization and interaction with results from the package is performed using
Jupyter Notebooks. Descriptions of each notebook and what they inspect are
provided below. Online interactive versions are provided for each notebook via
`BinderHub`_. To launch a BinderHub server click `here`_ or select a notebook
from the table below to load that notebook directly.

.. note:: Please note that the BinderHub server might take a while to
   initialize and may require patience. Refreshing the page may help somewhat.

+------------------------------------+------------------------------------------------------------------------------+
| Notebook                           | Description                                                                  |
+====================================+==============================================================================+
| `apply_classification.ipynb`_      | Applies photometric classifications to all supernovae in our data sample.    |
+------------------------------------+------------------------------------------------------------------------------+
| `classifying_single_target.ipynb`_ | Demonstrates the classification technique of Gonzalez-Gaitan+ 14 on          |
|                                    | a single supernova.                                                          |
+------------------------------------+------------------------------------------------------------------------------+
| `compare_fits_to_published.ipynb`_ | Compares fit results from our pipeline against published values.             |
+------------------------------------+------------------------------------------------------------------------------+
| `fit_inspection.ipynb`_            | Inspects fit results for individual light curves.                            |
+------------------------------------+------------------------------------------------------------------------------+
| `inspecting_91bg_model.ipynb`_     | Demonstrates the properties of the 91bg model we use for classification.     |
+------------------------------------+------------------------------------------------------------------------------+
| `sncosmo_chisq_bug.ipynb`_         | Outlines a bug in the calculation of chi-squared in SNCosmo and demonstrates |
|                                    | that our results do not suffer from this bug.                                |
+------------------------------------+------------------------------------------------------------------------------+

.. _BinderHub: https://binderhub.readthedocs.io/en/latest/
.. _here: https://mybinder.org/v2/gh/mwvgroup/Photometric-Classification/master?filepath=notebooks%2F
.. _apply_classification.ipynb: https://hub.gke.mybinder.org/user/mwvgroup-photom--classification-2y0bflpl/notebooks/notebooks/apply_classification.ipynb
.. _classifying_single_target.ipynb: https://hub.gke.mybinder.org/user/mwvgroup-photom--classification-2y0bflpl/notebooks/notebooks/classifying_single_target.ipynb
.. _compare_fits_to_published.ipynb: https://hub.gke.mybinder.org/user/mwvgroup-photom--classification-2y0bflpl/notebooks/notebooks/compare_fits_to_published.ipynb
.. _fit_inspection.ipynb: https://hub.gke.mybinder.org/user/mwvgroup-photom--classification-2y0bflpl/notebooks/notebooks/fit_inspection.ipynb
.. _inspecting_91bg_model.ipynb: https://hub.gke.mybinder.org/user/mwvgroup-photom--classification-2y0bflpl/notebooks/notebooks/inspecting_91bg_model.ipynb
.. _sncosmo_chisq_bug.ipynb: https://hub.gke.mybinder.org/user/mwvgroup-photom--classification-2y0bflpl/notebooks/notebooks/sncosmo_chisq_bug.ipynb
