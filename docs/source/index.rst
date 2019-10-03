Classification of Peculiar Supernovae
=====================================

.. |python| image:: https://img.shields.io/badge/python-3.7-success.svg
    :target: #

.. |travis| image:: https://travis-ci.com/mwvgroup/Photometric-Classification.svg?branch=master
    :target: https://travis-ci.com/mwvgroup/Photometric-Classification

.. |cover| image:: https://coveralls.io/repos/github/mwvgroup/Photometric-Classification/badge.svg?branch=master
    :target: https://coveralls.io/github/mwvgroup/Photometric-Classification?branch=master

.. |docs| image:: https://readthedocs.org/projects/photometric-classification/badge/?version=latest
    :target: https://photometric-classification.readthedocs.io/en/latest/?badge=latest

.. |binder| image:: https://mybinder.org/badge_logo.svg
    :target: https://mybinder.org/v2/gh/mwvgroup/Photometric-Classification/master?filepath=notebooks%2F

.. rst-class:: badges

   +-------------------------------------------+
   | |python| |travis| |cover| |docs| |binder| |
   +-------------------------------------------+

.. warning:: This is an ongoing project. Results and documentation may be
   subject to frequent change.

With the discovery of the accelerating expansion of the universe, type Ia
supernovae (SNe Ia) have been used increasingly to determine cosmological
parameters. The use of SNe Ia as cosmological probes relies on the fact that
SN Ia luminosities at the time of maximum are not only bright but also have
low intrinsic scatter. However, SNe are in-fact a highly heterogeneous
collection of objects spanning a diverse collection of subtypes.

One approach to this challenge is to identify SNe subtypes based on their
light-curve properties. In Gonzalez-Gaitan et al. 2014 (G14) a photometric
identification technique was introduced for discriminating SN 1991bg-like
objects in photometric samples. Using several low-redshift samples from the
literature, it was demonstrated that this method is not only capable of
identifying dim, fast-declining SNe, but can also identify other peculiar
transients such as SNe Iax-like, SN 2006bt-like, and super-Chandrasekhar
SNe Ia.

We here apply the same classification technique to a larger target sample with
two significant changes. The first is the use of a newer model for 91bg-like
SN that has been extended further into the near infra-red (NIR) and
ultraviolet (UV). By using this extended model we are able to apply the
classification to a larger, higher-redshift sample of SNe Ia. The second
change is the inclusion of improved covariance values between the model
parameters stretch and color.


Using this Documentation
------------------------

This documentation serves as a complete project writeup and is a complete
collection of the developers thoughts/decisions. The **Project Notes** section
documents this project from a scientific perspective. It it provided to ensure
reproducibility of the results and to clarify various design decisions.
The **API Reference** documents how to use the project's code base along with
various technical clarifications. Source code for this project
can be found online via `GitHub`_.

.. _GitHub: https://github.com/mwvgroup/sdss-classification/

Project Repository Structure
----------------------------

We provide summaries for key files in the project repository::

   project parent
   ├── config_files/
   ├── docs/
   ├── notebooks/
   ├── phot_class/
   ├── results/
   ├── tests/
   │
   ├── README.md
   ├── environment.yml
   ├── run_pipeline.py
   ├── run_pipeline.sh
   └── setup.py


+----------------------+------------------------------------------------------------------------------+
| File                 | Description                                                                  |
+======================+==============================================================================+
| *config_files/*      | Files specifying fitting arguments anf priors for different models/surveys.  |
+----------------------+------------------------------------------------------------------------------+
| *docs/*              | The project documentation source code (what you're reading right now).       |
+----------------------+------------------------------------------------------------------------------+
| *notebooks/*         | Notebooks that inspect results of the classification pipeline.               |
|                      | These contain mostly plotting code.                                          |
+----------------------+------------------------------------------------------------------------------+
| *phot_class/*        | The python package containing all of the logic for fitting and               |
|                      | classifying light curves.                                                    |
+----------------------+------------------------------------------------------------------------------+
| *results/*           | Classification results returned from the ``phot_class`` package.             |
|                      | This directory can be re-generated by running *run_pipeline.sh*              |
+----------------------+------------------------------------------------------------------------------+
| *tests/*             | Tests for the ``phot_class`` package.                                        |
+----------------------+------------------------------------------------------------------------------+
| *README.md*          | Github Landing document pointing readers to the online documentation.        |
+----------------------+------------------------------------------------------------------------------+
| *environment.yml*    | Installation requirements for ``phot_class``.                                |
+----------------------+------------------------------------------------------------------------------+
| *run_pipeline.py*    | A command line interface for running the ``phot_class`` package.             |
+----------------------+------------------------------------------------------------------------------+
| *run_pipeline.sh*    | Runs the command line interface for various combinations of arguments        |
+----------------------+------------------------------------------------------------------------------+
| *setup.py*           | Installation script for the ``phot_class`` package.                          |
+----------------------+------------------------------------------------------------------------------+

.. toctree::
   :hidden:
   :maxdepth: 1

   Overview<self>

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Project Notes:

   project_notes/classification_scheme
   project_notes/data
   project_notes/fitters
   project_notes/literature
   project_notes/notebooks

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Getting Running:

   getting_running/installation
   getting_running/quick_start


.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: API Reference:

   api_reference/overview
   api_reference/fit_funcs
   api_reference/classification
   api_reference/fom
   api_reference/utils
   api_reference/models
   api_reference/sncosmo_sims
