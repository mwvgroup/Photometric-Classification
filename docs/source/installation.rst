Installation
============

This project uses the ``phot_class`` package which is a
proprietary analysis pipeline that is **not** available via a package
manager. Source code can be downloaded to your local machine from `GitHub`_
or by using ``git``:

.. code:: bash

   git clone https://github.com/mwvgroup/Photometric-Classification

To ensure reproducibility and reduce the potential for programmatic errors,
it is recommended to work within a conda environment. The necessary
installation dependencies for creating this environment have been specified
in a configuration file for convenience:

.. code:: bash

   conda env create --name phot_class --file Photometric-Classification/environment.yml


The ``phot_class`` can then be installed as follows:

.. code:: bash

   # Activate the conda environment
   conda activate phot_class

   python Photometric-Classification/setup.py install --user

   # Optionally, you can then deactivate the conda environment when finished
   conda deactivate

.. _GitHub: https://github.com/mwvgroup/Photometric-Classification
