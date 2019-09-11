Installation
============

This project uses the ``classification_pipeline`` package which is a
proprietary analysis pipeline that is **not** available via a package
manager. Source code can be downloaded to your local machine from `GitHub`_
or by using ``git``:

.. code:: bash

   git clone https://github.com/mwvgroup/sdss-classification/

Dependencies for the ``classification_pipeline`` package in addition to the
package itself can then be installed by running the following from within the
project repository:

.. code:: bash

   pip install -r requirements.txt
   python setup.py install --user

.. _GitHub: https://github.com/mwvgroup/sdss-classification/
