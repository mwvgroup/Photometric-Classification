The Command Line Interface
==========================

The ``phot_class`` package can be accessed either from within a Python
environment or via the command line using *run_pipeline.py*. For interacting
within an environment, see the API documentation. We here outline how to use
the command line interface.

Creating a Config File
----------------------

Using the command line interface requires creating a configuration file.
Config files are used to define priors and keyword arguments used when fitting
each light-curve. Files should be saved in yaml format and follow the
template below. Notice that priors and keyword arguments are specified for
both the custom ``hsiao_x1`` and ``sn91bg`` models defined in the project's
analysis pipeline.

.. code-block::

    hsiao_x1:
      <object id>:
        kwargs:
          <keyword arguments>

        priors:
          <Initial guess for fit parameters>

    sn91bg:
      <object id>:
        kwargs:
          <keyword arguments>

        priors:
          <Initial guess for fit parameters>

Pre-made config files are available in the project's repository.


Running the Analysis
--------------------

.. argparse::
   :filename: ../../run_pipeline.py
   :func: create_cli_parser
   :prog: run_pipeline.py
