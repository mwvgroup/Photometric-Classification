The Command Line Interface
==========================

The ``phot_class`` package can be accessed either from within a Python
environment or via the command line using *run_pipeline.py*. For interacting
within an environment, see the API documentation. We here outline how to use
the command line interface.

Creating a Config File
----------------------

.. note:: Pre-made config files are available in the project's repository.

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

Specifying Priors
^^^^^^^^^^^^^^^^^

Valid prior values for both models include ``z``, ``t0``, ``amplitude``,
``x1``, and ``mwebv``. The ``sn91bg`` model also has an additional ``c``
parameter. It is important to understand how these parameters are handled by
the pipeline:

  - When fitting all bands simultaneously, ``t0``, ``amplitude``,
    ``x1``, and ``c`` are always varied. ``z`` is only varied if it is not
    specified by the prior.
  - When fitting bands independently, only  ``amplitude``, ``x1``, and ``c``
    are varied. ``z`` and ``t0`` are fixed to the value determined when fitting
    all bands simultaneously.
  - ``mwebv`` is never varied in any fit, and is fixed to the given value. If no
    value is given, the default value is 0.

See the :ref:`lc-fitting` section for more information on how we fit targets.

Specifying Kwargs
^^^^^^^^^^^^^^^^^

Valid keyword arguments depend slightly on the fitting function that is being
used, but are overwhelmingly uniform across the available functions. Each
fitting function is a simple wrapper around an ``sncosmo`` minimization
routine (The wrapping is just to avoid an argument mutation bug). For more
information on the available fitting routines, see :ref:`fit-functions`.

Running the Analysis
--------------------

.. argparse::
   :filename: ../../run_pipeline.py
   :func: create_cli_parser
   :prog: run_pipeline.py
