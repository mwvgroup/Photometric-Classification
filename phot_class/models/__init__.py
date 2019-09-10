#!/usr/bin/env python3.7
# -*- coding: UTF-8 -*-

"""The ``models`` module defines ``Source`` classes for modeling
91bg-like supernovae with ``sncosmo``. The model is based on the 91bg template
from Nugent et al. 2002 but is extended into the ultra-violet. The model we
use here was originally intended for use with the FORTRAN package ``SNANA``.
Care was taken to ensure the model was ported correctly into Python and that
the predicted fluxes, parameter covariances, etc. are the same. To validate
our efforts, the ``snana_sims`` and ``sncosmo_sims`` submodules simulate
fluxes using the original and ported models respectively.

.. note::
  - For more information on the Nugent template see
    `Nugent et al. 2002 <https://iopscience.iop.org/article/10.1086/341707>`_.
  - For more information on sncosmo and ``Source`` classes see the
    `sncosmo documentation <https://sncosmo.readthedocs.io/>`_.
  - You will need to to have the ``SNANA`` package install on your
    machine in order to use the ``snana_sims`` package. For more information on
    SNANA see the `SNANA documentation <http://snana.uchicago.edu>`_.

About the 91bg model
--------------------

The is based on 35 SED templates provided by S. Gonzalez-Gaitan. These SEDs
are based on the Nugent et al. 2002 templates but have been extended into to
the UV. They form a 7 by 5 grid covering 7 different stretch and 5 different
color values. The ranges and relations of color and stretch were obtained by
using `SiFTO <https://iopscience.iop.org/article/10.1086/588518/meta>`_ to fit
template to multiple 91bg light-curves at low-z.

Available Model Versions
------------------------

There are multiple versions of the 91bg model included in this package. For
information on how each version determines the flux for a given set of
parameters, see the documentation for the given source class. It should **NOT**
be assumed that flux is predicted the same way for different model versions.

+--------+---------------------------+--------------------------------------------------------+
| Name   | Version                   | Description                                            |
+========+===========================+========================================================+
| sn91bg | 'phase_limited' (Default) | 1991bg model restricted to a phase similar to Salt 2.4 |
+--------+---------------------------+--------------------------------------------------------+
| sn91bg | 'color_interpolation'     | Full 1991bg model that interpolates in color space     |
+--------+---------------------------+--------------------------------------------------------+

.. warning:: The ``color_interpolation`` version of the 91bg model does not
   pass our test suite and should be used with caution.

Usage Example
-------------

>>> import sncosmo
>>> from matplotlib import pyplot as plt
>>>
>>> from phot_class import models
>>>
>>> # Make sncosmo aware of the 1991bg models
>>> models.register_sources()
>>>
>>> # Initialize a model where the version is the model mass
>>> source = sncosmo.get_source('sn91bg', version='phase_limited')
>>> model = sncosmo.Model(source=source)
>>>
>>> # Get more information on how the source determines flux
>>> print(help(source))
>>>
>>> # run the fit
>>> data = sncosmo.load_example_data()
>>> result, fitted_model = sncosmo.fit_lc(
>>>     data, model,
>>>     ['z', 't0', 'x0'],  # parameters of model to vary
>>>     bounds={'z':(0.3, 0.7)})  # bounds on parameters (if any)
>>>
>>> # Plot results
>>> fig = sncosmo.plot_lc(data, model=fitted_model, errors=result.errors)
>>> plt.show()

Module Level Functions
----------------------
"""

from ._sources import SN91bg as _SN91bg
from ._sources import load_template


def register_sources(force=False):
    """Register SN 1991bg-like models with SNCosmo

    Versions include: 'phase_limited', 'color_interpolation'

    Args:
        force (bool): Whether to overwrite an already registered source
    """

    import sncosmo

    sncosmo.register_loader(
        data_class=sncosmo.Source,
        name='sn91bg',
        func=_SN91bg,
        version='phase_limited',
        force=force)

    sncosmo.register_loader(
        data_class=sncosmo.Source,
        name='sn91bg',
        func=_SN91bg,
        version='color_interpolation',
        force=force)