.. _fitters:

Implementation
==============

We here discuss various implementation details related to how we fit and
classify light-curves.


Light-Curve Fitters
-------------------

Light-curve fitters are an extremely common tool used in supernova analyses,
particularly in cosmological studies. Perhaps not surprisingly, there exist
a multitude of available fitting routines each with its own unique goals,
approaches, and implementations. It is important to distinguish between the
different options offered by each light curve fitter.

The photometric classification scheme that our project is based on was
originally implemented with SiFTO. We have chosen to instead work with SNCosmo
so that we can use a different set of models. We here outline the approach of
both fitters and how they work.


SiFTO (`Conley et al. 2008 <https://doi.org/10.1086/588518>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. important:: SiFTO is not used at any stage of our analysis. We consider it
   here because understanding it is necessary for completeness.

SiFTO is an empirical light-curve fitter that is similar to SALT2 in the sense
that it uses magnitude, light curve shape (i.e. stretch), and color to fit
supernovae. However, it varies distinctly in approach from SALT2 by using
band-specific stretch and color parameters as opposed to a set of global,
light-curve specific parameters.

**Stretch:** The shape of the SiFTO model is described by a single stretch
parameter measured in the rest-frame B-band. However, the model applies this
stretch value differently in each observed filter. The functional form of this
application depends both on the B-band stretch value and the effective wavelength
of the bandpass.

**Color:** SiFTO does not incorporate a dedicated color term, but instead
allows the scale factor of the template to vary independently for each band.
The resulting colors can be combined to form a single color parameter afterwords,
however, this calculation is performed when using fit results to estimate
distance and not as part of the fit itself.

When trained on the same data sets, the results of SiFTO are consistent with
SALT2, although the scatter in the SiFTO results is lower. Perhaps not
surprisingly, when the band-specific SiFTO parameters are generalized to a
single "stretch" and "color" value, there is a strong correlation between the
two model parameters. This correlation is roughly linear for stretch, but is
better described by a third-order polynomial for color.

We note the following warning from Conley+ 2008:

.. warning:: "SiFTO is designed for use with modern, well-measured SN data
   sets, and is not suited for the analysis of poorly sampled light curves.
   Specifically, if the number of observations in a filter in the range -20 to
   40 rest-frame days relative to B maximum is less than about 3, then the fit
   in that filter may not be reliable. This proviso only applies to historical
   SN Ia samples – even at the highest redshifts, current observations are of
   high enough quality that this is not a problem. Our approach is not well
   suited to the rest-frame near-IR (I or redder), where a stretch-like
   prescription does not work well."


SNCosmo (`sncosmo.readthedocs.io <https://sncosmo.readthedocs.io>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``sncosmo`` is a Python package built for supernova cosmology. At it's core,
the package is designed so users can customize their analysis with purpose
built models, filters, minimization routines, etc.. There are a number of
built in minimization routines that can be used to fit light-curves. However,
all of the built in routines are designed to simultaneously fit multiple
band-passes - a very different approach to SiFTO.

There is no single set of parameters used by ``sncosmo``. Generally, the
parameters are determined by the model that is being fit. Additional
parameters can also be added to any model to represent external effects like
extinction. These effects can be implimented at either the rest or observed
frame.


Models
------

We consider multiple light curve models in our analysis, including SALT2,
a modified version of the Hsiao model, and a custom 91bg model.


SALT2 (`Guy et al. 2007 <https://www.aanda.org/htbin/resolve?bibcode=2007A%26A...466...11GFUL>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

SALT2 is intended to fit light curves while accounting for intrinsic variations
in the supernova population. The original SALT model used the spectral sequence
of Nugent et al. 2002 in conjunction with the parameters phase, wavelength, and
stretch. SALT2 takes a similar approach but with several improvements, among
which include the use of spectroscopic data in the training sample. In addition
to improving the model's resolution, this allows k-corrections to be handled in
a consistent way for both photometric and spectroscopic observations.

The SALT2 model is similar to a principal component decomposition times an
overall color law :math:`CL(\lambda)`.

.. math::

    F = x_0 \times [M_0(phase, \lambda) + x_1 \times M_1(phase, \lambda) +  ...] \times exp[c \times CL(\lambda)]

Here :math:`x_i` are the components of the model, :math:`c` is a color term,
:math:`M_0` represents the average spectral template, and the remaining
:math:`M_i` describe the variability of the supernova population. At the time
it's inception, the data available did not have a sufficiently high sampling
to perform a full principle component analysis. Furthermore, it was determined
during the training of the model that the terms for :math:`i>1` could not be
well constrained.

A total of 176 supernovae were used in the SALT2 training set, Eeach having
at least two light-curves in different filters. The initial training value for
:math:`M_0` was chosen to be the original SALT model with a stretch of 1.
:math:`M_1` was chosen as the difference between :math:`M_0` and SALT with a
stretch of 1.1.

Guy et al. 2007 outlines three potential use cases for SALT2 or its potential
successors:

1. The classification of type Ia supernovae via a chi-squared cutoff
2. The photometric estimation of redshifts
3. Inspection of the intrinsic variability of supernovae by improving the data
   set to allow more terms in the principle component analysis. In doing so one
   could look for correlations between supernova properties and parameter
   values.


Custom Hsiao-like Model (`Hsiao et al. 2007 <https://doi.org/10.1086/518232>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Hsiao model has a similar light-curve profile to that of SALT2, but has
several major differences. First, the model is trained on a more diverse set
of SNe than the SALT2 model. This makes it a better representation of
spectroscopically normal SNe Ia since it includes more of that subtype's
intrinsic diversity. It also covers a much broader phase and wavelength range,
allowing us to fit using a larger subset of the measured photometry.

The Hsiao model built in to ``sncosmo`` also only has three parameters: ``z``,
``t0``, and ``amplitude``. We add an additional "stretch" parameter called
``x1`` for convenience. This parameter is used to stretch the light-curve phase
as

.. math::

    stretched phase = phase / (1 - x1)

We intrinsically bound :math:`-.5 \leq x1 \leq .5` within the model.


Custom sn1991bg-like Model (`Nugent et al. 2002 <https://iopscience.iop.org/article/10.1086/341707>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The 91bg model used for this project is based on the 91bg template from
Nugent et al. 2002 but is extended into the ultra-violet. This model was
originally formatted for use with the FORTRAN package ``SNANA``.
Care was taken to ensure the model was ported correctly into Python and that
the predicted fluxes, parameter covariances, etc. are the same.

The model works by interpolating from a grid of SED templates covering 7
stretch and 5 color values. The ranges and relations of color and stretch were
obtained by using `SiFTO <https://iopscience.iop.org/article/10.1086/588518/meta>`_
to fit the template to multiple 91bg light-curves at low-z.

The full phase range of the 91bg template extends from -18 to 100 days. When
comparing this model against other models (e.g. via a chi-squared value) it is
sometimes beneficial to limit the phase range of our model to more closely
resemble what it is being compared to. For this reason the model has been
specifically programmed so that the template can be arbitrarily limited in
phase space at instantiation.


.. _lc-fitting:

Fitting and Classification
--------------------------

Our implementation of the classification method is as follows.

  0. The Milky Way extinction is determined for each target using the
     `Schlegel, Finkbeiner & Davis (1998) <https://doi.org/10.1086/305772>`_
     dust map and the `Fitzpatrik 99 <https://doi.org/10.1086/316293>`_
     extinction law. This value is never varied in any fit, and is fixed to
     the given value. We can also optionally set the extinction to zero.
  1. Each light curve is fit using both models and all available band passes.
     At this step ``t0``, ``amplitude``, ``x1``, and ``c`` are always varied.
     ``z`` is only varied if it is not specified by a prior (i.e. if it is
     not available from a spectroscopic observation). Results from this fit
     are used to determine the characteristic parameters for the given
     light-curve (The values one might publish in a summary table).
  2. Each bandpass is fit independently using both models. Here, ``z`` and `
     `t0`` are fixed to the value determined when fitting all bands
     simultaneously.
  3. Any fits that fail are dropped from our sample.
  4. The bandpasses are separated into the rest-frame blue and red
     (blue/redward of 5500 Angstroms.)
  5. The chi-squard values from the band-by-band fits are summed for each
     model in both the red and blue bandpasses. These values are used to
     determine the position of each target on the classification plot
     (see the :ref:`classification` section).
