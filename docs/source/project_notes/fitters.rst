Light-Curve Fitters
===================

Light-curve fitters are an extremely common tool used in supernova analyses,
particularly in cosmological studies. Perhaps not surprisingly, there exist
a multitude of available fitting routines each with its own unique goals,
approaches, and implementations. Some examples of this diversity include:

 - MLCS/MLCS2k2 (Riess+ 1996, Jha+ 2007)
 - stretch (Perlmutter+ 1997; Goldhaber+ 2001; Knop+ 2003)
 - super stretch (Wang+ 2006)
 - DM15 (Hamuy+ 1996b; Prieto+ 2006)
 - BATM (Tonry+ 2003),
 - CMAGIC (Wang+ 2003),
 - SALT/SALT2 (Guy+ 2005, 2007)

It is important to distinguish between the different parameters offered by
each light curve fitter, and more important still to acknowledge that even
when two fitters use parameters with the same name, those two parameters are
not, in fact, the same.

The photometric classification scheme that our project is based on was
originally implemented with SiFTO. We here outline the approach of various
third party fitters relevant to our project and how they work.


SiFTO (`Conley et al. 2008 <https://doi.org/10.1086/588518>`_)
--------------------------------------------------------------

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


SALT2 (`Guy et al. 2007 <https://www.aanda.org/htbin/resolve?bibcode=2007A%26A...466...11GFUL>`_)
-------------------------------------------------------------------------------------------------

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
