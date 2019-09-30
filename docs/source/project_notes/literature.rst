Literature Search
=================

Notes on some *potentially* relevant papers to this project.

SN1991bg properties and peculiar classification
-----------------------------------------------

Comparative Analysis of Peculiar Type Ia 1991bg-like Supernovae Spectra (`Doull+ 2011 <https://ui.adsabs.harvard.edu/abs/2011PASP..123..765D/abstract>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Quantitative Classification of Type I Supernovae Using Spectroscopic Features at Maximum Brightness (`Fengwu+ 2006 <https://ui.adsabs.harvard.edu/abs/2017arXiv170702543S/abstract>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The evolution of SNe classification is based on unquantified spectroscopic
qualities near peak luminosity. As a quick recap on the evolution of Type I
classification:

 - Minkowski (1941): Introduces type I/II
 - Wheeler & Levreault (1985) and Elias et al. (1985): Introduces type Ib
 - Wheeler & Harkness (1990): Introduces type Ic

Fengwu+ presents a quantified set of classification criteria based on the
depths of spectroscopic features. The boundaries of each classification are
selected to represent the classifications of objects already available in the
literature. The chosen SNe sample includes:

 - 3 SNe 2002cx-like
 - 1 SN 1991T- like
 - 21 SNe 1991bg-like
 - 8 SNe 1999aa-like
 - the unclassified peculiar object SN2000cx
 - 12 SNe Ib
 - 19 SNe Ic
 - 4 SNe Ib/c

Classifications are based on the depths of the Si II (6150 A) and O I (7774 A)
features relative to the pseudo-continuum. The steps of determining these
values include 1. Smoothing the spectra to eliminate narrow features,
2. identifying the boundaries of the features, 3. removing non-intrinsic
features, and 4. measuring the feature depths.

The spectra are smoothed to eliminate narrower features using the prescription
of Savitzky & Golay 1964. In general, SNe should not have strong, intrinsic
narrow features due to their high ejecta velocity. This process helps increase
the size of the usable data set by improving the precision of measurements
performed on low SNR spectra that would otherwise have to be dropped. The
window size :math:`w` and polynomial order :math:`o` of the smoothing function
are determined using the average wavelength interval :math:`\Delta \lambda` as

..  math::

  w = s \\times int(\\frac{50}{\\Delta \\lambda}) + 1
  o = max(3, \\frac{w - 1}{2})

After smoothing the spectra, the boundaries for both features are determined by
eye. Identifying the feature boundaries can be difficult since the high ejecta
velocity (broad feature widths) causes some features to blend together
(e.g. Si II 6355 with He I 6678). THe paper argues that the sample size is
small enough, and the inspection duration is quick enough, that the
introduction of human bias is negligible so long as only one person is used.

Non-intrinsic features that overlap with the features we are interested in are
replaced with a psedo continuum (linear interpolation) between the
contaminating feature's start and end points. Since the paper is interesting
in feature depth and not width, the effects of this de-resolution are negligible.

Finally, a pseudo-continuum is adopted for the Si II and O I features in the
same way as for the contaminating features. If the spectrum is visually
determined to not be smooth enough within the feature, the spectrum withing
the feature is smoothed using a 9th degree polynomial. The line depth is then
calculated as:

.. math::

   a = max(1 - \frac{F_\\lambda}{F_{\\lambda, continuum}})

The paper finds that peculiar Type Ia's generally have shallower Si II 6355
lines. The same cannot be said for O I 7774, where the normal and combined
peculiar SNe follow a similar range and distribution. However, the 91bg and
99aa objects are distinguishable by O I. This indicates an intrinsic diversity
of O I optical depths in SNe Ia photospheres.

Although the paper struggles to confidently distinguishing the normal and
peculiar subsets, they are able to find significant differences between SNe
Ib and Ic using the ratio :math:`r = a(λ6150 A) / a(OIλ7774 A)`. The Ib and Ic
populations are entirely separated by a line near :math:`r=1`.

The concluded classification criteria is as follows:

 1. SNe Ia (including normal Ia, Ia-1991bg and Ia-1999aa): :math:`a(6150 A) > 0.35`
 2. SNe Ib: :math:`a(6150 A) > 0.35` and :math:`a(6150 A) / a(7774 A) > 1`
 3. SNe Ic (except for Ic-BL): :math:`a(6150 A)<0.35` and :math:`a(6150 A) / a(7774 A) < 1`


The Birth Rate of Subluminous and Overluminous Type Ia Supernovae (`meng+ 2011 <https://ui.adsabs.harvard.edu/abs/2011A%26A...525A.129M/abstract>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Comparative Direct Analysis of Type Ia Supernova Spectra. II. Maximum Light (`Branch+ 2006 <https://ui.adsabs.harvard.edu/abs/2006PASP..118..560B/abstract>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This paper identifies classifications of SNe Ia using spectrograph observations
taken within three days of maximum (the time cutoff is chosen such that
the spectral evolution will be minimal while ensuring there are a sufficient
number of targets to be scientifically interesting). To simplify the process
of feature comparison, spectra are first tilted by multiplying the flux by
:math:`\\lambda^\\alpha` where :math:`\\alpha` is chosen such that the peak
flux at 4600 and 6300 A are equal. The Equivalent widths are then plotted for
the feature at 5750 A vs the feature at 6100 A. After applying a nearest
neighbor algorithm, four groups emerged: shallow silicon, core-normal,
broad line, and cool (which includes SN 1991bg).

broad-line SNe Ia have absorption features at 6100 A absorptions that are
broader and deeper than core-normal SNe Ia. However, SNe in this category do
not appear to follow a simple one-dimensional sequence based on their distance
from the core-normal population.

The shallow silicon group are not (necessarily) very different from the core
normal group. Other than a narrower Si feature, they look remarkably similar.
The primary reason for the spectroscopic differences seems to be the lower
temperature, as indicated by low temperature ion signatures (e.g. Ti).
Otherwise, they have the same ions evident in their spectra, just at very
different optical depths. This aligns with their lower temperatures since "as
noted by Hatano et al. (2002) and Ho flich et al. (2002), there is a
temperature threshold below which, owing to abrupt changes in key ionization
ratios, line optical depths change abruptly (Hatano et al. 1999)."

The core-normal subgroup have a very high degree of similarity, suggesting
a standard, common physical mechanism involving no large inhomogeneities near
the characteristic photospheric velocity of 12,000 km/s.


A high peculiarity rate for Type Ia SNe (`Li+ 1999 <https://ui.adsabs.harvard.edu/abs/2000AIPC..522...91L/abstract>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Evidence for a Spectroscopic Sequence among Type 1a Supernovae  (`Nugent+ 1995 <https://ui.adsabs.harvard.edu/abs/1995ApJ...455L.147N/abstract>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Determining the Type, Redshift, and Age of a Supernova Spectrum  (`Blondin+ <https://ui.adsabs.harvard.edu/abs/2007ApJ...666.1024B/abstract>`_
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
