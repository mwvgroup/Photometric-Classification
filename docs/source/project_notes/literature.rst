Literature Search
=================

Notes on some *potentially* relevant papers to this project.

Properties of 91bg-like SNe
---------------------------

Comparative Analysis of Peculiar Type Ia 1991bg-like Supernovae Spectra (`Doull+ 2011 <https://ui.adsabs.harvard.edu/abs/2011PASP..123..765D/abstract>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: **TLDR:** Several 91bg-like spectra are fit with SYNOW and common
   ions are identified.

Several 91bg-like spectra are selected according to the Branch 2006 sub-typing
method and are fit with the SYNOW spectral fitter. These SNe have an unusually
deep and wide absorption around 4200 A and very strong Si II at 5972 A. Nine
Ions fit most of the observed features well at all epochs: O I, Na I, Mg II,
Si II, S II, Ca II, Ti II, Cr II, and Fe II. At early times Ca I also seems
to play a role.

The 4200 A absorption feature is well fit by a combination of Ti II and Mg II.
This is indicative of a cooler temperature. As the spectra evolve to around
20 days, the absorptions and emissions become more extreme, and the
photosphere velocity begins to decrees. For phases beyond this, assumptions
of the SYNOW fitter begin to play a noticeable effect.

Some of the spectral diversity may be due to differences in viewing angle on
an asymmetric explosion (see Maeda+ 2010). The spectral similarities indicate
SNe are a continuous distribution and it is possible the heterogeneity of SNe
Ia is not due to fundamental differences in the underlying physical mechanism.
Instead, there may be a set of primary, stable parameters for the explosion
followed by secondary sets of slightly less stable parameters.

This explanation would explain the structure we see in Branch style W-W plots.
The most stable parameters would result in the highest density regions
with the most common SNe Ia (the core-normals). The lowest density regions
would contain SNe Ia with less stable parameters. The SNe in between would be
transition regions.


The Birth Rate of Subluminous and Overluminous Type Ia Supernovae (`meng+ 2011 <https://ui.adsabs.harvard.edu/abs/2011A%26A...525A.129M/abstract>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

(`Li+ 2010 <>`_)
^^^^^^^^^^^^^^^^

Evidence for a Spectroscopic Sequence among Type Ia Supernovae  (`Nugent+ 1995 <https://ui.adsabs.harvard.edu/abs/1995ApJ...455L.147N/abstract>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Determining the Type, Redshift, and Age of a Supernova Spectrum  (`Blondin+ <https://ui.adsabs.harvard.edu/abs/2007ApJ...666.1024B/abstract>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


SNe Classification
------------------

Quantitative Classification of Type I Supernovae Using Spectroscopic Features at Maximum Brightness (`Fengwu+ 2006 <https://ui.adsabs.harvard.edu/abs/2017arXiv170702543S/abstract>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: **TLDR:** A quantitative definition of what feature strengths
   indicate a Ia, Ib, and Ic.

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

  w = s \times int(\frac{50}{\Delta \lambda}) + 1

  o = max(3, \frac{w - 1}{2})

After smoothing the spectra, the boundaries for both features are determined by
eye. Identifying the feature boundaries can be difficult since the high ejecta
velocity (broad feature widths) causes some features to blend together
(e.g. Si II 6355 with He I 6678). THe paper argues that the sample size is
small enough, and the inspection duration is quick enough, that the
introduction of human bias is negligible so long as only one person is used.

Non-intrinsic features that overlap with the features we are interested in are
replaced with a pseudo continuum (linear interpolation) between the
contaminating feature's start and end points. Since the paper is interesting
in feature depth and not width, the effects of this de-resolution are negligible.

Finally, a pseudo-continuum is adopted for the Si II and O I features in the
same way as for the contaminating features. If the spectrum is visually
determined to not be smooth enough within the feature, the spectrum withing
the feature is smoothed using a 9th degree polynomial. The line depth is then
calculated as:

.. math::

   a = max(1 - \frac{F_\lambda}{F_{\lambda, continuum}})

The paper finds that peculiar Type Ia's generally have shallower Si II 6355
lines. The same cannot be said for O I 7774, where the normal and combined
peculiar SNe follow a similar range and distribution. However, the 91bg and
99aa objects are distinguishable by O I. This indicates an intrinsic diversity
of O I optical depths in SNe Ia photospheres.

Although the paper struggles to confidently distinguishing the normal and
peculiar subsets, they are able to find significant differences between SNe
Ib and Ic using the ratio r = a(6150) / a(7774). The Ib and Ic
populations are entirely separated by a line near :math:`r=1`.

The concluded classification criteria is as follows:

 1. SNe Ia (including normal Ia, Ia-1991bg and Ia-1999aa): a(6150 A) > 0.35
 2. SNe Ib: a(6150 A) < 0.35 and a(6150) / a(7774) > 1
 3. SNe Ic (except for Ic-BL): a(6150)<0.35 and a(6150) / a(7774) < 1


Comparative Direct Analysis of Type Ia Supernova Spectra II. Maximum Light (`Branch+ 2006 <https://ui.adsabs.harvard.edu/abs/2006PASP..118..560B/abstract>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. note:: **TLDR:** SNe Ia are subclassed into shallow silicon, core-normal,
   broad line, and cool groups based on the strength 5750 A vs 6100 A.

This paper identifies classifications of SNe Ia using the width of the 5750
and 6100 features (usually attributed to Si ii at 5972 and 6355). To simplify
the process of feature comparison, spectra are first tilted by multiplying
the flux by :math:`\lambda^\alpha` where :math:`\alpha` is chosen such that
the peak flux near 4600 and 6300 A are equal. The Equivalent widths are then
plotted for the feature at 5750 A vs the feature at 6100 A. After applying a
nearest neighbor algorithm, four groups emerged: shallow silicon, core-normal,
broad line, and cool (which includes SN 1991bg).

Broad-line SNe Ia have absorption features at 6100 A absorptions that are
broader and deeper than core-normal SNe Ia. However, SNe in this category do
not appear to follow a simple one-dimensional sequence based on their distance
from the core-normal population.

The shallow silicon group are not (necessarily) very different from the core
normal group. Other than a narrower Si feature, they look remarkably similar.
The primary reason for the spectroscopic differences seems to be the lower
temperature, as indicated by low temperature ion signatures (e.g. Ti).
Otherwise, they have the same ions evident in their spectra, just at very
different optical depths. This aligns with their lower temperatures since "as
noted by Hatano et al. (2002) and Ho Flich et al. (2002), there is a
temperature threshold below which, owing to abrupt changes in key ionization
ratios, line optical depths change abruptly (Hatano et al. 1999)."

The core-normal subgroup have a very high degree of similarity, suggesting
a standard, common physical mechanism involving no large inhomogeneities near
the characteristic photosphere velocity of 12,000 km/s.


PELICAN: deeP architecture for the LIght Curve ANalysis (`Pasquet+ 2019 <https://ui.adsabs.harvard.edu/abs/2019A%26A...627A..21P/abstract>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


Previous works using SDSS
-------------------------

Line Profiles of Intermediate Redshift Type Ia Supernovae (`Konishi+ 2011 <https://ui.adsabs.harvard.edu/abs/2011arXiv1103.2497K/abstract>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Spectral properties of type Ia supernovae up to z ∼ 0.3 (`Nordin+ 2011a <https://ui.adsabs.harvard.edu/abs/2011A%26A...526A.119N/abstract>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Evidence For A Correlation Between The Si Ii Λ4000 Width And Type Ia Supernova Color (`Nordin+ 2011b <https://iopscience.iop.org/article/10.1088/0004-637X/734/1/42>`_)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


