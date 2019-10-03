Classification Scheme
=====================

We here discuss the classification we wish to apply, which is originally
presented in `González-Gaitán et al. 2014 <https://ui.adsabs.harvard.edu/abs/2014ApJ...795..142G/abstract>`_

Background
----------

Despite being surprisingly homogeneous, type Ia Supernovae (SNe Ia) are in fact
a heterogeneous collection of events who's variety is well parameterized by
empirical descriptions such as stretch and color. Although prevalent theories
attribute supernova explosions as being the thermonuclear runaway of an
exploding CO white dwarf, modeling efforts have seen little progress
in providing a definitive explanation. The goal of classification is to
provide a framework from which we can further build and explore theoretical
explanations.

Sub-types of Type Ia SNe include:

- SN 1991T-like (Filippenko+ 1992a; Phillips+ 1992; Maza+ 1994)
- SN 1991bg-like (Filippenko+ 1992b; Leibundgut+ 1993; Hamuy+ 1994)
- 2002cx-like SNe or “SNe Iax” (e.g. Li+ 2003; Foley+ 2013)
- super-Chandrasekhar mass candidates (e.g. Howell+ 2006; Scalzo+ 2012)
- SN 2000cx-like (Li+ 2001; Candia+ 2003; Silverman+ 2013a)
- SN 2006bt-like (Foley+ 2010b)
- SNe Ia with possible circumstellar material (CSM) interactions (Hamuy 2003; Dilday+ 2012; Silverman+ 2013b)
- Ca-rich SNe (Perets+ 2010, 2011b; Valenti+ 2013a; Kasliwal+ 2012)

.. note:: See `Gal-Yam 2017 <https://link.springer.com/referenceworkentry/10.1007/978-3-319-21846-5_35>`_
   for a comprehensive discussion.

We are primarily focused on the identification of SN 1991-bg like objects
(91bgs). These tend to be around 1.1 mag fainter and decline notably faster
than normal Ia's. Their redder colors and more prevalent Ti lines also
indicate that they are cooler. Since 91bgs are cooler, the recombination
of Fe III to Fe II also occurs sooner causing the peak brightness in the
redder bands to happen sooner than in the bluer bands.

González-Gaitán et al. 2014 (GG14) puts forth a classification technique by
which 91bg and a few other peculiar SNe can be identified (described below).
This approach is compared against existing approaches (SNID and GELATO2) and
found to be in good agreement.


The Original Data Set
---------------------

GG14 focuses on low-redshift SN samples (z < 0.1). In addition to a few
targets picked from the literature, the data is primarily taken from:

- The Caĺan/Tololo survey (Hamuy et al. 1996a)
- The Carnegie Supernova Project (CSP)
- The Center for Astrophysics (CfA) (Hicken et al. 2009, 2012)
- The Lick Observatory Supernova Search (Ganeshalingam et al. 2010)

No initial cuts are applied to the data.


The Approach
------------

91bg's can be photometrically distinguished from normal Ia's by the morphology
of their light-curve - particularly in then redder bands (as described above).
To identify 91bg like objects, the photometry for each target is first split
into the red and blue bands as separated by 5500 angstroms in the rest frame.
Both data sets are then fit using two light-curve models: one representing
normal Ia's and one representing 91bg's. Using the resulting chi-squared
values, targets are classified based on their position in the following phase
space:

.. math::

    x \def \chi^2_{blue}(Ia) - \chi^2_{blue}(91bg)
    y \def \chi^2_{red}(Ia) - \chi^2_{red}(91bg)

By construction of the above coordinates, we expect 91bg's to fall in the
upper right (first) quadrant while type Ia's should fall in the lower left
(third) quadrant.

In principle there is no physical reason why the quadrants should be separated
by intersecting lines at :math:`(0, 0)`. To determine where in this phase
space we should draw the boundaries separating each classification, we use
whatever boundaries are found to maximize the figure of merit (FOM)
which is given by:

.. math::

    FOM = \frac{N_{true}}{N_{total}} \times \frac{N_{true}}{N_{true} + N_{false}}

Where each term is defined as the following:

 - :math:`N_{total}` is the total number of objects with a given type (i.e. the "truth")
 - :math:`N_{true}` is the number object correctly classified as a given type
 - :math:`N_{false}` is the number of objects falsely classified a given type


Here *"a given type"* refers to 91bg-like.

The added work of fitting red and blue bands independently could be avoided
and a chi-squared for all bands could be used. Although this works, it does not
work as well since there may be some normal SNe with lower stretch or redder
colors. It is also possible to use a simple cut on the fitted color and/or
stretch values (using either the normal or 91bg template). However, this
sufferes from a similar problem where highly reddened, normal Ia's at low
stretch would look like Ia's.

Some Additional Details
-----------------------

- The normal Ia template used is the Hsiao+ 2007 template
- The 91bg template is the Nugent+ 2002 template
- Data cuts were implemented as follows:

  1. Exclude targets without observations in at least two filters, each with
     at least one data point between -15 and 0 days and one between 0 and 25 days.
  2. Data past 85 days was ignored
  3. If observations are available in duplicate filters (e.g. the same filter
     from different surveys) the filter with the most data points is used.


