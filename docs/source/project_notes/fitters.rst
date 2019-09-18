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
originally implemented with SiFTO. We here outline the approach of different
fitters relevant to our project and how they work.



SiFTO (`Conley+ 2008 <https://doi.org/10.1086/588518>`_)
--------------------------------------------------------

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

