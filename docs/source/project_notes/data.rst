Data
====

We primarily consider observations from the Sloan Digital Sky Survey (SDSS)
Supernova Survey. However, during our investigations we also occasionally
considered data from the Dark Energy Survey (DES) Year 3 Cosmology release
(SN3YR), and the third data release of the Carnegie Supernova Project
(CSP DR3). We note that just because a data set is considered here or
elsewhere in our preliminary exploration, it is not guaranteed to be used
in the published work.

We rely on the ``sndata`` python package for programmatic access to
each of these data sets. In principle, the pipeline we have built should
extend trivially to any data set ported into ``sndata``.
Some notes are included below for each data set. These are included for use
by developers, and should not be considered a formal reference.

.. note:: For more information on ``sndata`` see the
  `official documentation <https://sndata.readthedocs.io/en/latest/>`_

The SDSS-II SN Survey Data (`Sako et al. 2018 <https://iopscience.iop.org/article/10.1088/1538-3873/aab4e0/meta>`_)
-------------------------------------------------------------------------------------------------------------------

The SDSS-II SN Survey was implemented as an extension to the original SDSS and
comprised of three, three-month campaigns in Fall 2005, 2006, and 2007.
The data release contains 10,258 total sources: 1,988 objects classified as
SN Ia, 2,049 SN II, 3 SLSN, 78 objects classified as either SN Ib or Ic, and
4,131 sources that are either variable or AGNs. Additionally, there are 2,009
targets with light curves that were deemed too sparse or noisy to provide a
classification.

Out of the total target list, dedicated spectroscopic observations are only
available for 626 SN, resulting in 499 SN Ia (about 43%) being
spectroscopically confirmed. All objects were observed in two or more filters
and visually inspected to ensure they were not artifacts. SNe in the sample
were found to span a redshift range of 0.05 < z < 0.4 (Frieman et al. 2008).

Each target in the sample is assigned one of the following classifications:
    - **Unknown:** The light curve was too sparse and/or noisy to make a useful classification.
    - **Variable:** The source was observed in more than one observing season (and hence is not a supernova).
    - **AGN:** An optical spectrum was identified as having features associated with an active galaxy.
    - **SNII, SNIbc (either Ib or Ic), SNIa:**
      - SN classifications without prefixes are based on a spectrum (including a few non-SDSS spectra)
      - A prefix “p” indicates the redshift is unknown and that the identification was made with photometric data only
      - A prefix “z” indicates that a redshift is measured from its candidate host galaxy and the classification uses that redshift as a prior.

    - **SN Ia?:** classification is based on a spectrum that suggests a SN Ia but is inconclusive.

Some SN candidates have associated notes indicating targets that may have peculiar
features or candidates where the classification was obtained by other groups.
Although the survey did not perform a systematic search for peculiar features,
likely peculiar features were noted if found.

The photometric flags are given in Table 4 and are as follows:
    - **1** - SN typing based on spectra obtained by groups outside SDSS. The spectra used for typing are not included in the data release.
    - **2** - Peculiar type Ia SN possibly similar to sn91bg
    - **3** - Peculiar type Ia SN possibly similar to sn00cx
    - **4** - Peculiar type Ia SN possibly similar to sn02ci
    - **5** - Peculiar type Ia SN possibly similar to sn02cx

Light curve fits are performed using PSNID, SALT2, and MLCS2k2. Published fits
use the SALT2 model as implemented in SNANA version 10.31b and the spectral
templates and color law reported in Guy et al. (2010, G10). Fits were
performed using the Doi et al. 2010 SDSS filter response curves. These are
available online [here](http://www.ioa.s.u-tokyo.ac.jp/~doi/sdss/SDSSresponse.html).
The following cuts were applied to the photometric data.

  - Only data with a photometric flag less than 1024 is used for the fit.
    A full overview of the flag values is provided in Holtzman et al. 2008.
  - Epochs earlier than 15 days or later than 45 days in the rest frame are
    ignored. 152 epochs in 105 SN were also ignored after being visually
    classified as outliers.

.. important:: SDSS measures magnitudes using the inverse hyperbolic sine
   magnitude system. These differ from the standard AB system by an additive
   constant found in Table 7. The fluxes in the online files have already been
   converted to the AB system and are given in Micro-Janskies.
