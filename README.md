# SDSS-Classification

  [![Build Status](https://travis-ci.com/mwvgroup/SDSS-Classification.svg?token=MKWwaqNeMpyaNQ2HGxM7&branch=master)](https://travis-ci.com/mwvgroup/SDSS-Classification)

This repository applies the photometric classification technique from González-Gaitan et al. 2014 to SDSS, DES, and CSP data to identify peculiar Type Ia Supernovae (SNe Ia).  Work on the accompanying paper can be found [here](https://github.com/mwvgroup/91bg_paper).

1. [Project Todo List](#todo)
1. [Running the Analysis Pipeline](#running-the-analysis-pipeline)
1. [About the 91bg model](#about-the-91bg-model)
1. [Notes on SDSS Data](#notes-on-the-sdss-ii-sn-survey-data-sako-et-al-2018)
1. [Notes on DES Data](#notes-on-the-des-year-3-cosmology-data-brout-sako-et-al-2019-and-brout-scolnic-et-al-2019)
1. [Notes on CSP Data](#notes-on-the-csp-data)



## Todo:

- The chi-squared values calculated by SNCosmo include data points outside the phase range of the model. This results in unreasonbly (and for our purposes incorrectly) large chi-squared values. We need to calculate them manually.
- 4 parameter Satl2.4 fits have been visually inspected for CSP. This needs to be repeated for other models / surveys.
- The `num_points_<band>` column in the fit summary tables contain the number of data points per observer frame band. This needs to be changed to the number of points in the rest frame so that we can ignore lightcurves with too few data points.
- Our ported 91bg model does not always work with `sncosmo.nest_lc`.
- Cache and manually edit fit priors for each target.



## Running the Analysis Pipeline

The analysis pipeline runs light-curve fits for a combination of surveys, models, and number of fit parameters. By default the pipeline uses the Salt2.4 model to run 4 and 5 parameter fits (with / without redshift). To run the pipeline:

```bash
python run_pipeline.py -a fitting_params.yml -s csp
```

The `fitting_params.yml` file specified survey specific boundaries for each fit. To fit a different model or to only run fits with a certain number of parameters:

```bash
python run_pipeline.py -a fitting_params.yml -s csp -n 4 -m sn_91bg
```

A full description of available arguments can be retrieved by running `python run_pipeline.py -h`. 



The pipeline will automatically perform nested sampling to determine the initial fit parameters for each model / survey and save the results to file. If results for a given light curve fit have already been cached, then the pipeline will skip the sampling process and use the cahced values. Sampled values are determined on an object / model basis. This means the same priors are used for the 4 and 5 parameter fits and for all band pass collections. The Jupyter notebook `3_fit_inspection.ipynb` can be used to inspect light curve fits and manually adjust priors for any targets with poor fits.




## About the 91bg model

- The [91bg model](https://github.com/mwvgroup/SDSS-Classification/blob/issues/3/port_model/sncosmo_91bgmodel/_sncosmo_91bgmodel.py) is based on 35 [SED templates](https://github.com/mwvgroup/SDSS-Classification/tree/issues/3/port_model/snana_sims/91BG_SED) provided by S. Gonzalez-Gaitan. These SEDs are based on [Nugent's](https://iopscience.iop.org/article/10.1086/341707) 91bg but extending a bit to the UV. They form a 7\*5 grid with 7 different stretches and 5 different colors. The ranges and relations of colors and stretch to generate the different SNANA templates were obtained by fitting with [SiFTO](https://iopscience.iop.org/article/10.1086/588518/meta) this template to all 91bg at low-z.

- The 7 SEDs with different stretches can be generated from one single SED by stretching the phase (or time) axis: flux = template(t, &lambda;) &rArr; flux = template(t / stretch, &lambda;).

- Our approach is picking 5 SEDs with same stretch (st = 0.65) and different colors (c = [0, 0.25, 0.5, 1.0]), using template(t / (stretch/0.65), &lambda;) to fit for stretch and linearly interpolating over 5 SEDs to fit for color.

- For now, this model works well on recovering stretches and colors of snana simulated light curves but the fitted redshifts are slightly larger than the simulated value. We'll keep trying to further optimize this model.

  


## Notes on the SDSS-II SN Survey Data (Sako et al. 2018)

- Observed during three-month campaigns in Fall 2005, 2006, and 2007 as part of an extension to the original SDSS.
- SN span redshift range 0.05 < z < 0.4 (Frieman et al. 2008). 
- 10,258 Included sources. 3225 Variable, 499 classified as SN Ia, and 86 as Core collapse.
- All objects were observed in two or more filters and visually inspected to ensure they were not artifacts.
- Classification schemes:
  - **Unknown:** The light curve was too sparse and/or noisy to make a useful classification.
  - **Variable:** The source was observed in more than one observing season (and hence is not a supernova).
  - **AGN:** An optical spectrum was identified as having features associated with an active galaxy.
  - **SNII, SNIbc (either Ib or Ic), SNIa:** 
    - SN classifications without prefixes are based on a spectrum (including a few non-SDSS spectra)
    - A prefix “p” indicates the redshift is unknown and that the identification was made with photometric data only
    - A prefix “z” indicates that a redshift is measured from its candidate host galaxy and the classification uses that redshift as a prior. 
  - **SN Ia?:** classification is based on a spectrum that suggests a SN Ia but is inconclusive.
- Some SN candidates have associated notes indicating candidates that may have peculiar features or candidates where the typing spectrum was obtained by other groups. Quoting the data release: "We did not search for these peculiar features in a systematic way, but we have noted the likely peculiar features that were found." The flags are given in Table 4 and are as follows:
  - **1**  SN typing based on spectra obtained by groups outside SDSS. The spectra used for typing are not included in the data release. 
  - **2**  Peculiar type Ia SN possibly similar to sn91bg 
  - **3**  Peculiar type Ia SN possibly similar to sn00cx 
  - **4**  Peculiar type Ia SN possibly similar to sn02ci 
  - **5**  Peculiar type Ia SN possibly similar to sn02cx 
- Light curve fits are performed using PSNID, SALT2, and MLCS2k2
  - SALT2 fits are performed both with and without providing a spectroscopic redshift (where available). When fitting for redshift a flat prior is used.
  - Only data with a photometric flag less than 1024 is used for the fit. A full overview of the flag values is provided in Holtzman et al. 2008.
  - Epochs earlier than 15 days or later than 45 days in the rest frame are ignored. 152 epochs in 105 SN were also ignored after being visually classified as outliers.
  - Published fits use the SALT2 model as implemented in SNANA version 10.31b and the spectral templates and color law reported in Guy et al. (2010, G10)
  - Fits are performed using the Doi et al. 2010 SDSS filter response curves. These are available online [here](http://www.ioa.s.u-tokyo.ac.jp/~doi/sdss/SDSSresponse.html).
- The magnitudes in the online data release use the SDSS inverse hyperbolic sine magnitude system. These differ from the standard AB system by an additive constant found in Table 7. The fluxes in the online files have already been converted to the AB system and are given in Micro-Janskies.



## Notes on the DES Year 3 Cosmology Data (Brout, Sako et al. 2019 and Brout, Scolnic et al. 2019) 

- Five year survey from 2013 - 2018 in *ugriz* within 0.017 < z < 0.849
-  Average cadence of 7 days per filter
- These results contain the first three years of DES from Sept. 2013 to Feb. 2016
  - Discovered ∼12,000 transients
  - ∼3,000 identified as likely SNe Ia based on their light curves and out of 533 targeted for spectroscopic classification 251 were confirmed (D’Andrea et al. 2018).
  - These papers include a low-z SN sample from CfA3, CfA4, and CSP1. It's not clear if the numbers above include the low-z sample.
  
- The 31.00 zero-point is for internal DES use. The ZP in the public data files is 27.5.



## Notes on the CSP data

- The official CSP webpage for DR3 is [here](https://csp.obs.carnegiescience.edu/news-items/csp-dr3-photometry-released).
- Spectroscopic classification can be found in [Folatelli et al. (2013)](https://arxiv.org/abs/1305.6997). Data tables are available from [Vizier](http://cdsarc.u-strasbg.fr/viz-bin/cat/J/ApJ/773/53)
