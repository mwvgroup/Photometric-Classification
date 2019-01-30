- # SDSS-Classification

  [![Build Status](https://travis-ci.com/mwvgroup/SDSS-Classification.svg?token=MKWwaqNeMpyaNQ2HGxM7&branch=master)](https://travis-ci.com/mwvgroup/SDSS-Classification)

  This project applies the photometric classification technique from González-Gaitan 
  2014 on SDSS data to identify peculiar Type Ia Supernovae (SNe Ia). DES data is also
  considered, but is not a primary objective.

  ## Todo:

- 

- 

- SNCosmo will raise occasional warnings about poor fits, bad S/N, and dropped bands. These warnings are currently being logged but otherwise ignored and need to be handled correctly.

- If a redshift for a supernova is not available, SNCosmo is instructed to fit for the redshift. A prior needs to be provided for the redshift.

- The `num_points_<band>` column in the fit summary tables contain the number of data points per observer frame band. We are interested in the number of points in the rest frame.

- We are currently fitting data with SNCosmo's built-in 91bg model, but want to use our own custom model. This custom model is currently used to simulate SNANA light curves.

- Determine what K-correction was used in the published SDSS fit results and improve the quality of fits with large chi-squared compared to published results.

## Directory Overviews / File Lists

#### *data_access/* 

  A python 2.7 module for accessing SDSS and DES supernova data. Data is downloaded
  automatically if it is not locally available. An example of accessing SDSS data
  is provided below. Note that the DES interface is the same, except you would
  import `des_data` instead of `sdss_data`.

```python
  from data_access import sdss_data
  
  # Summary table of SDSS SN data
  print(sdss_data.master_table) 
  
  # Get data for a specific object
  print(sdss_data.get_data_for_id(685))
  
  # Iterable of SNCosmo input tables for each target
  for table in sdss_data.iter_sncosmo_input():
      print(table)
```

#### snana/ 

  Simulating 91bg light curves using the SNANA Fortran package

- **sim_91bg.input:** An input file to SNANA that generates simulation results
- **snana_results/:** Directory of simulation results generated by *sim_91bg.input*
- **SNANA_plots.ipynb:** Plots exploring 91bg simulation results.

  

#### SNCosmo/

  Fitting SDSS and DES data with the SNCosmo python package

- **fit_sdss.py:** Uses SNCosmo to fit SDSSS light curves with a normal sn Ia model and the builtin 91bg model in the `ug`, `riz`, and `ugriz` rest frame bands.
- **fit_des.py:** The same as *fit_sdss.py* except for DES data.
- **sncosmo_results.ipynb:** Compares fit results between the normal and 91bg models in different bands.
- **compare_to_published.ipynb:** Compares fit results with published results.

## Notes on the SDSS-II SN Survey Data (Seiko et al. 2018)

- Observed during three-month campains in Fall 2005, 2006, and 2007
- Redshift range 0.05 < z < 0.4 (Frieman et al. 2008)
- 10,258 Included sources. 3225 Variable, 499 classified as SN Ia, and 86 as Core collapse.
- All objects were observed in two or more filters and visually inspected for artifacts.
- *ugriz* Filters spaning 350 - 1000 nm (Fukugita et al. 1996)
- Classification schemes:
  - **Unknown:** The light curve was too sparse and/or noisy to make a useful classification
  - **Variable:** The source was observed in more than one observing season
  - **AGN:** An optical spectrum was identified as having features associated with an active galaxy
  - **SNII, SNIbc (either Ib or Ic), SNIa:** SN classifications without prefixes are based on a spectrum (including a few non-SDSS spectra)
  - A prefix “p” indicates the redshift is unknown and that the identification was made with photometric data only
  - A prefix “z” indicates that a redshift is measured from its candidate host galaxy and the classification uses that redshift as a prior. 
  - **SN Ia?:** classification is based on a spectrum that suggests a SN Ia but is inconclusive.
- Some SN candidates have associated notes indicating candidates that may have peculiar features or candidates where the typing spectrum was obtained by other groups. Quoting the data release: "We did not search for these peculiar features in a systematic way, but we have noted the likely peculiar features that were found." Note flags are as follows:
  - 1  SN typing based on spectra obtained by groups outside SDSS. The spectra used for typing are not included in the data release. 
  - 2  Peculiar type Ia SN possibly similar to sn91bg 
  - 3  Peculiar type Ia SN possibly similar to sn00cx 
  - 4  Peculiar type Ia SN possibly similar to sn02ci 
  - 5  Peculiar type Ia SN possibly similar to sn02cx 
- Light curve fits are perfromed using PSNID, SLAT2, and MLCS2k2
  - SALT2 fits are performed both without and with providing a spectroscopic redshift (where available)
  - Only fit data with a photometric flag less than 1024 (Holtzman et al. 2008)

- Quoting section7.2 of the data release:  "The results of the SALT2 fits depend on the version of the code used, the spectral templates, and the color law. Our fits use the SALT2 model as implemented in SNANA version 10.31b and the spectral templates and color law reported in Guy et al. (2010, G10).... For the SDSS data, the largest differences in the fitted param- eters arises from the difference in the color law between G07 and G10. The SDSS-II- SNLS joint light curve anal- ysis paper on cosmology (Betoule et al. 2014) releases a new version of the SALT2 model that is based on adding the full SDSS-II spectroscopically confirmed SN sample to the SALT2 training set." The relationship between the fit color parameters is shown to be linear. The G07 color law results in a value of the *c* parameter that is 20% higher than G10 on average
