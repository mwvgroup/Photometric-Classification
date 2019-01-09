# SDSS-Classification

This project applies the photometric classification technique from González-Gaitan 
2014 on SDSS data to identify peculiar Type Ia Supernovae (SNe Ia).

## Todo:
- SNCosmo will raise occasional warning about poor fits, bad S/N, and dropped
    bands. These warnings are currently being logged but otherwise ignored and
    need to be handled correctly.
- If a redshift for a supernova is not available, SNCosmo is instructed to fit
    for the redshift. A prior needs to be provided for the redshift.
- The `num_points_<band>` column in the fit summary tables contain the number
    of data points ber **observer frame** band. We are interested in the number
    of points in the **rest frame** 
- We are currently fitting data with SNCosmo's builtin 91bg model, but would
    like to (possibly) use our own custom model. This custom model is currently
    used to simulate SNANA light curves. 

## File List / Overview

#### *./* 

General utilities

 - **download_data.py:** This script downloads data from the SDSS supernova survey to a directory called *data/* .
 - **parse_sn_data.py:** This module parses SDSS data tables. Missing data is downloaded automatically.



#### snana/ 

Simulating 91bg light curves using the SNANA Fortran package

- **sim_91bg.input:** An input file to SNANA that generates simulation results
- **snana_results/:** Directory of simulation results generated by *sim_91bg.input*
- **SNANA_plots.ipynb:** Plots exploring 91bg simulation results.



#### SNCosmo/

Fitting SDSS data with the SNcosmo python package

- **fit_sncosmo.py:** Uses SNCosmo to fit SDSSS light curves with a normal sn Ia model and the builtin 91bg model in the `ug`, `riz`, and `ugriz` rest frame bands.
- **sncosmo_results.ipynb:** Compares fit results between the normal and 91bg models in different bands.
- **sncosmo_vs_sdss.ipynb:** Compares fit results with published SDSS results.



#### snoopy/

Fitting SDSS data with the Snoopy (`snpy`) python package. Development here is unfinished and abandoned after the SNCsomo package was found to be more user friendly.

- **fit_snoopy.py:** Uses Snoopy to fit SDSSS light curves with a normal sn Ia model.