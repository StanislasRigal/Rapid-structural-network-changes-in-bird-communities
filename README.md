# README for code and data of Rapid structural network changes in bird communities

## Introduction

This repository provides R script and data used to assess network changes in French common bird communities between 2001 and 2017, reported in Rigal, S., Devictor, V., & Dakos, V. (2025). Rapid structural network changes in bird communities. bioRxiv, 2025-02 (https://doi.org/10.1101/2025.02.14.638220).

## Prerequisites

The R scripts have been implemented on R version 4.2.2. All the required R packages with their version are listed in "packages.R".

## Contents 

### Folder structure

`output`: provides intermediate datasets constructed along the analysis  
`raw_data`: contains the initial datasets used to run the whole analysis that came from the French Breeding Bird Survey (STOC-EPS https://www.vigienature.fr/fr/suivi-temporel-des-oiseaux-communs-stoc)  
`function_XXX.R`: scripts of functions for the analysis  
`packages.R`: R packages required  
`main_code.Rmd`: R script with explanations

### File formats 

Raw data are provided both in `csv` (and zipped csv) and `rds`. Output data are in `rds` and can be used in R.

## Usage

Please refer to the`main_code.Rmd` file for usage details.

## License

This work is licensed under Creative Commons Attribution 4.0 International.

## Contact 

Please contact the corresponding author at stanislas.rigal[at]inrae.fr.