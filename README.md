# trait-ipm: trait-based models of plant population dynamics 
Scripts to load plant demographic and trait data from the COMPADRE and BIEN databases, impute missing trait data using probabilistic matrix factorization, and fit a functional data model that relates size-structured demographic vital rates to traits. These scripts are supporting material to Yen et al. (submitted) Predicting size-structured population dynamics from traits.

Copyright &copy; 2018, Jian Yen

*****

## License details
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

*****

## Overview
This repository contains scripts to run the analysis presented in:
Yen JDL, Golding N, Vesk P (submitted) Predicting size-structured population dynamics from traits. Submitted to J Ecol.

Maintainer: Jian Yen (jdl.yen@gmail.com)

Updated: 19 March 2018

## Usage
Several scripts and pre-compiled data files are provided:
-main.R runs the central analysis that relates size-structured vital rates to traits. This script requires pre-compiled data sets (.RData binaries) and requires a working installation of Python 3 with the GPflow library.  
-load-compadre-data.R loads the compadre data from a pre-downloaded .RData binary (COMPARE_v.x.y.z.RData; available at www.compadre-db.org/.  
-load-trait-data-bien.R downloads the BIEN trait data using the BIEN R package.  
-impute-trait-data.R imputes missing trait data using a hierarchical Bayesian variation of probabilistic matrix factorization (see Schrodt et al. 2015, GEB 24:1510-1521 for details of the method). This model was implemented in Stan, using the rstan R package.  
-several additional scripts contain helper functions for the above scripts and should not need to be run directly.  

## Data availability
All analyses using the COMPADRE and BIEN databases and all users should refer to their respective websites (www.compadre-db.org/ and biendata.org). For reproducibility purposes, users can re-run our analysis from scratch using the COMPADRE database v.4.0.1 and the included script to download BIEN trait data. Each step can take some time (>30 min), so we also provide pre-compiled .RData files with the trait and demographic data used in our analysis. These files are available at figshare.com/not_added_yet and are provided for reproducibility purposes only. Any users wishing to use COMPADRE or BIEN data should refer to the websites listed above and follow their respective data use policies.

