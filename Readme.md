# RXMISS

R code and data for experiments on **random-X prediction error under model misspecification**.

------------------------------------------------------------------------

## Repository structure

``` text
RXMISS/
  CODE/
    Realdata_B.R
    exchange_search_utils.R
    plot_graphic_representation.R
    smoothing_spline_utility.R
  DATA/
    diamonds.csv
    soil.csv
  OUTPUT/
    (generated figures and results)
  real_data_example.Rproj
```

------------------------------------------------------------------------

## Overview

This repository contains code to study subsampling strategies for linear regression under model misspecification, with a focus on comparing a proposed method against several alternatives.

The main script, `CODE/Realdata_B.R`, loads the data, fits full-sample reference estimators, generates subsamples using different methods, and compares their estimation and prediction performance.

The methods considered include: - uniform subsampling\
- leverage-based subsampling\
- shrinked leverage sampling\
- IBOSS\
- LowCon\
- the proposed MISS approach

------------------------------------------------------------------------

## Data

The repository includes two datasets:

-   `DATA/soil.csv`\
-   `DATA/diamonds.csv`

The script is currently configured to run the **soil dataset** by default. The diamonds example is included but commented out and can be activated by editing the script.

------------------------------------------------------------------------

## Main files

### `CODE/Realdata_B.R`

Main script for running the real-data experiments. It: - loads required packages\
- reads and preprocesses the data\
- computes benchmark estimators\
- runs repeated subsampling experiments\
- reports estimation error and MSPE

### `CODE/exchange_search_utils.R`

Utility functions for the exchange algorithm used in the proposed method.

### `CODE/smoothing_spline_utility.R`

Auxiliary functions for spline-based estimation.

### `CODE/plot_graphic_representation.R`

Functions for visualization of results.

------------------------------------------------------------------------

## Requirements

The code uses the following R packages:

```         
crmReg, MASS, gss, foreach, doParallel, randtoolbox, nabor, mvtnorm, future, progressr, doFuture
```

------------------------------------------------------------------------

## How to run

1.  Clone the repository:

``` bash
git clone https://github.com/alvarocia/RXMISS.git
cd RXMISS
```

2.  Open `real_data_example.Rproj` in RStudio, or start an R session in the project directory.

3.  Run the main script:

``` r
source("CODE/Realdata_B.R")
```

------------------------------------------------------------------------

## Notes

-   The script runs repeated simulations and may take time depending on the configuration.\
-   Parallel computation is used via `future` and `doFuture`.\
-   Parameters such as subsample size, number of iterations, and dataset can be modified directly in the main script.\
-   The `OUTPUT` folder contains generated results and can be ignored or regenerated.

------------------------------------------------------------------------

## Acknowledgment

Part of the code is adapted from the methodology and experimental framework introduced in the LowCon paper:

Meng, C., Xie, R., Mandal, A., Zhang, X., Zhong, W., and Ma, P. (2021). *LowCon: A design-based subsampling approach in a misspecified linear model*. Journal of Computational and Graphical Statistics, 30(3), 694–708.

The present implementation builds upon and modifies parts of that framework for the purposes of this project.

------------------------------------------------------------------------

## Author

**alvarocia**
