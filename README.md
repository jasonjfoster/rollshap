# rollshap

[![](https://github.com/jasonjfoster/rollshap/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/jasonjfoster/rollshap/actions/workflows/check-standard.yaml)

## Overview

`rollshap` is a package that provides analytical computation of rolling and expanding Shapley values for time-series data.

## Installation

Install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("jasonjfoster/rollshap")
```

## Usage

Load the package and supply a dataset:

``` r
library(rollshap) # roll (>= 1.1.7)

n <- 15
m <- 3
x <- matrix(rnorm(n * m), nrow = n, ncol = m)
y <- rnorm(n)
weights <- 0.9 ^ (n:1)
```
Then, to compute rolling and expanding Shapley values, use the `roll_shap` function:

```r
# rolling Shapley values with complete windows
roll_shap(x, y, width = 5)

# rolling Shapley values with partial windows
roll_shap(x, y, width = 5, min_obs = 1)

# expanding Shapley values with partial windows
roll_shap(x, y, width = n, min_obs = 1)

# expanding Shapley values with partial windows and weights
roll_shap(x, y, width = n, min_obs = 1, weights = weights)
```

Note that handling of missing values is supported as well (see the `min_obs`, `complete_obs`, and `na_restore` arguments).