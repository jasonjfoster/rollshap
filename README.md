# rollshap

[![](https://github.com/jasonjfoster/rollshap/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/jasonjfoster/rollshap/actions/workflows/check-standard.yaml)
[![](https://www.r-pkg.org/badges/version/rollshap)](https://cran.r-project.org/package=rollshap)
[![](https://codecov.io/gh/jasonjfoster/rollshap/graph/badge.svg)](https://app.codecov.io/github/jasonjfoster/rollshap)
[![](https://cranlogs.r-pkg.org/badges/rollshap?color=brightgreen)](https://www.r-pkg.org/pkg/rollshap)

## Overview

'rollshap' provides analytical computation of rolling and expanding Shapley values for time-series data.

The 'rollshap' package decomposes the coefficient of determination (R-squared) of a linear regression into nonnegative contributions from each explanatory variable using the Shapley value from cooperative game theory (Shapley, 1953, <doi:10.1515/9781400881970-018>). For each window, the exact Shapley value is computed by fitting all subsets of the explanatory variables and averaging the marginal contribution to R-squared across all orderings, which returns an order-invariant attribution that sums to the full-model R-squared. Use cases include:

* **Variable importance**: rolling decomposition of explanatory power across factors or features
* **Factor attribution**: quantifying the contribution of each factor to a model's fit through time
* **Feature selection**: identifying variables whose marginal contribution is persistent or transient

The package supports rolling and expanding windows, weights, and handling of missing values via the min_obs, complete_obs, and na_restore arguments. The implementation uses the online and offline algorithms from the 'roll' package to compute rolling and expanding cross-products efficiently, with parallelism across columns and windows provided by 'RcppParallel'.

## Installation

Install the released version from CRAN:

```r
install.packages("rollshap")
```

Or the development version from GitHub:

```r
# install.packages("pak")
pak::pak("jasonjfoster/rollshap")
```

## Usage

Load the package and supply a dataset:

```r
library(rollshap) # roll (>= 1.1.7)

n <- 15
m <- 3
x <- matrix(rnorm(n * m), nrow = n, ncol = m)
y <- rnorm(n)
weights <- 0.9 ^ (n:1)
```

Then, to compute rolling and expanding Shapley values, use the `roll_shap()` function:

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

Handling of missing values is also supported (see the `min_obs`, `complete_obs`, and `na_restore` arguments).

## References

Shapley, L.S. (1953). "A Value for n-Person Games." In *Contributions to the Theory of Games, Volume II*, edited by H.W. Kuhn and A.W. Tucker, 307-317. Princeton University Press. <doi:10.1515/9781400881970-018>
