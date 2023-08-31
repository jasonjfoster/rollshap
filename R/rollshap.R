##' Rolling Shapley Values
##'
##' A function for computing the rolling and expanding Shapley values of time-series data.
##' 
##' @param x vector or matrix. Rows are observations and columns are the independent variables.
##' @param y vector or matrix. Rows are observations and columns are the dependent variables.
##' @param width integer. Window size.
##' @param weights vector. Weights for each observation within a window.
##' @param intercept logical. Either \code{TRUE} to include or \code{FALSE} to remove the intercept.
##' @param min_obs integer. Minimum number of observations required to have a value within a window,
##' otherwise result is \code{NA}.
##' @param complete_obs	logical. If \code{TRUE} then rows containing any missing values are removed,
##' if \code{FALSE} then pairwise is used.
##' @param na_restore logical. Should missing values be restored?
##' @param online logical. Process observations using an online algorithm.
##' @return An object of the same class and dimension as \code{x} with the rolling and expanding
##' Shapley values.
##' @examples
##' n <- 15
##' m <- 3
##' x <- matrix(rnorm(n * m), nrow = n, ncol = m)
##' y <- rnorm(n)
##' weights <- 0.9 ^ (n:1)
##' 
##' # rolling Shapley values with complete windows
##' roll_shap(x, y, width = 5)
##' 
##' # rolling Shapley values with partial windows
##' roll_shap(x, y, width = 5, min_obs = 1)
##' 
##' # expanding Shapley values with partial windows
##' roll_shap(x, y, width = n, min_obs = 1)
##' 
##' # expanding Shapley values with partial windows and weights
##' roll_shap(x, y, width = n, min_obs = 1, weights = weights)
##' @export
roll_shap <- function(x, y, width, weights = rep(1, width), intercept = TRUE,
                      min_obs = width, complete_obs = TRUE, na_restore = FALSE,
                      online = TRUE) {
  return(.Call(`_rollshap_roll_shap`,
               x, y,
               as.integer(width),
               as.numeric(weights),
               as.logical(intercept),
               as.integer(min_obs),
               as.logical(complete_obs),
               as.logical(na_restore),
               as.logical(online)
  ))
}