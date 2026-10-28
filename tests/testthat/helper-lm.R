lm_rsq <- function(x, y, weights, intercept) {
  
  if (ncol(x) == 0) {
    return(0)
  }
  
  data <- as.data.frame(cbind(y, x))
  
  if (intercept) {
    
    fit <- lm(reformulate(termlabels = ".", response = as.name(names(data)[1])),
              weights = weights, data = data)
    tss <- sum(weights * (y - weighted.mean(y, weights)) ^ 2)
    
  } else {
    
    # uncentered total sum of squares (consistent with summary.lm)
    fit <- lm(reformulate(termlabels = ".-1", response = as.name(names(data)[1])),
              weights = weights, data = data)
    tss <- sum(weights * y ^ 2)
    
  }
  
  rss <- sum(weights * residuals(fit) ^ 2)
  result <- 1 - rss / tss
  
  return(result)
  
}

lm_shap <- function(x, y, weights, intercept) {
  
  x <- as.matrix(zoo::coredata(x))
  y <- as.numeric(zoo::coredata(y))
  
  n_cols_x <- ncol(x)
  n_combn <- 2 ^ n_cols_x
  
  # binary combinations of the columns of 'x'
  ix <- outer(0:(n_combn - 1), 2 ^ (0:(n_cols_x - 1)), bitwAnd) != 0
  n_ix <- rowSums(ix)
  
  # Shapley weights by number of columns: s! (p - s - 1)! / p!
  shap_weights <- factorial(0:(n_cols_x - 1)) * factorial(n_cols_x - 1 - 0:(n_cols_x - 1)) /
    factorial(n_cols_x)
  
  rsq <- rep(as.numeric(NA), n_combn)
  
  for (k in 1:n_combn) {
    rsq[k] <- lm_rsq(x[ , ix[k, ], drop = FALSE], y, weights, intercept)
  }
  
  result <- rep(as.numeric(NA), n_cols_x)
  
  # exact Shapley value for r-squared
  for (j in 1:n_cols_x) {
    
    ix_neg <- which(!ix[ , j])
    ix_pos <- ix_neg + 2 ^ (j - 1)
    
    result[j] <- sum(shap_weights[n_ix[ix_neg] + 1] * (rsq[ix_pos] - rsq[ix_neg]))
    
  }
  
  return(result)
  
}
