rollapplyr_shap <- function(f, x, y, width, weights, intercept) {
  
  x <- as.matrix(zoo::coredata(x))
  y <- as.matrix(zoo::coredata(y))
  
  n_rows <- nrow(x)
  n_cols <- ncol(x)
  result <- matrix(as.numeric(NA), n_rows, n_cols)
  
  # parameters to estimate: p features (+ 1 intercept)
  if (intercept) {
    df_fit <- n_cols + 1
  } else {
    df_fit <- n_cols
  }
  
  for (i in 1:n_rows) {
    
    x_subset <- x[max(1, i - width + 1):i, , drop = FALSE]
    y_subset <- y[max(1, i - width + 1):i, , drop = FALSE]
    weights_subset <- tail(weights, min(i, width))
    
    # roll::roll_lm: coefficients are estimable if n >= p (+ 1) and full rank
    if (intercept) {
      rank_fit <- qr(cbind(1, x_subset))$rank
    } else {
      rank_fit <- qr(x_subset)$rank
    }
    
    # roll::roll_lm: r-squared is undefined if 'y' is zero
    var_y <- crossprod(y_subset)
    
    if ((rank_fit == df_fit) && !is.na(var_y) && (var_y > .Machine$double.eps)) {
      
      status_fit <- tryCatch(f(x_subset, y_subset, weights_subset, intercept),
                             
                             error = function(x) {
                               rep(NA, n_cols)
                             })
      
      result[i, ] <- status_fit
      
    }
    
  }
  
  return(result)
  
}
