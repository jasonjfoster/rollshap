dimnames_lm_x <- function(dimnames, n_cols_x, intercept) {
  
  if (intercept && (length(dimnames) > 1)) {
    
    result <- list(dimnames[[1]], c("(Intercept)", dimnames[[2]]))
    
  } else if (!intercept && (length(dimnames) > 1)) {
    
    result <- list(dimnames[[1]], dimnames[[2]])
    
  } else if (intercept) {
    
    result <- list(NULL, c("(Intercept)", paste0("x", rep(1:(n_cols_x - 1)))))
    
  } else {
    
    result <- list(NULL, paste0("x", rep(1:n_cols_x)))
    
  }
  
  return(result)
  
}

rollapplyr_relimp <- function(x, y, width, intercept) {
  
  if (is.matrix(x) || is.matrix(y) || intercept) {
    
    if (!is.matrix(x)) {
      
      temp_attr <- attributes(x)
      x <- as.matrix(zoo::coredata(x))
      attr(x, "dimnames") <- NULL
      attr(x, "index") <- temp_attr[["index"]]
      attr(x, "class") <- temp_attr[["class"]]
      
    }
    
    if (!is.matrix(y)) {
      
      temp_attr <- attributes(y)
      y <- as.matrix(zoo::coredata(y))
      attr(y, "dimnames") <- NULL
      attr(y, "index") <- temp_attr[["index"]]
      attr(y, "class") <- temp_attr[["class"]]
      
    }
    
    n_rows_xy <- nrow(x)
    n_cols_x <- ncol(x)
    
    # if (intercept) {
    #   n_cols_x <- n_cols_x + 1
    # }
    
    result <- matrix(as.numeric(NA), n_rows_xy, n_cols_x)
    
    if (zoo::is.zoo(x)) {
      
      x_attr <- attributes(x)
      x_attr[["dim"]] <- NULL
      x_attr[["dimnames"]] <- NULL
      
    } else if (zoo::is.zoo(y)) {
      
      x_attr <- attributes(y)
      x_attr[["dim"]] <- NULL
      x_attr[["dimnames"]] <- NULL
      
    }
    
    for (i in 1:n_rows_xy) {
      
      x_subset <- x[max(1, i - width + 1):i, , drop = FALSE]
      y_subset <- y[max(1, i - width + 1):i, , drop = FALSE]
      data <- as.data.frame(cbind(y_subset, x_subset))
      
      if (intercept) {
        
        # "Unparseable 'response'; use is deprecated.  Use as.name(.) or `..`!"
        fit <- lm(reformulate(termlabels = ".", response = as.name(names(data)[1])), data = data)
        
      } else {
        fit <- lm(reformulate(termlabels = ".-1", response = as.name(names(data)[1])), data = data)
      }

      summary_fit <- summary(fit)
      summary_fit_coef <- coef(summary_fit)[ , "Estimate"]

      if ((n_cols_x == 1) && (nrow(data) - (n_cols_x + 1) > 0)) { # summary_fit$df[2] > 0)

        # relaimpo: relative importance equals r-squared for univariate case
        fit_rsq <- summary_fit$r.squared
        result[i] <- fit_rsq

      } else if (nrow(data) - (n_cols_x + 1) > 0) {

          fit_rsq <- tryCatch(relaimpo::calc.relimp(fit, type = "lmg", rela = FALSE)@lmg,
                              error = function(x) NA)
      
          result[i, ] <- as.numeric(fit_rsq)

      }

    }
    
    if (exists("x_attr")) {
      attributes(result) <- x_attr
    }
    
    attr(result, "dim") <- c(n_rows_xy, n_cols_x)
    
    x_dimnames <- dimnames(x)
    y_dimnames <- dimnames(y)
    
    attr(result, "dimnames") <- dimnames_lm_x(x_dimnames, n_cols_x, FALSE) # exclude intercept
  
  }
  
  # } else {
    
  #   n_rows_xy <- length(x)
  #   n_cols_x <- 1
    
  #   result <- rep(as.numeric(NA), n_rows_xy)
    
  #   if (zoo::is.zoo(x)) {

  #     x_attr <- attributes(x)
  #     x_attr[["dim"]] <- NULL
  #     x_attr[["dimnames"]] <- NULL

  #   } else if (zoo::is.zoo(y)) {

  #     x_attr <- attributes(y)
  #     x_attr[["dim"]] <- NULL
  #     x_attr[["dimnames"]] <- NULL

  #   }
    
  #   for (i in 1:n_rows_xy) {
      
  #     x_subset <- x[max(1, i - width + 1):i]
  #     y_subset <- y[max(1, i - width + 1):i]
  #     data <- as.data.frame(cbind(y_subset, x_subset))
      
  #     if (intercept) {
  #       fit <- lm(reformulate(termlabels = ".", response = names(data)[1]), data = data)
  #     } else {
  #       fit <- lm(reformulate(termlabels = ".-1", response = names(data)[1]), data = data)
  #     }

  #     summary_fit <- summary(fit)
  #     summary_fit_coef <- coef(summary_fit)[ , "Estimate"]
      
  #     if ((n_cols_x == 1) && (nrow(data) - (n_cols_x + 1) > 0)) { # summary_fit$df[2] > 0)

  #       # relaimpo: relative importance equals r-squared for univariate case
  #       fit_rsq <- summary_fit$r.squared
  #       result[i] <- fit_rsq

  #     } else if (nrow(data) - (n_cols_x + 1) > 0) {

  #         fit_rsq <- tryCatch(relaimpo::calc.relimp(fit, type = "lmg", rela = FALSE)@lmg,
  #                             error = function(x) NA)
      
  #         result[i, ] <- as.numeric(fit_rsq)

  #     }

  #   }
    
  #   if (exists("x_attr")) {
  #     attributes(result) <- x_attr
  #   }
    
  # }
  
  return(result)
  
}
