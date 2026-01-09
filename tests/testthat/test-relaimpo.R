test_that("equivalent to relaimpo::calc.relimp", {
  
  # skip("long-running test")
  
  packages <- c("zoo", "relaimpo")
  
  status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  
  if (!all(status)) {
    
    skip(paste0(paste0("'", paste0(names(status[!status]), collapse = ", '", sep = "'")),
                "' package(s) required for this test"))
    
  }
  
  # test data
  test_data_x <- c(lapply(test_ls, function(x){x[ , 1:3]}),
                   list("random vector with 0's" = test_ls[[1]][ , 1]))
  test_data_y <- c(lapply(test_ls, function(x){x[ , 4, drop = FALSE]}), # univariate 'y' for relaimpo::calc.relimp
                   list("random vector with 0's" = test_ls[[1]][ , 4]))
  
  for (ax in 1:length(test_data_x)) {
    for (b in 1:length(test_width)) {
      
      width <- test_width[b]     
      # test_weights <- list(rep(1, width))
      test_weights <- list(lambda ^ (width:1))

      for (ay in 1:length(test_data_y)) {
        
        result1 <- roll_shap(test_data_x[[ax]], test_data_y[[ay]],
                             width, test_weights[[1]],
                             test_intercept[1], test_min_obs[1])
        
        # relaimpo: "model must contain intercept"
        result2 <- rollapplyr_relimp(test_data_x[[ax]], test_data_y[[ay]],
                                     width = width, test_weights[[1]],
                                     test_intercept[1])
        
        # relaimpo: uses n > p + 1 instead of n >= p + 1
        # relaimpo: requires at least one degree of freedom for error term
        # parameters to estimate: p + 1 (p features + 1 intercept)
        # observations used: n
        # degrees of freedom for error: n - (p + 1)
        expect_equal(result1[n_vars:n_obs, , drop = FALSE],
                     result2[n_vars:n_obs, , drop = FALSE],
                     check.attributes = FALSE)
        
      }
      
    }
  }
  
})