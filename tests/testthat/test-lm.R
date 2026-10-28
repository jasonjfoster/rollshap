test_that("equivalent to stats::lm", {
  
  packages <- c("zoo")
  
  status <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  
  if (!all(status)) {
    
    skip(paste0(paste0("'", paste0(names(status[!status]), collapse = ", '", sep = "'")),
                "' package(s) required for this test"))
    
  }
  
  # test data (non-NA dataset only: lm drops NA rows differently than the C++ windows)
  test_data_x <- c(lapply(test_ls[1], function(x){x[ , 1:3]}),
                   list("random vector with 0's" = test_ls[[1]][ , 1]))
  test_data_y <- c(lapply(test_ls[1], function(x){x[ , 4, drop = FALSE]}), # univariate 'y' for lm_shap
                   list("random vector with 0's" = test_ls[[1]][ , 4]))
  
  for (ax in 1:length(test_data_x)) {
    for (b in 1:length(test_width)) {
      
      width <- test_width[b]     
      # test_weights <- list(rep(1, width))
      test_weights <- list(lambda ^ (width:1))

      for (ay in 1:length(test_data_y)) {
        for (g in 1:length(test_intercept)) {
          
          # as.matrix: vector 'x' returns a vector
          expect_equal(as.matrix(roll_shap(test_data_x[[ax]], test_data_y[[ay]],
                                           width, test_weights[[1]],
                                           test_intercept[g], test_min_obs[1]))[n_vars:n_obs, , drop = FALSE],
                       rollapplyr_shap(lm_shap, test_data_x[[ax]], test_data_y[[ay]],
                                       width, test_weights[[1]],
                                       test_intercept[g])[n_vars:n_obs, , drop = FALSE],
                       check.attributes = FALSE)
          
        }
      }
      
    }
  }
  
})
