test_that("equal to online algorithm", {
  
  # test data
  test_roll_x <- c(lapply(test_ls, function(x){x[ , 1:3]}),
                   list("random vector with 0's and NA's" = test_ls[[2]][ , 1]))
  test_roll_y <- c(lapply(test_ls, function(x){x[ , 4:5]}),
                   list("random vector with 0's and NA's" = test_ls[[2]][ , 4]))
  
  if (requireNamespace("zoo", quietly = TRUE)) {
    names(test_roll_x[[3]]) <- zoo::index(test_roll_x[[1]])
    names(test_roll_y[[3]]) <- zoo::index(test_roll_y[[1]])
  }
  
  for (ax in 1:length(test_roll_x)) {
    for (b in 1:length(test_width)) {
      
      width <- test_width[b]
      # test_weights <- list(rep(1, width))
      test_weights <- list(lambda ^ ((2 * width):1))
      
      for (c in 1:length(test_min_obs)) {
        for (d in 1:length(test_complete_obs)) {
          for (e in 1:length(test_na_restore)) {
            for (f in 1:length(test_weights)) {
              for (ay in 1:length(test_roll_y)) {
                for (g in 1:length(test_intercept)) {
                  
                  # "'complete_obs = FALSE' is not supported"
                  expect_equal(suppressWarnings(roll_shap(test_roll_x[[ax]], test_roll_y[[ay]],
                                                          width, test_weights[[f]],
                                                          test_intercept[g], test_min_obs[c],
                                                          test_complete_obs[d], test_na_restore[e],
                                                          test_online[1])),
                               suppressWarnings(roll_shap(test_roll_x[[ax]], test_roll_y[[ay]],
                                                          width, test_weights[[f]],
                                                          test_intercept[g], test_min_obs[c],
                                                          test_complete_obs[d], test_na_restore[e],
                                                          test_online[2])))
                  
                }
              }
            }
          }
        }
      }
      
    }
  }
  
})
