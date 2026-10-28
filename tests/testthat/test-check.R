test_that("input checks", {
  
  # test data
  test_check_x <- c(lapply(test_ls[1], function(x){x[ , 1:3]}),
                    list("random vector with 0's" = test_ls[[1]][ , 1]))
  test_check_y <- test_ls[[1]][ , 4]
  width <- test_width[2]
  
  for (ax in 1:length(test_check_x)) {
    
    # check_rows_equal
    expect_error(roll_shap(test_check_x[[ax]], test_check_y[-1], width),
                 "number of rows in 'x' must equal the number of rows in 'y'")
    
    # check_pos_int: 'width'
    expect_error(roll_shap(test_check_x[[ax]], test_check_y, width = 0),
                 "value of 'width' must be greater than or equal to one")
    
    # check_pos_int: 'min_obs'
    expect_error(roll_shap(test_check_x[[ax]], test_check_y, width, min_obs = 0),
                 "value of 'min_obs' must be greater than or equal to one")
    
    # check_weights
    expect_error(roll_shap(test_check_x[[ax]], test_check_y, width, rep(1, width - 1)),
                 "length of 'weights' must be greater than or equal to the number of rows in 'x' (and 'y') or 'width'",
                 fixed = TRUE)
    
  }
  
  # check_bounds_int: 'n_cols' (matrix 'x' only)
  expect_error(roll_shap(matrix(rnorm(n_obs * 21), nrow = n_obs, ncol = 21), test_check_y, width),
               "value of 'n_cols' must be between 1 and 20")
  
})

test_that("input warnings", {
  
  # test data
  test_check_x <- c(lapply(test_ls[1], function(x){x[ , 1:3]}),
                    list("random vector with 0's" = test_ls[[1]][ , 1]))
  test_check_y <- test_ls[[1]][ , 4]
  width <- test_width[2]
  
  for (ax in 1:length(test_check_x)) {
    
    # check_lambda
    expect_warning(roll_shap(test_check_x[[ax]], test_check_y, width, 1:width),
                   "'online' is only supported for equal or exponential decay 'weights'")
    expect_silent(roll_shap(test_check_x[[ax]], test_check_y, width, 1:width, online = FALSE))
    
    # 'complete_obs = FALSE'
    expect_warning(roll_shap(test_check_x[[ax]], test_check_y, width, complete_obs = FALSE),
                   "'complete_obs = FALSE' is not supported")
    
  }
  
})
