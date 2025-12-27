set.seed(5640)
n_vars <- 5
n_obs <- 15 # 15 # 1500
n_size <- n_obs * n_vars
lambda <- 0.9 # 0.9 # 1 / 0.9
dates <- rev(seq(Sys.Date(), length.out = n_obs, by = "-1 day"))

test_width <- c(1, 5, 10, 15)
test_intercept <- c(TRUE, FALSE)
test_min_obs <- c(1, 5, 10, 15)
test_complete_obs <- c(TRUE, FALSE)
test_na_restore <- c(TRUE, FALSE)
test_online <- c(TRUE, FALSE)

# test data
test_ls <- list("random matrix with 0's" =
                  matrix(sample(c(0, rnorm(n_size - 1)), n_size, replace = TRUE,
                                prob = rep(1 / n_size, n_size)),
                         nrow = n_obs, ncol = n_vars))

if (requireNamespace("zoo", quietly = TRUE)) {
  test_ls[[1]] <- zoo::zoo(test_ls[[1]], dates)
}

test_ls[[1]] <- setNames(test_ls[[1]], paste0("x", rep(1:n_vars)))

# colnames(test_ls[[1]]) <- paste0("x", rep(1:n_vars))