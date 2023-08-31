#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <roll.h>
using namespace Rcpp;
using namespace RcppParallel;

void check_width(const int& width) {
  
  if (width < 1) {
    stop("value of 'width' must be greater than zero");
  }
  
}

void check_weights_lm(const int& n_rows_xy, const int& width,
                      const arma::vec& weights) {
  
  if ((int)weights.size() < std::min(width, n_rows_xy)) {
    stop("length of 'weights' must equal either the number of rows in 'x' (and 'y') or 'width'");
  }
  
}

bool check_lambda(const arma::vec& weights, const int& n_rows_x,
                  const int& width, const bool& online) {
  
  // check if equal-weights
  bool status_eq = all(weights == weights[0]);
  bool status_exp = true;
  
  // check if exponential-weights
  if (!status_eq) {
    
    int i = 0;
    int n = weights.size();
    long double lambda = 0;
    long double lambda_prev = 0;
    
    // check if constant ratio
    while (status_exp && (i <= (n - 2))) {
      
      // ratio of weights
      lambda_prev = lambda;
      lambda = weights[n - i - 2] / weights[n - i - 1];
      
      // tolerance for consistency with R's all.equal
      if (((i > 0) && (std::abs(lambda - lambda_prev) > sqrt(arma::datum::eps))) ||
          ((weights[n - i - 2] > weights[n - i - 1]) && (width < n_rows_x)) ||
          (std::isnan(lambda) || (std::isinf(lambda)))) {
        
        status_exp = false;
        
      }
      
      i += 1;
      
    }
    
  }
  
  if (!status_exp && online) {
    warning("'online' is only supported for equal or exponential decay 'weights'");
  }
  
  return status_exp;
  
}

void check_min_obs(const int& min_obs) {
  
  if (min_obs < 1) {
    stop("value of 'min_obs' must be greater than zero");
  }
  
}

void check_lm(const int& n_rows_x, const int& n_rows_y) {
  
  if (n_rows_x != n_rows_y) {
    stop("number of rows in 'x' must equal the number of rows in 'y'");
  }
  
}

List dimnames_lm_x(const List& input, const int& n_cols_x,
                   const bool& intercept) {
  
  if (intercept && (input.size() > 1)) {
    
    CharacterVector dimnames_cols = input[1];
    CharacterVector result(n_cols_x);
    result(0) = "(Intercept)";
    
    std::copy(dimnames_cols.begin(), dimnames_cols.end(), result.begin() + 1);
    
    return List::create(input[0], result);
    
  } else if (!intercept && (input.size() > 1)) {
    
    return List::create(input[0], input[1]);
    
  } else if (intercept) {
    
    CharacterVector result(n_cols_x);
    result(0) = "(Intercept)";
    
    for (int i = 1; i < n_cols_x; i++) {
      
      result[i] = "x";
      result[i] += i;
      
    }
    
    return List::create(R_NilValue, result);
    
  } else {
    
    CharacterVector result(n_cols_x);
    
    for (int i = 0; i < n_cols_x; i++) {
      
      result[i] = "x";
      result[i] += i + 1;
      
    }
    
    return List::create(R_NilValue, result);
    
  }
  
}

CharacterVector dimnames_lm_y(const List& input, const int& n_cols_y) {
  
  if (input.size() > 1) {
    
    return input[1];
    
  } else {
    
    CharacterVector result(n_cols_y);
    
    for (int i = 0; i < n_cols_y; i++) {
      
      result[i] = "y";
      result[i] += i + 1;
      
    }
    
    return result;
    
  }
  
}

arma::uvec any_na_x(const NumericMatrix& x) {
  
  int n_rows_x = x.nrow();
  int n_cols_x = x.ncol();
  arma::uvec result(n_rows_x);
  
  for (int i = 0; i < n_rows_x; i++) {
    
    int any_na = 0;
    int j = 0;
    
    while ((any_na == 0) && (j < n_cols_x)) {
      if (std::isnan(x(i, j))) {
        any_na = 1;
      }
      j += 1;
    }
    
    result[i] = any_na;
    
  }
  
  return result;
  
} 

List roll_lm_z(const SEXP& x, const NumericVector& y,
               const int& width, const arma::vec& weights,
               const bool& intercept, const int& min_obs,
               const bool& complete_obs, const bool& na_restore,
               const bool& online) {
  
  if (Rf_isMatrix(x)) {
    
    NumericMatrix xx(x);
    int n = weights.size();
    int n_rows_xy = xx.nrow();
    int n_cols_x = xx.ncol() + 1;
    arma::uvec arma_any_na(n_rows_xy);
    arma::vec arma_n_obs(n_rows_xy);
    arma::vec arma_sum_w(n_rows_xy);
    arma::mat arma_mean(n_rows_xy, n_cols_x);
    arma::cube arma_cov(n_cols_x, n_cols_x, n_rows_xy);
    arma::vec arma_rsq(n_rows_xy);
    List result(3);
    
    // check 'x' and 'y' arguments for errors
    check_lm(n_rows_xy, y.size());
    
    // check 'width' argument for errors
    check_width(width);
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    check_weights_lm(n_rows_xy, width, weights);
    bool status = check_lambda(weights, n_rows_xy, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    check_min_obs(min_obs);
    
    // cbind x and y variables
    NumericMatrix data(n_rows_xy, n_cols_x);
    std::copy(xx.begin(), xx.end(), data.begin());
    std::copy(y.begin(), y.end(), data.begin() + n_rows_xy * (n_cols_x - 1));
    
    // default 'complete_obs' argument is 'true',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_x(data);
    } else {
      
      warning("'complete_obs = FALSE' is not supported");
      arma_any_na = any_na_x(data);
      
    }
    
    // compute rolling covariances
    if (status && online) {
      
      roll::RollCovOnlineMatLm roll_cov_online(data, n, n_rows_xy, n_cols_x, width,
                                               weights, intercept, min_obs,
                                               arma_any_na, na_restore,
                                               arma_n_obs, arma_sum_w, arma_mean,
                                               arma_cov);
      parallelFor(0, n_cols_x, roll_cov_online);
      
    } else {
      
      roll::RollCovOfflineMatLm roll_cov_offline(data, n, n_rows_xy, n_cols_x, width,
                                                 weights, intercept, min_obs,
                                                 arma_any_na, na_restore,
                                                 arma_n_obs, arma_sum_w, arma_mean,
                                                 arma_cov);
      parallelFor(0, n_rows_xy * n_cols_x * (n_cols_x + 1) / 2, roll_cov_offline);
      
    }
    
    // compute rolling linear models
    if (intercept) {
      
      arma::mat arma_coef(n_rows_xy, n_cols_x);
      arma::mat arma_se(n_rows_xy, n_cols_x);
      roll::RollLmMatInterceptTRUE roll_lm_slices(arma_cov, n, n_rows_xy, n_cols_x, width,
                                                  arma_n_obs, arma_sum_w, arma_mean,
                                                  arma_coef, arma_rsq, arma_se);
      parallelFor(0, n_rows_xy, roll_lm_slices);
      
      result = List::create(Named("coefficients") = arma_coef,
                            Named("r.squared") = arma_rsq,
                            Named("std.error") = arma_se);
      
    } else if (!intercept) {
      
      arma::mat arma_coef(n_rows_xy, n_cols_x - 1);
      arma::mat arma_se(n_rows_xy, n_cols_x - 1);
      roll::RollLmMatInterceptFALSE roll_lm_slices(arma_cov, n, n_rows_xy, n_cols_x, width,
                                                   arma_n_obs, arma_sum_w,
                                                   arma_coef, arma_rsq, arma_se);
      parallelFor(0, n_rows_xy, roll_lm_slices);
      
      // create and return a list
      result = List::create(Named("coefficients") = arma_coef,
                            Named("r.squared") = arma_rsq,
                            Named("std.error") = arma_se);
      
    }
    
    return result;
    
  } else {
    
    NumericVector xx(x);
    int n = weights.size();
    int n_rows_xy = xx.size();
    int n_cols_x = 1 + 1;
    arma::uvec arma_any_na(n_rows_xy);
    arma::vec arma_n_obs(n_rows_xy);
    arma::vec arma_sum_w(n_rows_xy);
    arma::mat arma_mean(n_rows_xy, n_cols_x);
    arma::cube arma_cov(n_cols_x, n_cols_x, n_rows_xy);
    arma::vec arma_rsq(n_rows_xy);
    List result(3);
    
    // check 'x' and 'y' arguments for errors
    check_lm(n_rows_xy, y.size());
    
    // check 'width' argument for errors
    check_width(width);
    
    // default 'weights' argument is equal-weighted,
    // otherwise check argument for errors
    check_weights_lm(n_rows_xy, width, weights);
    bool status = check_lambda(weights, n_rows_xy, width, online);
    
    // default 'min_obs' argument is 'width',
    // otherwise check argument for errors
    check_min_obs(min_obs);
    
    // cbind x and y variables
    NumericMatrix data(n_rows_xy, n_cols_x);
    std::copy(xx.begin(), xx.end(), data.begin());
    std::copy(y.begin(), y.end(), data.begin() + n_rows_xy * (n_cols_x - 1));
    
    // default 'complete_obs' argument is 'true',
    // otherwise check argument for errors
    if (complete_obs) {
      arma_any_na = any_na_x(data);
    } else {
      
      warning("'complete_obs = FALSE' is not supported");
      arma_any_na = any_na_x(data);
      
    }
    
    // compute rolling covariances
    if (status && online) {
      
      roll::RollCovOnlineMatLm roll_cov_online(data, n, n_rows_xy, n_cols_x, width,
                                               weights, intercept, min_obs,
                                               arma_any_na, na_restore,
                                               arma_n_obs, arma_sum_w, arma_mean,
                                               arma_cov);
      parallelFor(0, n_cols_x, roll_cov_online);
      
    } else {
      
      roll::RollCovOfflineMatLm roll_cov_offline(data, n, n_rows_xy, n_cols_x, width,
                                                 weights, intercept, min_obs,
                                                 arma_any_na, na_restore,
                                                 arma_n_obs, arma_sum_w, arma_mean,
                                                 arma_cov);
      parallelFor(0, n_rows_xy * n_cols_x * (n_cols_x + 1) / 2, roll_cov_offline);
      
    }
    
    // compute rolling linear models
    if (intercept) {
      
      arma::mat arma_coef(n_rows_xy, n_cols_x);
      arma::mat arma_se(n_rows_xy, n_cols_x);
      roll::RollLmMatInterceptTRUE roll_lm_slices(arma_cov, n, n_rows_xy, n_cols_x, width,
                                                  arma_n_obs, arma_sum_w, arma_mean,
                                                  arma_coef, arma_rsq, arma_se);
      parallelFor(0, n_rows_xy, roll_lm_slices);
      
      result = List::create(Named("coefficients") = arma_coef,
                            Named("r.squared") = arma_rsq,
                            Named("std.error") = arma_se);
      
    } else if (!intercept) {
      
      arma::vec arma_coef(n_rows_xy);
      arma::vec arma_se(n_rows_xy);
      roll::RollLmVecInterceptFALSE roll_lm_slices(arma_cov, n, n_rows_xy, width,
                                                   arma_n_obs, arma_sum_w,
                                                   arma_coef, arma_rsq, arma_se);
      parallelFor(0, n_rows_xy, roll_lm_slices);
      
      // create and return a list
      result = List::create(Named("coefficients") = arma_coef,
                            Named("r.squared") = arma_rsq,
                            Named("std.error") = arma_se);
      
    }
    
    return result;
    
  }
  
}

int factorial(int n) {
  
  int result = 1;
  
  for (int i = 2; i <= n; i++) {
    result *= i;
  }
  
  return result;
  
}

// [[Rcpp::export(.roll_shap)]]
SEXP roll_shap(const SEXP& x, const SEXP& y,
               const int& width, const arma::vec& weights,
               const bool& intercept, const int& min_obs,
               const bool& complete_obs, const bool& na_restore,
               const bool& online) {
  
  if (Rf_isMatrix(x) && Rf_isMatrix(y)) {
    
    NumericMatrix xx(x);
    NumericMatrix yy(y);
    
    int n = 0;
    int n_size = 0;
    int n_rows_xy = xx.nrow();
    int n_cols_x = xx.ncol();
    int n_cols_y = yy.ncol();
    int n_combn = pow((long double)2.0, n_cols_x);
    arma::mat arma_x = arma::mat(xx.begin(), n_rows_xy, n_cols_x);
    arma::ivec arma_n(n_combn);
    arma::umat arma_ix(n_cols_x, n_combn);
    arma::mat arma_rsq(n_rows_xy, n_combn);
    arma::mat arma_rsq_sum(n_rows_xy, n_cols_x);
    List result_rsq(n_cols_y);
    List result_z(3);
    
    // create a list of matrices,
    // otherwise a list of lists
    if (n_cols_y == 1) {
      
      // find all possible combinations of binary values
      for (int k = 0; k < n_combn; k++) {
        
        n = 0;
        n_size = k;
        
        for (int j = 0; j < n_cols_x; j++) {
          
          if (n_size % 2 == 0) {
            
            n += 1;
            
            arma_ix(j, k) = j + 1;
            
          }
          
          n_size /= 2;
          
        }
        
        arma_n[k] = n;
        
        if (n > 0) {
          
          arma::uvec arma_ix_subset = find(arma_ix.col(k));
          arma::mat arma_x_subset = arma_x.cols(arma_ix_subset);
          NumericMatrix x_subset(wrap(arma_x_subset));
          
          result_z = roll_lm_z(x_subset, yy(_, 0), width,
                               weights, intercept, min_obs,
                               complete_obs, na_restore,
                               online);
          
          arma::mat arma_rsq_z = result_z[1];
          arma_rsq.col(k) = arma_rsq_z;
          
        }
        
      }
      
      // calculate the exact Shapley value for r-squared
      for (int j = 0; j < n_cols_x; j++) {
        
        arma::uvec arma_ix_pos = find(arma_ix.row(j));
        arma::uvec arma_ix_neg = find(arma_ix.row(j) == 0);
        arma::ivec arma_ix_n = arma_n(arma_ix_neg);
        arma::mat arma_rsq_diff = arma_rsq.cols(arma_ix_pos) - arma_rsq.cols(arma_ix_neg);
        
        for (int k = 0; k < n_combn / 2; k++) {
          
          int s = arma_ix_n[k];
          long double weight = (factorial(s) * factorial(n_cols_x - s - 1)) /
            (long double)factorial(n_cols_x);
          
          arma_rsq_sum.col(j) += weight * arma_rsq_diff.col(k);
          
        }
        
      }
      
      // create and return a matrix or xts object for Shapley values
      NumericMatrix rsq(wrap(arma_rsq_sum));
      List dimnames_x = xx.attr("dimnames");
      rsq.attr("dimnames") = dimnames_x;
      rsq.attr("index") = xx.attr("index");
      rsq.attr(".indexCLASS") = xx.attr(".indexCLASS");
      rsq.attr(".indexTZ") = xx.attr(".indexTZ");
      rsq.attr("tclass") = xx.attr("tclass");
      rsq.attr("tzone") = xx.attr("tzone");
      rsq.attr("class") = xx.attr("class");
      
      return rsq;
      
    } else {
      
      for (int z = 0; z < n_cols_y; z++) {
        
        // find all possible combinations of binary values
        for (int k = 0; k < n_combn; k++) {
          
          n = 0;
          n_size = k;
          
          for (int j = 0; j < n_cols_x; j++) {
            
            if (n_size % 2 == 0) {
              
              n += 1;
              
              arma_ix(j, k) = j + 1;
              
            }
            
            n_size /= 2;
            
          }
          
          arma_n[k] = n;
          
          if (n > 0) {
            
            arma::uvec arma_ix_subset = find(arma_ix.col(k));
            arma::mat arma_x_subset = arma_x.cols(arma_ix_subset);
            NumericMatrix x_subset(wrap(arma_x_subset));
            
            result_z = roll_lm_z(x_subset, yy(_, z), width,
                                 weights, intercept, min_obs,
                                 complete_obs, na_restore,
                                 online);
            
            arma::mat arma_rsq_z = result_z[1];
            arma_rsq.col(k) = arma_rsq_z;
            
          }
          
        }
        
        // calculate the exact Shapley value for r-squared
        for (int j = 0; j < n_cols_x; j++) {
          
          arma::uvec arma_ix_pos = find(arma_ix.row(j));
          arma::uvec arma_ix_neg = find(arma_ix.row(j) == 0);
          arma::ivec arma_ix_n = arma_n(arma_ix_neg);
          arma::mat arma_rsq_diff = arma_rsq.cols(arma_ix_pos) - arma_rsq.cols(arma_ix_neg);
          
          for (int k = 0; k < n_combn / 2; k++) {
            
            int s = arma_ix_n[k];
            long double weight = (factorial(s) * factorial(n_cols_x - s - 1)) /
              (long double)factorial(n_cols_x);
            
            arma_rsq_sum.col(j) += weight * arma_rsq_diff.col(k);
            
          }
          
        }
        
        // create and return a matrix or xts object for Shapley values
        NumericMatrix rsq(wrap(arma_rsq_sum));
        List dimnames_x = xx.attr("dimnames");
        rsq.attr("dimnames") = dimnames_x;
        rsq.attr("index") = xx.attr("index");
        rsq.attr(".indexCLASS") = xx.attr(".indexCLASS");
        rsq.attr(".indexTZ") = xx.attr(".indexTZ");
        rsq.attr("tclass") = xx.attr("tclass");
        rsq.attr("tzone") = xx.attr("tzone");
        rsq.attr("class") = xx.attr("class");
        
        result_rsq(z) = rsq;
        
      }
      
      // add names to each list
      List dimnames_y = yy.attr("dimnames");
      result_rsq.attr("names") = dimnames_lm_y(dimnames_y, n_cols_y);
      
      return result_rsq;
      
    }
    
  } else if (Rf_isMatrix(x)) {
    
    NumericMatrix xx(x);
    NumericVector yy(y);
    
    int n = 0;
    int n_size = 0;
    int n_rows_xy = xx.nrow();
    int n_cols_x = xx.ncol();
    int n_combn = pow((long double)2.0, n_cols_x);
    arma::mat arma_x = arma::mat(xx.begin(), n_rows_xy, n_cols_x);
    arma::ivec arma_n(n_combn);
    arma::umat arma_ix(n_cols_x, n_combn);
    arma::mat arma_rsq(n_rows_xy, n_combn);
    arma::mat arma_rsq_sum(n_rows_xy, n_cols_x);
    List result_rsq(1);
    List result_z(3);
    
    // find all possible combinations of binary values
    for (int k = 0; k < n_combn; k++) {
      
      n = 0;
      n_size = k;
      
      for (int j = 0; j < n_cols_x; j++) {
        
        if (n_size % 2 == 0) {
          
          n += 1;
          
          arma_ix(j, k) = j + 1;
          
        }
        
        n_size /= 2;
        
      }
      
      arma_n[k] = n;
      
      if (n > 0) {
        
        arma::uvec arma_ix_subset = find(arma_ix.col(k));
        arma::mat arma_x_subset = arma_x.cols(arma_ix_subset);
        NumericMatrix x_subset(wrap(arma_x_subset));
        
        result_z = roll_lm_z(x_subset, yy, width,
                             weights, intercept, min_obs,
                             complete_obs, na_restore,
                             online);
        
        arma::mat arma_rsq_z = result_z[1];
        arma_rsq.col(k) = arma_rsq_z;
        
      }
      
    }
    
    // calculate the exact Shapley value for r-squared
    for (int j = 0; j < n_cols_x; j++) {
      
      arma::uvec arma_ix_pos = find(arma_ix.row(j));
      arma::uvec arma_ix_neg = find(arma_ix.row(j) == 0);
      arma::ivec arma_ix_n = arma_n(arma_ix_neg);
      arma::mat arma_rsq_diff = arma_rsq.cols(arma_ix_pos) - arma_rsq.cols(arma_ix_neg);
      
      for (int k = 0; k < n_combn / 2; k++) {
        
        int s = arma_ix_n[k];
        long double weight = (factorial(s) * factorial(n_cols_x - s - 1)) /
          (long double)factorial(n_cols_x);
        
        arma_rsq_sum.col(j) += weight * arma_rsq_diff.col(k);
        
      }
      
    }
    
    // create and return a matrix or xts object for Shapley values
    NumericMatrix rsq(wrap(arma_rsq_sum));
    List dimnames_x = xx.attr("dimnames");
    rsq.attr("dimnames") = dimnames_x;
    rsq.attr("index") = xx.attr("index");
    rsq.attr(".indexCLASS") = xx.attr(".indexCLASS");
    rsq.attr(".indexTZ") = xx.attr(".indexTZ");
    rsq.attr("tclass") = xx.attr("tclass");
    rsq.attr("tzone") = xx.attr("tzone");
    rsq.attr("class") = xx.attr("class");
    
    return rsq;
    
  } else if (Rf_isMatrix(y)) {
    
    NumericVector xx(x);
    NumericMatrix yy(y);
    NumericMatrix xxx(xx.size(), 1, xx.begin()); // consistent with roll's roll_lm
    
    int n_rows_xy = xxx.nrow();
    int n_cols_x = xxx.ncol();
    int n_cols_y = yy.ncol();
    List result_rsq(n_cols_y);
    List result_z(3);
    
    // create a list of matrices,
    // otherwise a list of lists
    if (n_cols_y == 1) {
      
      result_z = roll_lm_z(xxx, yy(_, 0), width,
                           weights, intercept, min_obs,
                           complete_obs, na_restore,
                           online);
      
      arma::vec arma_rsq_z = result_z[1];
      
      // create and return a vector object for Shapley values
      NumericVector rsq(wrap(arma_rsq_z));
      rsq.attr("dim") = IntegerVector::create(n_rows_xy, 1);
      List dimnames_x = xx.attr("dimnames");
      rsq.attr("dimnames") = dimnames_x;
      rsq.attr("index") = yy.attr("index");
      rsq.attr(".indexCLASS") = yy.attr(".indexCLASS");
      rsq.attr(".indexTZ") = yy.attr(".indexTZ");
      rsq.attr("tclass") = yy.attr("tclass");
      rsq.attr("tzone") = yy.attr("tzone");
      rsq.attr("class") = yy.attr("class");
      
      return rsq;
      
    } else {
      
      for (int z = 0; z < n_cols_y; z++) {
        
        result_z = roll_lm_z(xxx, yy(_, z), width,
                             weights, intercept, min_obs,
                             complete_obs, na_restore,
                             online);
        
        arma::vec arma_rsq_z = result_z[1];
        
        // create and return a vector object for Shapley values
        NumericVector rsq(wrap(arma_rsq_z));
        rsq.attr("dim") = IntegerVector::create(n_rows_xy, 1);
        List dimnames_x = xx.attr("dimnames");
        rsq.attr("dimnames") = dimnames_x;
        rsq.attr("index") = yy.attr("index");
        rsq.attr(".indexCLASS") = yy.attr(".indexCLASS");
        rsq.attr(".indexTZ") = yy.attr(".indexTZ");
        rsq.attr("tclass") = yy.attr("tclass");
        rsq.attr("tzone") = yy.attr("tzone");
        rsq.attr("class") = yy.attr("class");
        
        result_rsq(z) = rsq;
        
      }
      
      // add names to each list
      List dimnames_y = yy.attr("dimnames");
      result_rsq.attr("names") = dimnames_lm_y(dimnames_y, n_cols_y);
      
      return result_rsq;
      
    }
    
  } else {
    
    NumericVector xx(x);
    NumericVector yy(y);
    
    int n_rows_xy = xx.size();
    int n_cols_x = 1;
    int n_cols_y = 1;
    List result_rsq(n_cols_y);
    List result_z(3);
    
    // create a list of matrices
    result_z = roll_lm_z(xx, yy, width,
                         weights, intercept, min_obs,
                         complete_obs, na_restore,
                         online);
    
    arma::vec arma_rsq_z = result_z[1];
    
    if (intercept) {
      
      // create and return a vector object for Shapley values
      NumericVector rsq(wrap(arma_rsq_z));
      rsq.attr("dim") = IntegerVector::create(n_rows_xy, 1);
      List dimnames_x = xx.attr("dimnames");
      rsq.attr("dimnames") = dimnames_x;
      rsq.attr("index") = xx.attr("index");
      rsq.attr(".indexCLASS") = xx.attr(".indexCLASS");
      rsq.attr(".indexTZ") = xx.attr(".indexTZ");
      rsq.attr("tclass") = xx.attr("tclass");
      rsq.attr("tzone") = xx.attr("tzone");
      rsq.attr("class") = xx.attr("class");
      
      return rsq;
      
    } else {
      
      // create and return a vector object for Shapley values
      NumericVector rsq(wrap(arma_rsq_z));
      List names = xx.attr("names");
      rsq.attr("dim") = R_NilValue;
      if (names.size() > 0) {
        rsq.attr("names") = names;
      }
      rsq.attr("index") = xx.attr("index");
      rsq.attr("class") = xx.attr("class");
      
      return rsq;
      
    }
    
  }
  
}