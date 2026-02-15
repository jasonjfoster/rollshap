#ifndef ROLLSHAP_CHECK_H
#define ROLLSHAP_CHECK_H

#include <RcppArmadillo.h>
using namespace Rcpp;

// scalar checks: bounded integer [lower, upper]
inline void check_bounds_int(const int& value, const int& lower,
                      const int& upper, const char* name) {

  if ((value < lower) || (value > upper)) {
    stop("value of '%s' must be between %d and %d", name, lower, upper);
  }

}

#endif