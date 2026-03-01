// todo (roll >= 1.2.1): use namespace roll

#ifndef ROLLSHAP_ATTR_H
#define ROLLSHAP_ATTR_H

#include <RcppArmadillo.h>
using namespace Rcpp;

namespace rollshap {

// xts attributes: copy dimnames, index, .indexCLASS, .indexTZ, tclass, tzone, class
template <typename T, typename S>
inline void xts_attr(T& target, const S& source) {

  target.attr("dimnames") = source.attr("dimnames");
  target.attr("index") = source.attr("index");
  target.attr(".indexCLASS") = source.attr(".indexCLASS");
  target.attr(".indexTZ") = source.attr(".indexTZ");
  target.attr("tclass") = source.attr("tclass");
  target.attr("tzone") = source.attr("tzone");
  target.attr("class") = source.attr("class");

}

// xts attributes: copy xts attrs from source with supplied dimnames override
template <typename T, typename S>
inline void xts_attr(T& target, const S& source,
                     SEXP dimnames) {

  target.attr("dimnames") = dimnames;
  target.attr("index") = source.attr("index");
  target.attr(".indexCLASS") = source.attr(".indexCLASS");
  target.attr(".indexTZ") = source.attr(".indexTZ");
  target.attr("tclass") = source.attr("tclass");
  target.attr("tzone") = source.attr("tzone");
  target.attr("class") = source.attr("class");

}

// vector attributes: strip dim, conditionally copy names, copy index and class
template <typename T, typename S>
inline void vec_attr(T& target, const S& source) {

  target.attr("dim") = R_NilValue;
  List names = source.attr("names");

  if (names.size() > 0) {
    target.attr("names") = names;
  }

  target.attr("index") = source.attr("index");
  target.attr("class") = source.attr("class");

}

}

#endif