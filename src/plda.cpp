// src/plda.cpp
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Within-class standard deviation per feature, divisor n. C++-internal (returns arma::vec).
static arma::vec wcsd_impl(const arma::mat& x, const arma::ivec& y, int G) {
  const arma::uword n = x.n_rows, p = x.n_cols;
  arma::vec ss(p, arma::fill::zeros);
  for (int g = 1; g <= G; ++g) {
    arma::uvec idx = arma::find(y == g);
    if (idx.n_elem == 0) continue;
    arma::mat xg = x.rows(idx);
    arma::rowvec mg = arma::mean(xg, 0);
    xg.each_row() -= mg;
    ss += arma::sum(arma::square(xg), 0).t();
  }
  return arma::sqrt(ss / (double) n);
}

// [[Rcpp::export]]
Rcpp::NumericVector plda_wcsd_cpp(const arma::mat& x, const arma::ivec& y, int G) {
  arma::vec v = wcsd_impl(x, y, G);
  return Rcpp::NumericVector(v.begin(), v.end());
}

// Element-wise soft-threshold: sign(u) * max(|u| - lam, 0).
static inline arma::vec soft_threshold(const arma::vec& u, double lam) {
  return arma::sign(u) % arma::clamp(arma::abs(u) - lam, 0.0, arma::datum::inf);
}

// L2-normalize; returns the zero vector unchanged.
static inline arma::vec normalize_l2(const arma::vec& v) {
  double nv = arma::norm(v, 2);
  return (nv > 0.0) ? arma::vec(v / nv) : v;
}

// [[Rcpp::export]]
Rcpp::NumericVector plda_softthresh_cpp(const arma::vec& u, double lam) {
  arma::vec v = soft_threshold(u, lam);
  return Rcpp::NumericVector(v.begin(), v.end());
}
