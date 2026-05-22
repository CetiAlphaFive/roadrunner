// roadrunner -- src/ols.cpp
//
// Ordinary / weighted least squares engine.
//
// Exports:
//   ols_fit_cpp  -- weighted economy-QR solve. Returns coefficients,
//                   fitted values, residuals, sigma2, residual df,
//                   (X'WX)^-1, hat-diagonal, and the numerical rank.
//                   Errors clearly on a rank-deficient design matrix.
//   ols_vcov_cpp -- classical and HC0/HC1/HC2/HC3 sandwich
//                   variance-covariance matrices.
//
// Conventions:
//   - The weighted problem is solved on the transformed design
//     Xw = sqrt(w) .* X (rowwise), yw = sqrt(w) .* y. With w == 1 this
//     reduces to ordinary least squares.
//   - All routines are double precision; the cubic-cost piece (the QR
//     factorisation) is dispatched through Armadillo (LAPACK/BLAS).
//   - Determinism: every routine is a single deterministic linear-algebra
//     pass with no threading, so output is byte-identical run to run.

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// Weighted economy-QR least-squares fit.
//
// X  -- n x p design matrix (intercept column included by the caller).
// y  -- length-n response.
// w  -- length-n strictly positive weights (pass a vector of ones for OLS).
//
// Forms Xw = diag(sqrt(w)) X, yw = diag(sqrt(w)) y; takes the economy QR
// Xw = Q R; solves coef = R^-1 Q' yw. Rank deficiency is detected from the
// magnitude of the smallest |diag(R)| relative to the largest; a deficient
// design is a hard error -- the caller is told to drop collinear columns.
//
// [[Rcpp::export]]
Rcpp::List ols_fit_cpp(const arma::mat& X, const arma::vec& y,
                       const arma::vec& w) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  if (y.n_elem != n)
    Rcpp::stop("ols_fit_cpp: length(y) must equal nrow(X).");
  if (w.n_elem != n)
    Rcpp::stop("ols_fit_cpp: length(w) must equal nrow(X).");
  if (n <= p)
    Rcpp::stop("ols_fit_cpp: need more rows than columns "
               "(n = %u, p = %u).", n, p);

  const arma::vec sw = arma::sqrt(w);
  arma::mat Xw = X;
  Xw.each_col() %= sw;                  // rowwise sqrt(w) scaling.
  const arma::vec yw = y % sw;

  arma::mat Q, R;
  if (!arma::qr_econ(Q, R, Xw))
    Rcpp::stop("ols_fit_cpp: QR factorisation failed.");

  // Rank detection from the diagonal of R. For an economy QR of an n x p
  // matrix R is p x p upper-triangular; a near-zero pivot signals an
  // (almost) collinear column.
  const arma::vec rdiag = arma::abs(R.diag());
  const double rmax = rdiag.max();
  const double rtol = rmax * std::sqrt(arma::datum::eps) *
                      static_cast<double>(std::max<arma::uword>(n, p));
  arma::uword rank = 0;
  for (arma::uword j = 0; j < p; ++j)
    if (rdiag(j) > rtol) ++rank;
  if (rank < p)
    Rcpp::stop("ols_fit_cpp: rank-deficient model matrix "
               "(rank %u < %u columns); drop collinear columns.",
               rank, p);

  // coef = R^-1 (Q' yw). R is upper-triangular -> a triangular solve.
  const arma::vec Qty = Q.t() * yw;
  const arma::vec coef = arma::solve(arma::trimatu(R), Qty);

  // Fitted values / residuals on the ORIGINAL scale.
  const arma::vec fitted = X * coef;
  const arma::vec resid  = y - fitted;

  // sigma2 = sum(w .* resid^2) / (n - p).
  const arma::uword df = n - p;
  const double rss = arma::dot(w, arma::square(resid));
  const double sigma2 = rss / static_cast<double>(df);

  // (X'WX)^-1 = R^-1 R^-T. Rinv is upper-triangular.
  const arma::mat Rinv = arma::inv(arma::trimatu(R));
  const arma::mat XtXinv = Rinv * Rinv.t();

  // Hat diagonal of the WEIGHTED problem: h_ii = rowSums(Q^2).
  const arma::vec hatdiag = arma::sum(arma::square(Q), 1);

  return Rcpp::List::create(
    Rcpp::Named("coefficients") = coef,
    Rcpp::Named("fitted")       = fitted,
    Rcpp::Named("residuals")    = resid,
    Rcpp::Named("sigma2")       = sigma2,
    Rcpp::Named("rss")          = rss,
    Rcpp::Named("df")           = static_cast<int>(df),
    Rcpp::Named("XtXinv")       = XtXinv,
    Rcpp::Named("hatdiag")      = hatdiag,
    Rcpp::Named("rank")         = static_cast<int>(rank));
}

// Variance-covariance matrix for a weighted-least-squares fit.
//
// X       -- n x p design matrix.
// w       -- length-n weights.
// resid   -- length-n residuals on the ORIGINAL scale (y - X coef).
// XtXinv  -- (X'WX)^-1 (the bread; from ols_fit_cpp).
// hatdiag -- hat diagonal of the weighted problem (from ols_fit_cpp).
// sigma2  -- residual variance (from ols_fit_cpp).
// type    -- 0 classical, 1 HC0, 2 HC1, 3 HC2, 4 HC3.
//
// Classical: sigma2 * (X'WX)^-1.
// Sandwich:  B * M * B with bread B = (X'WX)^-1 and meat
//   M = sum_i (Xw_i Xw_i') ew_i^2 c_i,
// where Xw_i = sqrt(w_i) X_i, ew_i = sqrt(w_i) resid_i is the residual of
// the transformed problem, and the leverage correction c_i is
//   HC0: 1                HC1: n/(n-p)
//   HC2: 1/(1-h_ii)       HC3: 1/(1-h_ii)^2.
//
// [[Rcpp::export]]
arma::mat ols_vcov_cpp(const arma::mat& X, const arma::vec& w,
                       const arma::vec& resid, const arma::mat& XtXinv,
                       const arma::vec& hatdiag, double sigma2,
                       int type) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  if (type == 0)
    return sigma2 * XtXinv;

  // Transformed design and residuals: the sandwich is built on the
  // sqrt(w)-weighted problem so the HC corrections match sandwich::vcovHC
  // on a weighted lm.
  const arma::vec sw = arma::sqrt(w);
  arma::mat Xw = X;
  Xw.each_col() %= sw;
  const arma::vec ew = resid % sw;          // transformed residuals.

  arma::vec cfac(n, arma::fill::ones);
  if (type == 2) {                          // HC1
    cfac.fill(static_cast<double>(n) /
              static_cast<double>(n - p));
  } else if (type == 3) {                   // HC2
    cfac = 1.0 / (1.0 - hatdiag);
  } else if (type == 4) {                   // HC3
    cfac = 1.0 / arma::square(1.0 - hatdiag);
  } else if (type != 1) {                   // HC0 leaves cfac == 1
    Rcpp::stop("ols_vcov_cpp: unknown type code %d.", type);
  }

  // Diagonal meat weights: omega_i = ew_i^2 * c_i. The meat is then
  // Xw' diag(omega) Xw.
  const arma::vec omega = arma::square(ew) % cfac;
  arma::mat Xom = Xw;
  Xom.each_col() %= omega;                  // rowwise omega scaling.
  const arma::mat meat = Xw.t() * Xom;      // p x p.

  return XtXinv * meat * XtXinv;
}
