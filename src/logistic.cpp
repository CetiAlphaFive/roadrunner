// roadrunner -- src/logistic.cpp
//
// Binary logistic regression engine (iteratively reweighted least squares).
//
// Exports:
//   logreg_fit_cpp  -- IRLS / Fisher scoring solve. Returns coefficients,
//                      fitted probabilities, linear predictor, IRLS working
//                      weights, deviance, null deviance, iteration count,
//                      a convergence flag, (X'WX)^-1, and the hat-diagonal.
//                      Non-convergence (e.g. perfect separation) is NOT an
//                      error -- it returns converged = false and the R side
//                      warns.
//   logreg_vcov_cpp -- classical and HC0/HC1/HC2/HC3 sandwich
//                      variance-covariance matrices.
//
// Conventions:
//   - Each IRLS step solves a weighted least squares problem on the
//     working response z = eta + (y - mu) / (mu (1 - mu)) with working
//     weights W = w .* mu (1 - mu). The step is solved via an economy QR
//     of diag(sqrt(W)) X, mirroring the OLS engine.
//   - Fitted probabilities mu are clamped to [eps, 1 - eps] so the
//     working weights and the deviance stay finite under separation.
//   - All routines are double precision; the cubic-cost piece (the QR
//     factorisation) is dispatched through Armadillo (LAPACK/BLAS).
//   - Determinism: every routine is a single deterministic linear-algebra
//     pass with no threading, so output is byte-identical run to run and
//     across thread counts.

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

namespace {

// Numerically safe logistic CDF: mu = 1 / (1 + exp(-eta)), evaluated in a
// branchy form that avoids overflow for large |eta|.
inline double plogis_one(double eta) {
  if (eta >= 0.0) {
    const double z = std::exp(-eta);
    return 1.0 / (1.0 + z);
  }
  const double z = std::exp(eta);
  return z / (1.0 + z);
}

// Binomial deviance contribution clamped against log(0): the caller passes
// mu already clamped to [eps, 1 - eps], so each log is finite.
inline double dev_resid_sq(double yi, double mui, double wi) {
  double term = 0.0;
  if (yi > 0.0)        term += yi * std::log(yi / mui);
  if (yi < 1.0)        term += (1.0 - yi) * std::log((1.0 - yi) / (1.0 - mui));
  return 2.0 * wi * term;
}

}  // namespace

// Binary logistic regression via iteratively reweighted least squares.
//
// X     -- n x p design matrix (intercept column included by the caller).
// y     -- length-n response, each element in {0, 1}.
// w     -- length-n strictly positive prior weights (ones for the
//          unweighted fit).
// maxit -- maximum number of IRLS iterations.
// tol   -- convergence tolerance on the relative change in deviance.
//
// Fisher scoring: at each step form eta = X beta, mu = plogis(eta) clamped
// to [eps, 1 - eps], working weights W = w .* mu (1 - mu) and working
// response z = eta + (y - mu) / (mu (1 - mu)); update beta by a weighted
// economy-QR least-squares solve of diag(sqrt(W)) X against
// diag(sqrt(W)) z. Iterate until |dev - dev_old| / (|dev| + 0.1) < tol or
// maxit is reached. Non-convergence is reported via converged = false.
//
// [[Rcpp::export]]
Rcpp::List logreg_fit_cpp(const arma::mat& X, const arma::vec& y,
                          const arma::vec& w, int maxit, double tol) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  if (y.n_elem != n)
    Rcpp::stop("logreg_fit_cpp: length(y) must equal nrow(X).");
  if (w.n_elem != n)
    Rcpp::stop("logreg_fit_cpp: length(w) must equal nrow(X).");
  if (n <= p)
    Rcpp::stop("logreg_fit_cpp: need more rows than columns "
               "(n = %u, p = %u).", n, p);
  if (maxit < 1)
    Rcpp::stop("logreg_fit_cpp: maxit must be at least 1.");

  const double eps = 1e-10;          // probability clamp.
  const double wsum = arma::accu(w);

  // Null deviance: the intercept-only model fits the constant probability
  // mu0 = weighted mean of y (clamped). Used as the deviance baseline.
  double mu0 = arma::dot(w, y) / wsum;
  mu0 = std::min(std::max(mu0, eps), 1.0 - eps);
  double null_deviance = 0.0;
  for (arma::uword i = 0; i < n; ++i)
    null_deviance += dev_resid_sq(y(i), mu0, w(i));

  // IRLS state. beta starts at zero -> eta == 0, mu == 0.5 everywhere.
  arma::vec beta(p, arma::fill::zeros);
  arma::vec eta(n, arma::fill::zeros);
  arma::vec mu(n);
  mu.fill(0.5);
  arma::vec Wirls(n);                // working weights w .* mu (1 - mu).

  double deviance = arma::datum::inf;
  bool converged = false;
  int iter = 0;
  arma::mat R;                       // last accepted QR R factor (p x p).
  arma::mat Q;                       // last accepted QR Q factor (n x p).

  for (iter = 1; iter <= maxit; ++iter) {
    // Working weights and working response at the current iterate.
    const arma::vec varmu = mu % (1.0 - mu);
    Wirls = w % varmu;
    const arma::vec z = eta + (y - mu) / varmu;

    // Weighted economy-QR solve: diag(sqrt(W)) X beta = diag(sqrt(W)) z.
    const arma::vec sw = arma::sqrt(Wirls);
    arma::mat Xw = X;
    Xw.each_col() %= sw;
    const arma::vec zw = z % sw;

    arma::mat Qi, Ri;
    if (!arma::qr_econ(Qi, Ri, Xw))
      Rcpp::stop("logreg_fit_cpp: QR factorisation failed at iteration %d.",
                 iter);

    // Rank check on the working design. A near-zero pivot signals a
    // collinear column (independent of separation); a hard error.
    const arma::vec rdiag = arma::abs(Ri.diag());
    const double rmax = rdiag.max();
    const double rtol = rmax * std::sqrt(arma::datum::eps) *
                        static_cast<double>(std::max<arma::uword>(n, p));
    for (arma::uword j = 0; j < p; ++j)
      if (rdiag(j) <= rtol)
        Rcpp::stop("logreg_fit_cpp: rank-deficient model matrix; "
                   "drop collinear columns.");

    beta = arma::solve(arma::trimatu(Ri), Qi.t() * zw);
    Q = Qi;
    R = Ri;

    // Refresh eta / mu / deviance at the new beta.
    eta = X * beta;
    for (arma::uword i = 0; i < n; ++i)
      mu(i) = std::min(std::max(plogis_one(eta(i)), eps), 1.0 - eps);
    double dev_new = 0.0;
    for (arma::uword i = 0; i < n; ++i)
      dev_new += dev_resid_sq(y(i), mu(i), w(i));

    if (std::abs(dev_new - deviance) / (std::abs(dev_new) + 0.1) < tol) {
      deviance = dev_new;
      converged = true;
      break;
    }
    deviance = dev_new;
  }
  if (iter > maxit) iter = maxit;

  // Final Fisher information at the converged iterate. R and Q above come
  // from the QR of the working design BEFORE the last beta update; the
  // observed covariance and hat values must use the working weights at
  // the converged mu, so re-factor diag(sqrt(W_final)) X. This matches
  // the (X'WX)^-1 that glm() reports from its last iteration.
  {
    const arma::vec varmu = mu % (1.0 - mu);
    const arma::vec sw = arma::sqrt(w % varmu);
    arma::mat Xw = X;
    Xw.each_col() %= sw;
    arma::mat Qf, Rf;
    if (!arma::qr_econ(Qf, Rf, Xw))
      Rcpp::stop("logreg_fit_cpp: final QR factorisation failed.");
    Q = Qf;
    R = Rf;
  }

  // (X'WX)^-1 from the final QR: R^-1 R^-T. R is upper-triangular.
  const arma::mat Rinv = arma::inv(arma::trimatu(R));
  const arma::mat XtWXinv = Rinv * Rinv.t();

  // Hat diagonal of the final weighted problem: h_ii = rowSums(Q^2).
  const arma::vec hatdiag = arma::sum(arma::square(Q), 1);

  return Rcpp::List::create(
    Rcpp::Named("coefficients")  = beta,
    Rcpp::Named("fitted")        = mu,
    Rcpp::Named("eta")           = eta,
    Rcpp::Named("working")       = Wirls,
    Rcpp::Named("deviance")      = deviance,
    Rcpp::Named("null.deviance") = null_deviance,
    Rcpp::Named("iter")          = iter,
    Rcpp::Named("converged")     = converged,
    Rcpp::Named("XtWXinv")       = XtWXinv,
    Rcpp::Named("hatdiag")       = hatdiag);
}

// Variance-covariance matrix for a logistic-regression fit.
//
// X       -- n x p design matrix.
// w       -- length-n prior weights.
// y       -- length-n response in {0, 1}.
// mu      -- length-n fitted probabilities (from logreg_fit_cpp).
// XtWXinv -- (X'WX)^-1 (the bread; from logreg_fit_cpp).
// hatdiag -- hat diagonal of the final IRLS problem (from logreg_fit_cpp).
// type    -- 0 classical, 1 HC0, 2 HC1, 3 HC2, 4 HC3.
//
// Classical: (X'WX)^-1, where W = w .* mu (1 - mu) is the Fisher
// information weight; this is the maximum-likelihood asymptotic covariance.
// Sandwich:  B * M * B with bread B = (X'WX)^-1 and meat
//   M = sum_i (X_i X_i') s_i^2 c_i,
// where s_i = w_i (y_i - mu_i) is the score residual and the leverage
// correction c_i is
//   HC0: 1                HC1: n/(n-p)
//   HC2: 1/(1-h_ii)       HC3: 1/(1-h_ii)^2.
//
// [[Rcpp::export]]
arma::mat logreg_vcov_cpp(const arma::mat& X, const arma::vec& w,
                          const arma::vec& y, const arma::vec& mu,
                          const arma::mat& XtWXinv, const arma::vec& hatdiag,
                          int type) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  if (type == 0)
    return XtWXinv;

  // Score residuals s_i = w_i (y_i - mu_i). The meat is built directly on
  // the (unweighted-design) score contributions so HC corrections match
  // sandwich::vcovHC on a glm fit.
  const arma::vec score = w % (y - mu);

  arma::vec cfac(n, arma::fill::ones);
  if (type == 2) {                          // HC1
    cfac.fill(static_cast<double>(n) /
              static_cast<double>(n - p));
  } else if (type == 3) {                   // HC2
    cfac = 1.0 / (1.0 - hatdiag);
  } else if (type == 4) {                   // HC3
    cfac = 1.0 / arma::square(1.0 - hatdiag);
  } else if (type != 1) {                   // HC0 leaves cfac == 1
    Rcpp::stop("logreg_vcov_cpp: unknown type code %d.", type);
  }

  // Diagonal meat weights: omega_i = s_i^2 * c_i. The meat is then
  // X' diag(omega) X.
  const arma::vec omega = arma::square(score) % cfac;
  arma::mat Xom = X;
  Xom.each_col() %= omega;                  // rowwise omega scaling.
  const arma::mat meat = X.t() * Xom;       // p x p.

  return XtWXinv * meat * XtWXinv;
}
