// roadrunner -- src/krls.cpp
//
// Kernel Regularized Least Squares (Hainmueller & Hazlett 2014) engine.
// Mirrors the algorithm of KRLS::krls (KRLS 1.1.0) so fits agree to
// numerical precision when sigma + lambda are matched.
//
// Exports:
//   krls_kernel_cpp        -- TBB-parallel Gaussian kernel build (n x n).
//   krls_kernel_pred_cpp   -- TBB-parallel test-vs-train kernel (m x n).
//   krls_eig_cpp           -- LAPACK-backed symmetric eigendecomposition,
//                             returned in DESCENDING order to match
//                             base::eigen() conventions used by KRLS.
//   krls_solve_cpp         -- closed-form coefs / LOOe / diag(G^-1)
//                             from precomputed (V, Vsq, V'y, d).
//   krls_loo_loss_cpp      -- LOOe only (fast inner objective for the
//                             golden-section lambda search).
//   krls_deriv_cpp         -- vectorised pointwise marginal effects.
//   krls_avg_deriv_var_cpp -- variance of average marginal effects,
//                             avoiding the explicit n x n L matrix.
//
// Conventions:
//   - X is assumed already standardised (column-centred, unit sd) before
//     calling these functions.
//   - sigma is the bandwidth parameter such that K_ij =
//     exp(-||x_i - x_j||^2 / sigma).  Default sigma = ncol(X) (KRLS).
//   - All routines are single-precision-free; the linear algebra is
//     dispatched through Armadillo (LAPACK/BLAS) for the cubic-cost
//     pieces (eigendecomposition, V'y, V * alpha).

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <cmath>
#include <cstddef>

using namespace Rcpp;
using namespace RcppParallel;

// -------------------------------------------------------------------
// Gaussian-kernel build (train-train, n x n symmetric).
// -------------------------------------------------------------------
struct KrlsKernelWorker : public Worker {
  const RMatrix<double> X;
  RMatrix<double> K;
  const double sigma;

  KrlsKernelWorker(const NumericMatrix& X_, NumericMatrix& K_, double sigma_)
    : X(X_), K(K_), sigma(sigma_) {}

  void operator()(std::size_t begin, std::size_t end) {
    const std::size_t n = X.nrow();
    const std::size_t p = X.ncol();
    const double inv_sigma = 1.0 / sigma;
    for (std::size_t i = begin; i < end; ++i) {
      K(i, i) = 1.0;
      for (std::size_t j = i + 1; j < n; ++j) {
        double s = 0.0;
        for (std::size_t k = 0; k < p; ++k) {
          const double d = X(i, k) - X(j, k);
          s += d * d;
        }
        const double v = std::exp(-s * inv_sigma);
        K(i, j) = v;
        K(j, i) = v;
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix krls_kernel_cpp(const NumericMatrix& X, double sigma) {
  const int n = X.nrow();
  NumericMatrix K(n, n);
  KrlsKernelWorker w(X, K, sigma);
  parallelFor(0, n, w);
  return K;
}

// -------------------------------------------------------------------
// Test-vs-train Gaussian kernel (m_new x n_train).
// -------------------------------------------------------------------
struct KrlsKernelPredWorker : public Worker {
  const RMatrix<double> Xnew;
  const RMatrix<double> Xtr;
  RMatrix<double> K;
  const double sigma;

  KrlsKernelPredWorker(const NumericMatrix& Xnew_, const NumericMatrix& Xtr_,
                       NumericMatrix& K_, double sigma_)
    : Xnew(Xnew_), Xtr(Xtr_), K(K_), sigma(sigma_) {}

  void operator()(std::size_t begin, std::size_t end) {
    const std::size_t n = Xtr.nrow();
    const std::size_t p = Xtr.ncol();
    const double inv_sigma = 1.0 / sigma;
    for (std::size_t i = begin; i < end; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        double s = 0.0;
        for (std::size_t k = 0; k < p; ++k) {
          const double d = Xnew(i, k) - Xtr(j, k);
          s += d * d;
        }
        K(i, j) = std::exp(-s * inv_sigma);
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix krls_kernel_pred_cpp(const NumericMatrix& Xnew,
                                   const NumericMatrix& Xtrain,
                                   double sigma) {
  const int m = Xnew.nrow();
  const int n = Xtrain.nrow();
  NumericMatrix K(m, n);
  KrlsKernelPredWorker w(Xnew, Xtrain, K, sigma);
  parallelFor(0, m, w);
  return K;
}

// -------------------------------------------------------------------
// Symmetric eigendecomposition, descending order (matches base::eigen).
// -------------------------------------------------------------------
// [[Rcpp::export]]
List krls_eig_cpp(const arma::mat& K) {
  arma::vec vals;
  arma::mat vecs;
  // dc divide-and-conquer is fast + accurate for dense symmetric.
  const bool ok = arma::eig_sym(vals, vecs, K, "dc");
  if (!ok) {
    Rcpp::stop("krls_eig_cpp: eig_sym failed to converge");
  }
  // eig_sym returns ascending; flip to match KRLS (uses base::eigen which
  // returns descending eigenvalues).
  vals = arma::flipud(vals);
  vecs = arma::fliplr(vecs);
  return List::create(_["values"]  = vals,
                      _["vectors"] = vecs);
}

// -------------------------------------------------------------------
// Closed-form ridge solve given a precomputed eigendecomposition.
//
// Inputs (all already in standardised space):
//   d    : eigenvalues, length n.
//   V    : eigenvectors, n x n  (col i = eigenvector for d_i).
//   Vsq  : V .* V (n x n) -- precompute once outside the lambda loop.
//   Vty  : V' y, length n.
//   y2   : ||y||^2 (only used when caller wants it; unused here).
//   lambda : ridge penalty.
//
// Returns coeffs c = G^{-1} y where G = K + lambda I, the diagonal of
// G^{-1}, and the KRLS LOO error sum.
//
//   c     = V diag(1/(d+lambda)) V' y
//   diagG = sum_k V_{ik}^2 / (d_k + lambda)
//   Le    = sum_i (c_i / diagG_i)^2
// -------------------------------------------------------------------
// [[Rcpp::export]]
List krls_solve_cpp(const arma::vec& d, const arma::mat& V,
                    const arma::mat& Vsq, const arma::vec& Vty,
                    double lambda) {
  arma::vec inv = 1.0 / (d + lambda);
  arma::vec c   = V * (Vty % inv);
  arma::vec g   = Vsq * inv;
  arma::vec r   = c / g;
  return List::create(_["coeffs"]   = c,
                      _["Le"]       = arma::dot(r, r),
                      _["diagGinv"] = g);
}

// LOO loss only (avoids list allocation on the hot path).
// [[Rcpp::export]]
double krls_loo_loss_cpp(const arma::vec& d, const arma::mat& V,
                         const arma::mat& Vsq, const arma::vec& Vty,
                         double lambda) {
  arma::vec inv = 1.0 / (d + lambda);
  arma::vec c   = V * (Vty % inv);
  arma::vec g   = Vsq * inv;
  arma::vec r   = c / g;
  return arma::dot(r, r);
}

// -------------------------------------------------------------------
// Pointwise marginal effects.
//
//   derivmat[i, k] = (-2 / sigma) * sum_j (X_ik - X_jk) * K_ij * c_j
//                  = (-2 / sigma) * ( X_ik * (K c)_i - (K diag(c) X)_ik )
//
// All operations are matrix-multiplies dispatched to BLAS; no explicit
// per-variable n x n L matrix is constructed.
// -------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat krls_deriv_cpp(const arma::mat& X, const arma::mat& K,
                         const arma::vec& c, double sigma) {
  const arma::vec s1 = K * c;                  // n x 1, equals yhat_std.
  arma::mat Xc       = X;                      // n x p, will be scaled.
  Xc.each_col() %= c;                          // row j scaled by c_j.
  const arma::mat s2 = K * Xc;                 // n x p.
  arma::mat Xs1      = X;
  Xs1.each_col() %= s1;                        // row i scaled by s1_i.
  return (-2.0 / sigma) * (Xs1 - s2);
}

// -------------------------------------------------------------------
// Variance of the AVERAGE marginal effect, per variable.
//
// KRLS computes
//   var_avg_k = (1 / n^2) * (4 / sigma^2) * sum( t(L_k) %*% V_c %*% L_k )
// with L_k = distk_k * K  (n x n) and V_c the coefficient covariance
//   V_c = sigma2 * V diag((d + lambda)^-2) V'.
//
// Because sum(t(A) M A) = 1' t(A) M A 1 = (A 1)' M (A 1), only the row
// sum u_k = L_k * 1 is needed.  In closed form
//   u_{ik} = X_ik * (K * 1)_i - (K * X)_{ik}.
// Then
//   var_avg_k = (4 / sigma^2 / n^2) * u_k' V_c u_k
//             = (4 / sigma^2 / n^2) * sum_i w_i * (V' u_k)_{i}^2,
// where w = sigma2 * (d + lambda)^-2.
// -------------------------------------------------------------------
// [[Rcpp::export]]
arma::vec krls_avg_deriv_var_cpp(const arma::mat& X, const arma::mat& K,
                                 const arma::mat& V, const arma::vec& d,
                                 double sigma, double lambda,
                                 double sigma2) {
  const arma::uword n = X.n_rows;
  const arma::vec Krow = K * arma::ones<arma::vec>(n);   // n x 1.
  const arma::mat KX   = K * X;                          // n x p.
  arma::mat U          = X;                              // n x p (column scratch).
  U.each_col() %= Krow;                                  // X_ik * Krow_i.
  U -= KX;                                               // u_{ik}.
  const arma::mat VtU  = V.t() * U;                      // n x p.
  const arma::vec w    = sigma2 / arma::square(d + lambda); // n x 1.
  arma::mat VtU2       = VtU % VtU;                      // n x p.
  // sum_i w_i * VtU2_{ik} for each column k = (VtU2' * w)_k.
  arma::vec out = VtU2.t() * w;                          // p x 1.
  const double scale = 4.0 / (sigma * sigma) / static_cast<double>(n * n);
  return scale * out;
}

// -------------------------------------------------------------------
// Convenience: V .* V (used to compute diag(G^{-1}) at each lambda).
// Exposed so R-side can precompute once outside the lambda search loop.
// -------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat krls_vsq_cpp(const arma::mat& V) {
  return V % V;
}

// ---------------------------------------------------------------------------
// Phase 1 (v0.0.0.9042): shared pairwise squared-distance matrix.
//
// D[i, j] = ||X_a[i, ] - X_b[j, ]||^2
//
// Identity: D = ||X_a||^2 (row sums) + ||X_b||^2^T - 2 X_a X_b^T
// One dgemm dominates; tiny negatives clamped to 0 to absorb FP noise.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
arma::mat krls_pairwise_sqdist_cpp(const arma::mat& X_a, const arma::mat& X_b) {
  if (X_a.n_cols != X_b.n_cols) {
    Rcpp::stop("krls_pairwise_sqdist_cpp: X_a and X_b must have same n_cols");
  }
  arma::vec ra = arma::sum(arma::square(X_a), 1);   // n_a x 1
  arma::vec rb = arma::sum(arma::square(X_b), 1);   // n_b x 1
  arma::mat D = -2.0 * (X_a * X_b.t());             // n_a x n_b
  D.each_col() += ra;
  D.each_row() += rb.t();
  D.elem(arma::find(D < 0.0)).zeros();              // FP-noise floor
  return D;
}

// ---------------------------------------------------------------------------
// Phase 1 (v0.0.0.9042): inner autotune loop in C++.
//
// For each sigma s in sigma_grid:
//   K_tr = exp(-D_tr / s)
//   K_te = exp(-D_te / s)
//   (V, d) = eig_sym(K_tr)
//   lambda = golden-section LOO search (tol=1e-6, log-step L)
//   alpha  = V * ((Vty / (d + lambda)))    (closed-form ridge in eigenbasis)
//   mse_s  = mean((y_te - K_te * alpha)^2)
//
// THIS TASK: sequential only (nthreads ignored). Task 4 adds parallelFor.
// ---------------------------------------------------------------------------

static double krls_loo_loss_inline(const arma::vec& d, const arma::mat& V,
                                   const arma::mat& Vsq, const arma::vec& Vty,
                                   double lambda) {
  arma::vec inv = 1.0 / (d + lambda);
  arma::vec c   = V * (Vty % inv);
  arma::vec g   = Vsq * inv;
  return arma::dot(c / g, c / g);
}

static void krls_one_sigma(const arma::mat& D_tr, const arma::mat& D_te,
                           const arma::vec& y_tr, const arma::vec& y_te,
                           double s,
                           double lam_tol, double L0, double L_step,
                           double& mse_out, double& lam_out) {
  arma::mat K_tr = arma::exp(-D_tr / s);
  arma::mat K_te = arma::exp(-D_te / s);

  arma::vec d;
  arma::mat V;
  arma::eig_sym(d, V, K_tr);

  arma::vec Vty = V.t() * y_tr;
  arma::mat Vsq = V % V;

  // q = which.min(|d - max(d)/1000|), as 1-based R index in the
  // DESCENDING-ordered eigenvalue vector (matches krls_eig_cpp + R reference).
  // arma::eig_sym returns ASCENDING, so an ascending 0-based index q_idx
  // maps to descending 1-based position (n - q_idx).
  double dmax  = d.max();
  double target = dmax / 1000.0;
  arma::uword q_idx = 0;
  double mind = std::abs(d(0) - target);
  for (arma::uword i = 1; i < d.n_elem; ++i) {
    double dd = std::abs(d(i) - target);
    if (dd < mind) { mind = dd; q_idx = i; }
  }
  double q_R = (double)(d.n_elem - q_idx);  // descending 1-based index

  // Lower bracket: grow L until sum(d/(d+L)) <= q_R
  double L = L0;
  while (arma::sum(d / (d + L)) > q_R) {
    L *= L_step;
    if (L > 1e30) break;
  }

  // Upper bracket: shrink U from n until sum(d/(d+U)) >= 1
  double U = (double) d.n_elem;
  while (arma::sum(d / (d + U)) < 1.0) {
    U -= 1.0;
    if (U <= L) { U = L * 10.0; break; }
  }

  // Golden section
  const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
  double X1 = L + (1.0 - gr) * (U - L);
  double X2 = L + gr * (U - L);
  double S1 = krls_loo_loss_inline(d, V, Vsq, Vty, X1);
  double S2 = krls_loo_loss_inline(d, V, Vsq, Vty, X2);
  int iter = 0;
  while (std::abs(S1 - S2) > lam_tol && iter < 200) {
    if (S1 < S2) {
      U = X2; X2 = X1; X1 = L + (1.0 - gr) * (U - L);
      S2 = S1; S1 = krls_loo_loss_inline(d, V, Vsq, Vty, X1);
    } else {
      L = X1; X1 = X2; X2 = L + gr * (U - L);
      S1 = S2; S2 = krls_loo_loss_inline(d, V, Vsq, Vty, X2);
    }
    ++iter;
  }
  double lam = 0.5 * (X1 + X2);

  // Closed-form ridge solution in eigenbasis: alpha = V * (Vty / (d + lam))
  arma::vec alpha = V * (Vty / (d + lam));
  arma::vec yhat  = K_te * alpha;
  mse_out = arma::mean(arma::square(y_te - yhat));
  lam_out = lam;
}

// [[Rcpp::export]]
Rcpp::List krls_autotune_inner_cpp(const arma::mat& D_tr, const arma::mat& D_te,
                                   const arma::vec& y_tr, const arma::vec& y_te,
                                   const arma::vec& sigma_grid,
                                   Rcpp::List lambda_args,
                                   int nthreads) {
  const arma::uword nsigma = sigma_grid.n_elem;
  arma::vec mse(nsigma, arma::fill::zeros);
  arma::vec lam(nsigma, arma::fill::zeros);

  double lam_tol = Rcpp::as<double>(lambda_args["tol"]);
  double L0      = Rcpp::as<double>(lambda_args["L0"]);
  double L_step  = Rcpp::as<double>(lambda_args["L_step"]);

  // Sequential loop. Task 4 swaps in parallelFor.
  for (arma::uword i = 0; i < nsigma; ++i) {
    double m = 0.0, l = 0.0;
    krls_one_sigma(D_tr, D_te, y_tr, y_te, sigma_grid(i),
                   lam_tol, L0, L_step, m, l);
    mse(i) = m;
    lam(i) = l;
  }

  // Wrap as plain NumericVector (no dim attribute) so tests comparing
  // against R-side numeric vectors pass without dim mismatches.
  Rcpp::NumericVector mse_out(nsigma), lam_out(nsigma);
  for (arma::uword i = 0; i < nsigma; ++i) {
    mse_out[i] = mse(i);
    lam_out[i] = lam(i);
  }

  return Rcpp::List::create(
    Rcpp::Named("mse_per_sigma")    = mse_out,
    Rcpp::Named("lambda_per_sigma") = lam_out,
    Rcpp::Named("nthreads_used")    = 1
  );
}
