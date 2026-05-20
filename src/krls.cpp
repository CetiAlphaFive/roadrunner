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
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

// -------------------------------------------------------------------
// Kernel-type code (Phase Q2, v0.0.0.9049):
//   0 = Gaussian (default; back-compat path).
//   1 = Linear: K_ij = x_i' x_j.
//   2..5 = Polynomial degree 1..4: K_ij = (x_i' x_j + c)^d.
//
// The Gaussian path is preserved byte-for-byte (the existing
// constant-sigma fast loop and the ARD per-feature-divide loop both
// execute when kernel_type == 0).
// -------------------------------------------------------------------

// Compute (x_i . x_j) for the linear/poly kernels.
static inline double krls_dot_ij(const RMatrix<double>& X,
                                 std::size_t i, std::size_t j,
                                 std::size_t p) {
  double s = 0.0;
  for (std::size_t k = 0; k < p; ++k) {
    s += X(i, k) * X(j, k);
  }
  return s;
}
static inline double krls_dot_ij_xy(const RMatrix<double>& Xa,
                                    const RMatrix<double>& Xb,
                                    std::size_t i, std::size_t j,
                                    std::size_t p) {
  double s = 0.0;
  for (std::size_t k = 0; k < p; ++k) {
    s += Xa(i, k) * Xb(j, k);
  }
  return s;
}

// -------------------------------------------------------------------
// Gaussian-kernel build (train-train, n x n symmetric).
//
// Phase 2a (v0.0.0.9045): accepts a length-p `sigma_vec` of per-feature
// lengthscales. K_ij = exp(-sum_k (x_ik - x_jk)^2 / sigma_k).
//
// FP-determinism back-compat: when `sigma_vec` is bit-exact constant
// (every element == sigma_vec[0]), the worker dispatches to the
// EXACT v0.0.0.9044 inner loop (accumulate sum of squared d, then
// exp(-s * inv_sigma)) so the scalar-sigma path is byte-identical to
// pre-9045 fits. Non-constant sigma_vec uses per-element divides.
//
// Phase Q2 (v0.0.0.9049): polymorphic dispatch on kernel_type. When
// kernel_type == 0 (Gaussian) the existing fast loops are unchanged.
// kernel_type 1 = Linear; 2..5 = Poly1..4 with constant `kernel_c`.
// -------------------------------------------------------------------
struct KrlsKernelWorker : public Worker {
  const RMatrix<double> X;
  RMatrix<double> K;
  const std::vector<double> inv_sig;   // 1 / sigma_vec[k]
  const bool sigma_constant;
  const double inv_sigma_scalar;       // valid only when sigma_constant
  const int kernel_type;
  const double kernel_c;

  KrlsKernelWorker(const NumericMatrix& X_, NumericMatrix& K_,
                   const std::vector<double>& inv_sig_,
                   bool sigma_constant_, double inv_sigma_scalar_,
                   int kernel_type_, double kernel_c_)
    : X(X_), K(K_), inv_sig(inv_sig_),
      sigma_constant(sigma_constant_),
      inv_sigma_scalar(inv_sigma_scalar_),
      kernel_type(kernel_type_), kernel_c(kernel_c_) {}

  void operator()(std::size_t begin, std::size_t end) {
    const std::size_t n = X.nrow();
    const std::size_t p = X.ncol();
    if (kernel_type == 0) {
      if (sigma_constant) {
        // EXACT v0.0.0.9044 hot loop -- preserves FP byte-identity.
        const double inv_sigma = inv_sigma_scalar;
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
      } else {
        // ARD path -- per-feature divide.
        for (std::size_t i = begin; i < end; ++i) {
          K(i, i) = 1.0;
          for (std::size_t j = i + 1; j < n; ++j) {
            double s = 0.0;
            for (std::size_t k = 0; k < p; ++k) {
              const double d = X(i, k) - X(j, k);
              s += d * d * inv_sig[k];
            }
            const double v = std::exp(-s);
            K(i, j) = v;
            K(j, i) = v;
          }
        }
      }
    } else if (kernel_type == 1) {
      // Linear: K_ij = x_i . x_j (symmetric).
      for (std::size_t i = begin; i < end; ++i) {
        // diag: K_ii = x_i . x_i
        double si = 0.0;
        for (std::size_t k = 0; k < p; ++k) si += X(i, k) * X(i, k);
        K(i, i) = si;
        for (std::size_t j = i + 1; j < n; ++j) {
          const double v = krls_dot_ij(X, i, j, p);
          K(i, j) = v;
          K(j, i) = v;
        }
      }
    } else {
      // Polynomial 1..4: K_ij = (x_i . x_j + c)^d  with d = kernel_type - 1.
      const int deg = kernel_type - 1;
      for (std::size_t i = begin; i < end; ++i) {
        double si = 0.0;
        for (std::size_t k = 0; k < p; ++k) si += X(i, k) * X(i, k);
        double base_ii = si + kernel_c;
        double v_ii = 1.0;
        for (int e = 0; e < deg; ++e) v_ii *= base_ii;
        K(i, i) = v_ii;
        for (std::size_t j = i + 1; j < n; ++j) {
          const double base = krls_dot_ij(X, i, j, p) + kernel_c;
          double v = 1.0;
          for (int e = 0; e < deg; ++e) v *= base;
          K(i, j) = v;
          K(j, i) = v;
        }
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix krls_kernel_cpp(const NumericMatrix& X,
                              const NumericVector& sigma_vec,
                              int kernel_type = 0,
                              double kernel_c = 1.0) {
  const int n = X.nrow();
  const int p = X.ncol();
  if (kernel_type == 0 && sigma_vec.size() != p) {
    Rcpp::stop("krls_kernel_cpp: length(sigma_vec) (%d) must equal ncol(X) (%d)",
               sigma_vec.size(), p);
  }
  if (kernel_type < 0 || kernel_type > 5) {
    Rcpp::stop("krls_kernel_cpp: kernel_type must be in 0..5 (got %d)",
               kernel_type);
  }
  // Detect bit-exact constant sigma_vec (broadcast from scalar).
  // Only meaningful for Gaussian; harmless to compute for others.
  bool sigma_constant = true;
  double s0 = (sigma_vec.size() > 0) ? sigma_vec[0] : 1.0;
  for (int k = 1; k < sigma_vec.size(); ++k) {
    if (sigma_vec[k] != s0) { sigma_constant = false; break; }
  }
  std::vector<double> inv_sig(sigma_vec.size());
  for (int k = 0; k < sigma_vec.size(); ++k)
    inv_sig[k] = 1.0 / sigma_vec[k];
  const double inv_sigma_scalar = (sigma_vec.size() > 0) ? (1.0 / s0) : 0.0;

  NumericMatrix K(n, n);
  KrlsKernelWorker w(X, K, inv_sig, sigma_constant, inv_sigma_scalar,
                     kernel_type, kernel_c);
  parallelFor(0, n, w);
  return K;
}

// -------------------------------------------------------------------
// Test-vs-train Gaussian kernel (m_new x n_train).
//
// Same FP-determinism back-compat strategy as the train-train worker:
// constant-vector fast path mirrors v0.0.0.9044 byte-for-byte.
// -------------------------------------------------------------------
struct KrlsKernelPredWorker : public Worker {
  const RMatrix<double> Xnew;
  const RMatrix<double> Xtr;
  RMatrix<double> K;
  const std::vector<double> inv_sig;
  const bool sigma_constant;
  const double inv_sigma_scalar;
  const int kernel_type;
  const double kernel_c;

  KrlsKernelPredWorker(const NumericMatrix& Xnew_, const NumericMatrix& Xtr_,
                       NumericMatrix& K_,
                       const std::vector<double>& inv_sig_,
                       bool sigma_constant_, double inv_sigma_scalar_,
                       int kernel_type_, double kernel_c_)
    : Xnew(Xnew_), Xtr(Xtr_), K(K_), inv_sig(inv_sig_),
      sigma_constant(sigma_constant_),
      inv_sigma_scalar(inv_sigma_scalar_),
      kernel_type(kernel_type_), kernel_c(kernel_c_) {}

  void operator()(std::size_t begin, std::size_t end) {
    const std::size_t n = Xtr.nrow();
    const std::size_t p = Xtr.ncol();
    if (kernel_type == 0) {
      if (sigma_constant) {
        const double inv_sigma = inv_sigma_scalar;
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
      } else {
        for (std::size_t i = begin; i < end; ++i) {
          for (std::size_t j = 0; j < n; ++j) {
            double s = 0.0;
            for (std::size_t k = 0; k < p; ++k) {
              const double d = Xnew(i, k) - Xtr(j, k);
              s += d * d * inv_sig[k];
            }
            K(i, j) = std::exp(-s);
          }
        }
      }
    } else if (kernel_type == 1) {
      for (std::size_t i = begin; i < end; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
          K(i, j) = krls_dot_ij_xy(Xnew, Xtr, i, j, p);
        }
      }
    } else {
      const int deg = kernel_type - 1;
      for (std::size_t i = begin; i < end; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
          const double base = krls_dot_ij_xy(Xnew, Xtr, i, j, p) + kernel_c;
          double v = 1.0;
          for (int e = 0; e < deg; ++e) v *= base;
          K(i, j) = v;
        }
      }
    }
  }
};

// [[Rcpp::export]]
NumericMatrix krls_kernel_pred_cpp(const NumericMatrix& Xnew,
                                   const NumericMatrix& Xtrain,
                                   const NumericVector& sigma_vec,
                                   int kernel_type = 0,
                                   double kernel_c = 1.0) {
  const int m = Xnew.nrow();
  const int n = Xtrain.nrow();
  const int p = Xtrain.ncol();
  if (kernel_type == 0 && sigma_vec.size() != p) {
    Rcpp::stop("krls_kernel_pred_cpp: length(sigma_vec) (%d) must equal ncol(Xtrain) (%d)",
               sigma_vec.size(), p);
  }
  if (kernel_type < 0 || kernel_type > 5) {
    Rcpp::stop("krls_kernel_pred_cpp: kernel_type must be in 0..5 (got %d)",
               kernel_type);
  }
  bool sigma_constant = true;
  double s0 = (sigma_vec.size() > 0) ? sigma_vec[0] : 1.0;
  for (int k = 1; k < sigma_vec.size(); ++k) {
    if (sigma_vec[k] != s0) { sigma_constant = false; break; }
  }
  std::vector<double> inv_sig(sigma_vec.size());
  for (int k = 0; k < sigma_vec.size(); ++k)
    inv_sig[k] = 1.0 / sigma_vec[k];
  const double inv_sigma_scalar = (sigma_vec.size() > 0) ? (1.0 / s0) : 0.0;

  NumericMatrix K(m, n);
  KrlsKernelPredWorker w(Xnew, Xtrain, K, inv_sig, sigma_constant,
                         inv_sigma_scalar, kernel_type, kernel_c);
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

// GCV loss (Craven & Wahba 1979) on the SAME eigenbasis used by LOO.
//
//   trH(lam)   = sum d / (d + lam)
//   num_kept   = sum( (lam / (d + lam))^2 * Vty^2 )
//   resid_drop = yty - sum(Vty^2)
//   GCV        = (num_kept + resid_drop) / n /
//                pow( max(1.0 - trH/n, DEN_FLOOR), 2 )
//
// DEN_FLOOR clamps the denominator away from 0 in the interpolation
// limit (lam -> 0). Picked as 1e-8.
// [[Rcpp::export]]
double krls_gcv_loss_cpp(const arma::vec& d, const arma::vec& Vty,
                         double yty, int n_y, double lambda) {
  const arma::vec invshrink = lambda / (d + lambda);
  const double trH   = arma::sum(d / (d + lambda));
  const double num_k = arma::dot(arma::square(invshrink), arma::square(Vty));
  const double resid_drop = std::max(0.0, yty - arma::dot(Vty, Vty));
  const double n     = static_cast<double>(n_y);
  const double denom = std::max(1.0 - trH / n, 1.0e-8);
  return (num_k + resid_drop) / n / (denom * denom);
}

// -------------------------------------------------------------------
// Pointwise marginal effects.
//
//   derivmat[i, k] = (-2 / sigma_k) * sum_j (X_ik - X_jk) * K_ij * c_j
//                  = (-2 / sigma_k) * ( X_ik * (K c)_i - (K diag(c) X)_ik )
//
// Phase 2a (v0.0.0.9045): per-feature scale via `sigma_vec`. When
// `sigma_vec` is bit-exact constant (broadcast from scalar), the
// previous scalar BLAS-1 multiply `(-2.0 / sigma) * (...)` is used so
// the scalar path is byte-identical to v0.0.0.9044.
//
// All operations are matrix-multiplies dispatched to BLAS; no explicit
// per-variable n x n L matrix is constructed.
// -------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat krls_deriv_cpp(const arma::mat& X, const arma::mat& K,
                         const arma::vec& c, const arma::vec& sigma_vec,
                         int kernel_type = 0, double kernel_c = 1.0) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  if (kernel_type == 0) {
    if (sigma_vec.n_elem != p) {
      Rcpp::stop("krls_deriv_cpp: length(sigma_vec) must equal ncol(X)");
    }
    const arma::vec s1 = K * c;                  // n x 1, equals yhat_std.
    arma::mat Xc       = X;                      // n x p, will be scaled.
    Xc.each_col() %= c;                          // row j scaled by c_j.
    const arma::mat s2 = K * Xc;                 // n x p.
    arma::mat Xs1      = X;
    Xs1.each_col() %= s1;                        // row i scaled by s1_i.
    arma::mat D = Xs1 - s2;
    // Constant-vector fast path: preserves FP byte-identity vs v9044.
    bool sigma_constant = true;
    const double s0 = sigma_vec(0);
    for (arma::uword k = 1; k < p; ++k) {
      if (sigma_vec(k) != s0) { sigma_constant = false; break; }
    }
    if (sigma_constant) {
      return (-2.0 / s0) * D;
    }
    arma::rowvec scl = (-2.0 / sigma_vec).t();
    D.each_row() %= scl;
    return D;
  }
  if (kernel_type == 1) {
    // Linear: D[i, k] = sum_j x_{j,k} c_j = (X' c)_k. Constant per col.
    arma::rowvec colvec = (X.t() * c).t();        // 1 x p
    arma::mat D(n, p);
    for (arma::uword k = 0; k < p; ++k) D.col(k).fill(colvec(k));
    return D;
  }
  // Polynomial degree d = kernel_type - 1.
  // K_ij = (g_ij + c)^d  =>  H_ij = (g_ij + c)^(d-1)  with g_ij = x_i . x_j
  // dK_ij/dx_ik = d * H_ij * x_{j,k}
  // D[i, k] = sum_j d * H_ij * x_{j,k} * c_j = d * (H @ diag(c) @ X)[i, k].
  const int deg = kernel_type - 1;
  if (deg <= 0) {
    Rcpp::stop("krls_deriv_cpp: polynomial degree must be >= 1");
  }
  // Build H = (G + c)^(d-1) where G = X X'.
  arma::mat H;
  if (deg == 1) {
    // (d-1) = 0  =>  H is all ones; D[i, k] = sum_j x_{j,k} c_j (== linear case).
    arma::rowvec colvec = (X.t() * c).t();
    arma::mat D(n, p);
    for (arma::uword k = 0; k < p; ++k) D.col(k).fill(colvec(k));
    return D;
  }
  arma::mat G = X * X.t();
  arma::mat base = G + kernel_c;
  H = base;
  for (int e = 1; e < deg - 1; ++e) {
    H = H % base;   // elementwise power
  }
  // Two dgemm-style products: Xc (row j scaled by c_j), then H @ Xc.
  arma::mat Xc = X;
  Xc.each_col() %= c;
  arma::mat D = H * Xc;
  return static_cast<double>(deg) * D;
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
                                 const arma::vec& sigma_vec, double lambda,
                                 double sigma2,
                                 int kernel_type = 0,
                                 double kernel_c = 1.0) {
  const arma::uword n = X.n_rows;
  const arma::uword p = X.n_cols;
  if (kernel_type != 0) {
    // Phase Q2: derivation of avg-deriv variance under non-Gaussian
    // kernels is deferred to Q6. Return all-NA so the R-side caller can
    // emit a one-shot warning + leave var.avgderivatives = NA.
    arma::vec out(p);
    out.fill(NA_REAL);
    return out;
  }
  if (sigma_vec.n_elem != p) {
    Rcpp::stop("krls_avg_deriv_var_cpp: length(sigma_vec) must equal ncol(X)");
  }
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
  // Constant-vector fast path: preserves FP byte-identity vs v9044.
  bool sigma_constant = true;
  const double s0 = sigma_vec(0);
  for (arma::uword k = 1; k < p; ++k) {
    if (sigma_vec(k) != s0) { sigma_constant = false; break; }
  }
  if (sigma_constant) {
    const double scale = 4.0 / (s0 * s0) /
                         static_cast<double>(n * n);
    return scale * out;
  }
  arma::vec scale_k = 4.0 / arma::square(sigma_vec) /
                      static_cast<double>(n) /
                      static_cast<double>(n);
  return scale_k % out;
}

// -------------------------------------------------------------------
// Convenience: V .* V (used to compute diag(G^{-1}) at each lambda).
// Exposed so R-side can precompute once outside the lambda search loop.
// -------------------------------------------------------------------
// [[Rcpp::export]]
arma::mat krls_vsq_cpp(const arma::mat& V) {
  return V % V;
}

// -------------------------------------------------------------------
// Phase Q2 (v0.0.0.9049): GP posterior variance for KRLS predictions.
//
// For a fit (V, d, lambda) on training kernel K, the posterior variance
// at a new point x* is
//   var(f(x*)) = K(x*, x*) - K*' (K + lambda I)^{-1} K*
// where K* = K(X_tr, x*).  In the eigen-basis of K:
//   (K + lambda I)^{-1} = V diag(1/(d + lambda)) V'.
// Returns a length-n_new vector of variances on the *standardised* y
// scale.  Floors at 0 to absorb FP cancellation.
//
// kernel_type / kernel_c determine the kernel.  For Gaussian (0) the
// diagonal K(x*, x*) = 1; for linear, K(x*, x*) = ||x*||^2; for poly,
// K(x*, x*) = (||x*||^2 + c)^d.
// -------------------------------------------------------------------
// [[Rcpp::export]]
arma::vec krls_posterior_var_cpp(const arma::mat& Xnew,
                                 const arma::mat& Xtr,
                                 const arma::mat& V,
                                 const arma::vec& d,
                                 double lambda,
                                 int kernel_type,
                                 const arma::vec& sigma_vec,
                                 double kernel_c) {
  if (kernel_type < 0 || kernel_type > 5) {
    Rcpp::stop("krls_posterior_var_cpp: kernel_type must be 0..5");
  }
  const arma::uword m = Xnew.n_rows;
  const arma::uword n = Xtr.n_rows;
  if (V.n_rows != n) {
    Rcpp::stop("krls_posterior_var_cpp: V must have nrow(Xtr) rows");
  }
  // Build K* = K(Xnew, Xtr) (m x n).
  arma::mat Kstar(m, n);
  if (kernel_type == 0) {
    // Gaussian: K_ij = exp(- sum_k (xnew_ik - xtr_jk)^2 / sigma_k).
    bool sigma_constant = true;
    const double s0 = sigma_vec(0);
    for (arma::uword k = 1; k < sigma_vec.n_elem; ++k) {
      if (sigma_vec(k) != s0) { sigma_constant = false; break; }
    }
    if (sigma_constant) {
      const double inv_sigma = 1.0 / s0;
      for (arma::uword i = 0; i < m; ++i) {
        for (arma::uword j = 0; j < n; ++j) {
          double s = 0.0;
          for (arma::uword k = 0; k < Xnew.n_cols; ++k) {
            const double dxn = Xnew(i, k) - Xtr(j, k);
            s += dxn * dxn;
          }
          Kstar(i, j) = std::exp(-s * inv_sigma);
        }
      }
    } else {
      for (arma::uword i = 0; i < m; ++i) {
        for (arma::uword j = 0; j < n; ++j) {
          double s = 0.0;
          for (arma::uword k = 0; k < Xnew.n_cols; ++k) {
            const double dxn = Xnew(i, k) - Xtr(j, k);
            s += dxn * dxn / sigma_vec(k);
          }
          Kstar(i, j) = std::exp(-s);
        }
      }
    }
  } else if (kernel_type == 1) {
    Kstar = Xnew * Xtr.t();
  } else {
    const int deg = kernel_type - 1;
    arma::mat base = Xnew * Xtr.t() + kernel_c;
    Kstar = base;
    for (int e = 1; e < deg; ++e) Kstar = Kstar % base;
  }

  // Compute diag(K* V diag(1/(d+lam)) V' K*'):
  //   M = K* V              (m x n)
  //   row_i sum_j M_ij^2 * 1/(d_j + lambda)
  arma::mat M = Kstar * V;                          // m x n
  arma::vec inv = 1.0 / (d + lambda);
  arma::mat M2 = M % M;
  arma::vec quad = M2 * inv;                        // m x 1

  // Build diagonal of K(x*, x*).
  arma::vec kss(m);
  if (kernel_type == 0) {
    kss.fill(1.0);
  } else if (kernel_type == 1) {
    for (arma::uword i = 0; i < m; ++i) {
      double s = 0.0;
      for (arma::uword k = 0; k < Xnew.n_cols; ++k) {
        s += Xnew(i, k) * Xnew(i, k);
      }
      kss(i) = s;
    }
  } else {
    const int deg = kernel_type - 1;
    for (arma::uword i = 0; i < m; ++i) {
      double s = 0.0;
      for (arma::uword k = 0; k < Xnew.n_cols; ++k) {
        s += Xnew(i, k) * Xnew(i, k);
      }
      double v = 1.0;
      double base = s + kernel_c;
      for (int e = 0; e < deg; ++e) v *= base;
      kss(i) = v;
    }
  }
  arma::vec var = kss - quad;
  // FP-cancellation floor at 0.
  var.elem(arma::find(var < 0.0)).zeros();
  return var;
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

// Worker: each task is one sigma cell. Writes into pre-sized output
// buffers indexed by the sigma's slot in sigma_grid. No cross-thread
// shared mutable state inside krls_one_sigma -- all temporaries are
// stack-local arma::mat / arma::vec.
struct KrlsSigmaWorker : public RcppParallel::Worker {
  const arma::mat& D_tr;
  const arma::mat& D_te;
  const arma::vec& y_tr;
  const arma::vec& y_te;
  const arma::vec& sigma_grid;
  double lam_tol, L0, L_step;
  arma::vec& mse;
  arma::vec& lam;

  KrlsSigmaWorker(const arma::mat& D_tr_, const arma::mat& D_te_,
                  const arma::vec& y_tr_, const arma::vec& y_te_,
                  const arma::vec& sigma_grid_,
                  double lam_tol_, double L0_, double L_step_,
                  arma::vec& mse_, arma::vec& lam_)
    : D_tr(D_tr_), D_te(D_te_), y_tr(y_tr_), y_te(y_te_),
      sigma_grid(sigma_grid_), lam_tol(lam_tol_), L0(L0_), L_step(L_step_),
      mse(mse_), lam(lam_) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      double m = 0.0, l = 0.0;
      krls_one_sigma(D_tr, D_te, y_tr, y_te, sigma_grid(i),
                     lam_tol, L0, L_step, m, l);
      mse(i) = m;
      lam(i) = l;
    }
  }
};

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

  // Clamp worker count: at least 1, no more than nsigma (no idle workers).
  int n_workers = nthreads;
  if (n_workers < 1) n_workers = 1;
  if ((arma::uword) n_workers > nsigma) n_workers = (int) nsigma;

  KrlsSigmaWorker worker(D_tr, D_te, y_tr, y_te, sigma_grid,
                         lam_tol, L0, L_step, mse, lam);
  // Constrain the parallel_for to n_workers by passing numThreads directly
  // (mirrors the Phase 2 Nystrom inner). Avoids tbb::task_arena, which
  // depends on TBB symbols not shipped by Rtools45 on Windows.
  RcppParallel::parallelFor(0, nsigma, worker, /*grainSize=*/1,
                            /*numThreads=*/n_workers);

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
    Rcpp::Named("nthreads_used")    = n_workers
  );
}

// ---------------------------------------------------------------------------
// Phase 2 (v0.0.0.9043): Nystrom single-fit at fixed sigma.
//
// Builds C = K(X, Z), W = K(Z, Z), eig_sym(W) -> (U, D), regularizes via
// nystrom_eps relative ridge, builds Phi = C @ U @ diag(D_reg^{-1/2}),
// SVD of Phi for the LOO lambda objective, golden-section search,
// recovers alpha = U @ (Dinvsqrt * beta).
//
// Reference: /tmp/krls-reference/R/nystrom.R::.fit_krls_nystrom
//
// Bracket scheme: LINEAR `L_lam += 0.05` step matches the reference R
// `.nystrom_lambda_bounds`. Required so the parity test (NYS-FIT-1) can
// hit `tol = 1e-7` on coeffs against the reference implementation; using a
// log-scale step would bracket lambda from a different starting point and
// drift the golden-section minimum by O(1e-5).
// ---------------------------------------------------------------------------

static double nystrom_loo_loss(const arma::mat& U_phi,
                               const arma::vec& Sigma2,
                               const arma::vec& y,
                               double lambda) {
  arma::vec w     = Sigma2 / (Sigma2 + lambda);
  arma::vec Uty   = U_phi.t() * y;
  arma::vec yfit  = U_phi * (w % Uty);
  // diagS_i = sum_k U_phi(i,k)^2 * w_k
  arma::mat Usq = arma::square(U_phi);
  arma::vec diagS = Usq * w;
  arma::vec denom = 1.0 - diagS;
  if (arma::any(denom <= std::numeric_limits<double>::epsilon())) {
    return std::numeric_limits<double>::max();
  }
  arma::vec resid = (y - yfit) / denom;
  return arma::dot(resid, resid);
}

// [[Rcpp::export]]
Rcpp::List krls_nystrom_fit_cpp(const arma::mat& X_tr, const arma::mat& Z,
                                const arma::vec& y_tr, double sigma,
                                Rcpp::List lambda_args, double eps,
                                bool compute_vcov) {
  const arma::uword n = X_tr.n_rows;
  const arma::uword m = Z.n_rows;

  // C = K(X_tr, Z), W = K(Z, Z)
  arma::mat D_tr_Z = krls_pairwise_sqdist_cpp(X_tr, Z);
  arma::mat D_Z_Z  = krls_pairwise_sqdist_cpp(Z,    Z);
  arma::mat C = arma::exp(-D_tr_Z / sigma);
  arma::mat W = arma::exp(-D_Z_Z  / sigma);

  // eig_sym(W) -> (U, D); arma returns ASCENDING. Reference uses base::eigen
  // (descending) but operations downstream are invariant to eigenvector
  // ordering (we apply Dinvsqrt componentwise and combine via U * (...)),
  // so leaving ASCENDING is numerically equivalent.
  arma::vec D;
  arma::mat U;
  arma::eig_sym(D, U, W);
  D = arma::clamp(D, 0.0, arma::datum::inf);
  double Dmax = D.max();
  if (Dmax <= 0.0) {
    Rcpp::stop("krls_nystrom_fit_cpp: anchor W numerically zero; try larger sigma");
  }
  arma::vec D_reg   = arma::clamp(D, eps * Dmax, arma::datum::inf);
  arma::vec Dinv2   = 1.0 / arma::sqrt(D_reg);
  arma::uword floored_count = arma::sum(arma::conv_to<arma::uvec>::from(D < eps * Dmax));
  double D_min_raw = D.min();

  // Phi = C @ U @ diag(Dinv2)   (n x m)
  arma::mat Phi = C * U;
  Phi.each_row() %= Dinv2.t();

  // SVD of Phi (economy form to match base::svd() in the reference: for
  // Phi n x m with n >= m we want U n x m, V m x m, d length m).
  arma::mat U_phi, V_phi;
  arma::vec sigma_phi;
  bool ok = arma::svd_econ(U_phi, sigma_phi, V_phi, Phi);
  if (!ok) {
    Rcpp::stop("krls_nystrom_fit_cpp: SVD of Phi failed");
  }
  arma::vec Sigma2 = arma::square(sigma_phi);
  if (Sigma2.max() <= 0.0) {
    Rcpp::stop("krls_nystrom_fit_cpp: Nystrom feature spectrum is zero; "
               "increase sigma or change landmarks");
  }

  // Golden-section lambda search bounds (anchored at Sigma2)
  double L0     = Rcpp::as<double>(lambda_args["L0"]);
  double tol    = Rcpp::as<double>(lambda_args["tol"]);

  double Smax = Sigma2.max();
  double U_lam = Smax;
  int iter = 0;
  while (arma::sum(Sigma2 / (Sigma2 + U_lam)) >= 1.0 && iter < 200) {
    U_lam *= 2.0;
    ++iter;
  }
  // L: index q = which.min(|Sigma2 - Smax/1000|), 1-based equivalent.
  // Matches the reference (which.min returns the 1-based position of the
  // smallest |Sigma2 - target|; here we use 0-based q_idx then add 1).
  double tgt = Smax / 1000.0;
  arma::uword q_idx = 0;
  double mind = std::abs(Sigma2(0) - tgt);
  for (arma::uword i = 1; i < Sigma2.n_elem; ++i) {
    double dd = std::abs(Sigma2(i) - tgt);
    if (dd < mind) { mind = dd; q_idx = i; }
  }
  double q_R = (double)(q_idx + 1);
  double L_lam = L0;
  // LINEAR step matching reference R `.nystrom_lambda_bounds` (L + 0.05).
  // Per the user/spec decision (Task 3 critical numerical detail), the
  // C++ uses the same linear climb so NYS-FIT-1 reaches tol=1e-7 on coeffs.
  while (arma::sum(Sigma2 / (Sigma2 + L_lam)) > q_R && L_lam < Smax) {
    L_lam += 0.05;
  }
  if (!(L_lam < U_lam)) {
    Rcpp::stop("krls_nystrom_fit_cpp: lambda bounds collapsed; "
               "supply L and U explicitly");
  }

  // Golden-section search (mirrors reference: X1 = L + gr*(U-L),
  // X2 = U - gr*(U-L), with gr = 0.381966 = (3 - sqrt(5)) / 2).
  const double gr = (3.0 - std::sqrt(5.0)) / 2.0;
  double X1 = L_lam + gr * (U_lam - L_lam);
  double X2 = U_lam - gr * (U_lam - L_lam);
  double S1 = nystrom_loo_loss(U_phi, Sigma2, y_tr, X1);
  double S2 = nystrom_loo_loss(U_phi, Sigma2, y_tr, X2);
  int it = 0;
  while (std::abs(S1 - S2) > tol && it < 1000) {
    if (S1 < S2) {
      U_lam = X2; X2 = X1; X1 = L_lam + gr * (U_lam - L_lam);
      S2 = S1; S1 = nystrom_loo_loss(U_phi, Sigma2, y_tr, X1);
    } else {
      L_lam = X1; X1 = X2; X2 = U_lam - gr * (U_lam - L_lam);
      S1 = S2; S2 = nystrom_loo_loss(U_phi, Sigma2, y_tr, X2);
    }
    ++it;
  }
  // Reference returns whichever of X1/X2 has the smaller objective; not
  // the midpoint. Matching that to keep the parity test tight.
  double lambda = (S1 < S2) ? X1 : X2;

  // Solve beta + alpha
  arma::vec Uty   = U_phi.t() * y_tr;
  arma::vec w     = Sigma2 / (Sigma2 + lambda);
  arma::vec yfit  = U_phi * (w % Uty);
  arma::vec beta_coef = V_phi * ((sigma_phi % Uty) / (Sigma2 + lambda));
  arma::vec alpha = U * (Dinv2 % beta_coef);

  // Optional vcov
  arma::mat vcov_alpha;
  if (compute_vcov) {
    double sigmasq = arma::dot(y_tr - yfit, y_tr - yfit) / (double) n;
    arma::vec beta_weights = sigmasq * Sigma2 / arma::square(Sigma2 + lambda);
    arma::mat Vw = V_phi;
    Vw.each_row() %= beta_weights.t();
    arma::mat vcov_beta = Vw * V_phi.t();
    arma::mat alpha_map = U;
    alpha_map.each_row() %= Dinv2.t();
    vcov_alpha = alpha_map * vcov_beta * alpha_map.t();
  }

  Rcpp::List W_eigen = Rcpp::List::create(
    Rcpp::Named("values")  = D,
    Rcpp::Named("vectors") = U
  );

  // coeffs as length-m column matrix (matches reference shape)
  arma::mat coeffs_mat(alpha.n_elem, 1);
  coeffs_mat.col(0) = alpha;

  return Rcpp::List::create(
    Rcpp::Named("coeffs")           = coeffs_mat,
    Rcpp::Named("fitted_std")       = yfit,
    Rcpp::Named("lambda")           = lambda,
    Rcpp::Named("landmarks")        = Z,
    Rcpp::Named("W_eigen")          = W_eigen,
    Rcpp::Named("Dinvsqrt")         = Dinv2,
    Rcpp::Named("Sigma2")           = Sigma2,
    Rcpp::Named("vcov_alpha")       = vcov_alpha,
    Rcpp::Named("floored_count")    = (int) floored_count,
    Rcpp::Named("D_min_raw")        = D_min_raw,
    Rcpp::Named("D_max_raw")        = Dmax,
    Rcpp::Named("nystrom_m")        = (int) m
  );
}

// ---------------------------------------------------------------------------
// Phase 2 (v0.0.0.9043): Nystrom predict.
//
//   yhat = exp(-||X_new - Z||^2 / sigma) @ alpha
//
// Uses the dgemm identity in krls_pairwise_sqdist_cpp for the n_new x m
// cross-distance matrix, then elementwise exp + matrix-vector product.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
Rcpp::NumericVector krls_nystrom_predict_cpp(const arma::mat& X_new,
                                             const arma::mat& Z,
                                             const arma::vec& alpha,
                                             double sigma) {
  if (X_new.n_cols != Z.n_cols) {
    Rcpp::stop("krls_nystrom_predict_cpp: X_new and Z must have same n_cols");
  }
  if (alpha.n_elem != Z.n_rows) {
    Rcpp::stop("krls_nystrom_predict_cpp: length(alpha) must equal nrow(Z)");
  }
  arma::mat D = krls_pairwise_sqdist_cpp(X_new, Z);   // n_new x m
  arma::mat K = arma::exp(-D / sigma);
  arma::vec yhat = K * alpha;
  // Return as a plain numeric vector (no dim attribute) so callers see a
  // length-n_new vector consistent with predict.krls_rr's existing shape.
  Rcpp::NumericVector out(yhat.n_elem);
  std::copy(yhat.begin(), yhat.end(), out.begin());
  return out;
}

// ---------------------------------------------------------------------------
// Phase 2 (v0.0.0.9043): Nystrom autotune inner — parallel sweep over sigma.
//
// For a single CV fold, evaluate the full sigma grid in parallel. Each task
// builds C/W/K_te on the fold's (X_tr, Z, X_te), runs the Nystrom-LOO
// golden-section lambda search, computes alpha, and scores the held-out MSE.
// Reuses nystrom_loo_loss() from the single-fit path so the LOO objective is
// byte-identical to krls_nystrom_fit_cpp.
//
// Bracket scheme: LINEAR `L_lam += 0.05` step — matches Task 3 so the
// per-fold lambda choice (and downstream alpha + MSE) lines up with what the
// non-autotune path would emit at the same fold partition.
//
// Threading: parallelFor receives numThreads directly (no tbb::task_arena —
// Rtools45 on Windows does not ship the task_arena_base symbols). Worker
// count clamped to [1, nsigma].
// ---------------------------------------------------------------------------

struct KrlsNystromSigmaWorker : public RcppParallel::Worker {
  const arma::mat& D_tr_Z;
  const arma::mat& D_Z_Z;
  const arma::mat& D_te_Z;
  const arma::vec& y_tr;
  const arma::vec& y_te;
  const arma::vec& sigma_grid;
  double L0, tol, eps;
  arma::vec& mse;
  arma::vec& lam;

  KrlsNystromSigmaWorker(const arma::mat& D_tr_Z_, const arma::mat& D_Z_Z_,
                         const arma::mat& D_te_Z_,
                         const arma::vec& y_tr_, const arma::vec& y_te_,
                         const arma::vec& sigma_grid_,
                         double L0_, double tol_, double eps_,
                         arma::vec& mse_, arma::vec& lam_)
    : D_tr_Z(D_tr_Z_), D_Z_Z(D_Z_Z_), D_te_Z(D_te_Z_),
      y_tr(y_tr_), y_te(y_te_), sigma_grid(sigma_grid_),
      L0(L0_), tol(tol_), eps(eps_),
      mse(mse_), lam(lam_) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      double s = sigma_grid(i);
      arma::mat C    = arma::exp(-D_tr_Z / s);
      arma::mat W    = arma::exp(-D_Z_Z  / s);
      arma::mat K_te = arma::exp(-D_te_Z / s);

      arma::vec D; arma::mat U;
      bool eig_ok = arma::eig_sym(D, U, W);
      if (!eig_ok) {
        mse(i) = std::numeric_limits<double>::max(); lam(i) = 0.0; continue;
      }
      D = arma::clamp(D, 0.0, arma::datum::inf);
      double Dmax = D.max();
      if (Dmax <= 0.0) {
        mse(i) = std::numeric_limits<double>::max(); lam(i) = 0.0; continue;
      }
      arma::vec D_reg = arma::clamp(D, eps * Dmax, arma::datum::inf);
      arma::vec Dinv2 = 1.0 / arma::sqrt(D_reg);

      arma::mat Phi = C * U;
      Phi.each_row() %= Dinv2.t();

      arma::mat U_phi, V_phi; arma::vec sigma_phi;
      bool ok = arma::svd_econ(U_phi, sigma_phi, V_phi, Phi);
      if (!ok) {
        mse(i) = std::numeric_limits<double>::max(); lam(i) = 0.0; continue;
      }
      arma::vec Sigma2 = arma::square(sigma_phi);
      double Smax = Sigma2.max();
      if (Smax <= 0.0) {
        mse(i) = std::numeric_limits<double>::max(); lam(i) = 0.0; continue;
      }

      // Lambda bracket — mirror krls_nystrom_fit_cpp (linear L += 0.05).
      double U_lam = Smax;
      int it = 0;
      while (arma::sum(Sigma2 / (Sigma2 + U_lam)) >= 1.0 && it < 200) {
        U_lam *= 2.0; ++it;
      }
      double tgt = Smax / 1000.0;
      arma::uword q_idx = 0;
      double mind = std::abs(Sigma2(0) - tgt);
      for (arma::uword j = 1; j < Sigma2.n_elem; ++j) {
        double dd = std::abs(Sigma2(j) - tgt);
        if (dd < mind) { mind = dd; q_idx = j; }
      }
      double q_R = (double)(q_idx + 1);
      double L_lam = L0;
      while (arma::sum(Sigma2 / (Sigma2 + L_lam)) > q_R && L_lam < Smax) {
        L_lam += 0.05;
      }
      if (!(L_lam < U_lam)) {
        mse(i) = std::numeric_limits<double>::max(); lam(i) = 0.0; continue;
      }

      // Golden-section search — same gr and X1/X2 layout as single-fit path.
      const double gr = (3.0 - std::sqrt(5.0)) / 2.0;
      double X1 = L_lam + gr * (U_lam - L_lam);
      double X2 = U_lam - gr * (U_lam - L_lam);
      double S1 = nystrom_loo_loss(U_phi, Sigma2, y_tr, X1);
      double S2 = nystrom_loo_loss(U_phi, Sigma2, y_tr, X2);
      int it2 = 0;
      while (std::abs(S1 - S2) > tol && it2 < 1000) {
        if (S1 < S2) {
          U_lam = X2; X2 = X1; X1 = L_lam + gr * (U_lam - L_lam);
          S2 = S1; S1 = nystrom_loo_loss(U_phi, Sigma2, y_tr, X1);
        } else {
          L_lam = X1; X1 = X2; X2 = U_lam - gr * (U_lam - L_lam);
          S1 = S2; S2 = nystrom_loo_loss(U_phi, Sigma2, y_tr, X2);
        }
        ++it2;
      }
      double lambda = (S1 < S2) ? X1 : X2;

      // alpha + held-out MSE
      arma::vec Uty       = U_phi.t() * y_tr;
      arma::vec beta_coef = V_phi * ((sigma_phi % Uty) / (Sigma2 + lambda));
      arma::vec alpha     = U * (Dinv2 % beta_coef);
      arma::vec yhat_te   = K_te * alpha;

      mse(i) = arma::mean(arma::square(y_te - yhat_te));
      lam(i) = lambda;
    }
  }
};

// [[Rcpp::export]]
Rcpp::List krls_nystrom_autotune_inner_cpp(
    const arma::mat& X_tr, const arma::mat& Z, const arma::mat& X_te,
    const arma::vec& y_tr, const arma::vec& y_te,
    const arma::vec& sigma_grid, Rcpp::List lambda_args, double eps,
    int nthreads) {
  const arma::uword nsigma = sigma_grid.n_elem;
  arma::vec mse(nsigma, arma::fill::zeros);
  arma::vec lam(nsigma, arma::fill::zeros);

  arma::mat D_tr_Z = krls_pairwise_sqdist_cpp(X_tr, Z);
  arma::mat D_Z_Z  = krls_pairwise_sqdist_cpp(Z,    Z);
  arma::mat D_te_Z = krls_pairwise_sqdist_cpp(X_te, Z);

  double L0  = Rcpp::as<double>(lambda_args["L0"]);
  double tol = Rcpp::as<double>(lambda_args["tol"]);

  int n_workers = nthreads;
  if (n_workers < 1) n_workers = 1;
  if ((arma::uword) n_workers > nsigma) n_workers = (int) nsigma;

  KrlsNystromSigmaWorker worker(D_tr_Z, D_Z_Z, D_te_Z, y_tr, y_te, sigma_grid,
                                L0, tol, eps, mse, lam);
  RcppParallel::parallelFor(0, nsigma, worker, /*grainSize=*/1,
                            /*numThreads=*/n_workers);

  Rcpp::NumericVector mse_out(nsigma), lam_out(nsigma);
  for (arma::uword i = 0; i < nsigma; ++i) {
    mse_out[i] = mse(i);
    lam_out[i] = lam(i);
  }

  return Rcpp::List::create(
    Rcpp::Named("mse_per_sigma")    = mse_out,
    Rcpp::Named("lambda_per_sigma") = lam_out,
    Rcpp::Named("nthreads_used")    = n_workers
  );
}

// ---------------------------------------------------------------------------
// Phase Q5 (v0.0.0.9050): logistic-loss IRLS path.
//
// Optimises penalised binomial deviance for KRLS in coefficient space:
//   minimise    -2 * sum( y * log(p) + (1-y) * log(1-p) ) + lambda * c' K c
// where p = 1 / (1 + exp(-K c)) and lambda is the ridge penalty.
//
// Per-iteration Newton step (option ii in spec): solve
//   (K W K + lambda I) c = K W z
// with weights W = p*(1-p) and working response z = eta + (y - p) / W.
// Direct Cholesky-backed solve via Armadillo `solve_opts::likely_sympd`.
//
// Numerical safety:
//   - W floored at 1e-8 (prevents singular step under saturation).
//   - eta clipped to [-30, 30] (prevents `exp` overflow).
//   - Step-halving on deviance increase, up to 5 halvings per iter.
//   - Perfect-separation detection: warn + return last finite iterate
//     when ||c||_inf > 1e6 OR penalised dev < 1e-10.
//   - 50 outer iteration cap.
// ---------------------------------------------------------------------------

static inline arma::vec krls_clamp_eta(const arma::vec& eta, double lo, double hi) {
  arma::vec out = eta;
  out.elem(arma::find(out > hi)).fill(hi);
  out.elem(arma::find(out < lo)).fill(lo);
  return out;
}

static inline arma::vec krls_sigmoid_safe(const arma::vec& eta) {
  arma::vec e = krls_clamp_eta(eta, -30.0, 30.0);
  return 1.0 / (1.0 + arma::exp(-e));
}

// Penalised binomial deviance (the IRLS objective):
//   -2 * sum( w_i [ y_i log(p_i) + (1-y_i) log(1-p_i) ] ) + lambda * c' K c.
static double krls_pen_binomial_dev(const arma::vec& y, const arma::vec& p,
                                    const arma::vec& Kc, const arma::vec& c,
                                    const arma::vec& w_vec,
                                    double lambda, bool weighted) {
  const arma::uword n = y.n_elem;
  const double eps = 1e-15;
  double dev = 0.0;
  for (arma::uword i = 0; i < n; ++i) {
    const double pi = std::min(std::max(p(i), eps), 1.0 - eps);
    const double term = y(i) * std::log(pi) + (1.0 - y(i)) * std::log(1.0 - pi);
    const double wi = weighted ? w_vec(i) : 1.0;
    dev += wi * term;
  }
  const double pen = lambda * arma::dot(c, Kc);   // c' K c
  return -2.0 * dev + pen;
}

// [[Rcpp::export]]
Rcpp::List krls_irls_logistic_cpp(const arma::mat& K, const arma::vec& y,
                                  double lambda,
                                  Rcpp::Nullable<Rcpp::NumericVector> w_obs = R_NilValue,
                                  double tol = 1e-6, int max_iter = 50,
                                  int max_halve = 5, int trace = 0) {
  const arma::uword n = K.n_rows;
  if (K.n_cols != n) Rcpp::stop("krls_irls_logistic_cpp: K must be square");
  if (y.n_elem != n) Rcpp::stop("krls_irls_logistic_cpp: length(y) must equal nrow(K)");
  if (!(lambda > 0)) Rcpp::stop("krls_irls_logistic_cpp: lambda must be positive");

  const bool weighted = w_obs.isNotNull();
  arma::vec w_vec(n, arma::fill::ones);
  if (weighted) {
    Rcpp::NumericVector ww(w_obs);
    if ((arma::uword) ww.size() != n)
      Rcpp::stop("krls_irls_logistic_cpp: length(w_obs) must equal n");
    for (arma::uword i = 0; i < n; ++i) w_vec(i) = ww[i];
  }

  arma::vec c(n, arma::fill::zeros);
  arma::vec eta(n, arma::fill::zeros);
  arma::vec p = krls_sigmoid_safe(eta);
  arma::vec Kc = K * c;                             // == 0 initially
  double dev = krls_pen_binomial_dev(y, p, Kc, c, w_vec, lambda, weighted);

  arma::vec c_last_finite = c;
  arma::vec eta_last_finite = eta;
  arma::vec p_last_finite = p;
  double dev_last_finite = dev;
  int iter_last_finite = 0;

  const arma::mat I_n = arma::eye<arma::mat>(n, n);
  bool converged = false;
  bool separated = false;
  int iter_done = 0;

  for (int it = 1; it <= max_iter; ++it) {
    // W = p*(1-p) (with case weights baked in for the LHS), floored.
    arma::vec W = p % (1.0 - p);
    if (weighted) W = W % w_vec;
    W.transform([](double v) { return v < 1e-8 ? 1e-8 : v; });

    // z = eta + (y - p) / (p (1-p))  (unweighted denom; weights enter W).
    arma::vec denom = p % (1.0 - p);
    denom.transform([](double v) { return v < 1e-8 ? 1e-8 : v; });
    arma::vec z = eta + (y - p) / denom;

    // Solve (K W K + lambda I) c_new = K W z
    arma::mat KW = K.each_row() % W.t();           // K * diag(W); n x n
    arma::mat A  = KW * K + lambda * I_n;
    arma::vec b  = KW * z;

    arma::vec c_new;
    bool solved = arma::solve(c_new, A, b, arma::solve_opts::likely_sympd);
    if (!solved) {
      solved = arma::solve(c_new, A, b);
      if (!solved) { separated = true; break; }
    }

    // Step-halving on penalised deviance increase.
    arma::vec c_try = c_new;
    arma::vec eta_try, p_try, Kc_try;
    double dev_try = std::numeric_limits<double>::infinity();
    bool any_finite = false;
    bool accepted = false;
    arma::vec c_best = c;
    arma::vec eta_best = eta;
    arma::vec p_best = p;
    arma::vec Kc_best = Kc;
    double dev_best = dev;
    for (int h = 0; h <= max_halve; ++h) {
      Kc_try  = K * c_try;
      eta_try = krls_clamp_eta(Kc_try, -30.0, 30.0);
      p_try   = krls_sigmoid_safe(Kc_try);
      dev_try = krls_pen_binomial_dev(y, p_try, Kc_try, c_try, w_vec, lambda, weighted);
      if (std::isfinite(dev_try)) {
        any_finite = true;
        // Track the best (lowest dev) finite candidate seen so far.
        if (dev_try < dev_best) {
          c_best = c_try; eta_best = eta_try; p_best = p_try;
          Kc_best = Kc_try; dev_best = dev_try;
        }
        if (dev_try <= dev + 1e-12) { accepted = true; break; }
      }
      c_try = 0.5 * (c_try + c);
    }
    if (!any_finite) { separated = true; break; }
    // If no halve achieved descent, fall back to the best candidate seen.
    // If even that is no better than current dev, we've stalled — declare
    // convergence at the current iterate.
    if (!accepted) {
      if (dev_best >= dev) {
        converged = true;
        iter_done = it - 1;
        break;
      }
      c_try = c_best; eta_try = eta_best; p_try = p_best;
      Kc_try = Kc_best; dev_try = dev_best;
    }

    // Perfect-separation detection
    const double cn = arma::norm(c_try, "inf");
    if (cn > 1e6 || dev_try < 1e-10) {
      c   = c_try;
      eta = eta_try;
      p   = p_try;
      dev = dev_try;
      iter_done = it;
      separated = true;
      break;
    }

    // Convergence checks (pre-commit)
    const double c_norm_old = std::max(arma::norm(c, 2),
                                       std::numeric_limits<double>::epsilon());
    const double c_change   = arma::norm(c_try - c, 2) / c_norm_old;
    const double dev_change_rel = std::abs(dev_try - dev) /
                                    std::max(std::abs(dev), 1.0);

    // commit
    c   = c_try;
    eta = eta_try;
    p   = p_try;
    Kc  = Kc_try;
    dev = dev_try;

    // last-finite snapshot
    c_last_finite   = c;
    eta_last_finite = eta;
    p_last_finite   = p;
    dev_last_finite = dev;
    iter_last_finite = it;
    iter_done = it;

    if (trace > 0) {
      Rcpp::Rcout << "[IRLS] iter " << it
                  << "  dev=" << dev
                  << "  d_dev=" << dev_change_rel
                  << "  d_c="   << c_change << std::endl;
    }
    if (c_change < tol || dev_change_rel < tol) {
      converged = true;
      break;
    }
  }

  if (separated) {
    Rcpp::warning(
      "krls(loss='logistic'): possible perfect separation detected; "
      "returning last finite iterate. Try larger lambda or fewer features.");
    c   = c_last_finite;
    eta = eta_last_finite;
    p   = p_last_finite;
    dev = dev_last_finite;
    iter_done = iter_last_finite;
  }

  // Final W (for downstream LOO / vcov)
  arma::vec W_final = p % (1.0 - p);
  if (weighted) W_final = W_final % w_vec;
  W_final.transform([](double v) { return v < 1e-8 ? 1e-8 : v; });

  // H_diag = diag( K (K W K + lambda I)^{-1} K )
  arma::mat KW_f = K.each_row() % W_final.t();
  arma::mat A_f  = KW_f * K + lambda * I_n;
  arma::mat Ainv_K;
  bool ok_h = arma::solve(Ainv_K, A_f, K, arma::solve_opts::likely_sympd);
  if (!ok_h) ok_h = arma::solve(Ainv_K, A_f, K);
  arma::vec H_diag(n, arma::fill::zeros);
  if (ok_h) {
    arma::mat H = K * Ainv_K;
    H_diag = arma::diagvec(H);
  } else {
    H_diag.fill(NA_REAL);
  }

  return Rcpp::List::create(
    Rcpp::Named("coeffs")    = c,
    Rcpp::Named("eta")       = eta,
    Rcpp::Named("p")         = p,
    Rcpp::Named("W")         = W_final,
    Rcpp::Named("deviance")  = dev,
    Rcpp::Named("iter")      = iter_done,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("separated") = separated,
    Rcpp::Named("H_diag")    = H_diag
  );
}

// ---------------------------------------------------------------------------
// CT-2008 closed-form LOO deviance under logistic.
//
//   eta_loo[i] = eta[i] - (y[i] - p[i]) * H_diag[i] / (1 - W[i] * H_diag[i])
//   p_loo[i]   = plogis(eta_loo[i])
//   loo_dev    = -2 * sum( y * log(p_loo) + (1-y) * log(1-p_loo) )
// W = p*(1-p) (with floor); H_diag = diag( K (K W K + lambda I)^{-1} K ).
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double krls_logistic_loo_loss_cpp(const arma::vec& eta, const arma::vec& y,
                                  const arma::vec& p, const arma::vec& W,
                                  const arma::vec& H_diag) {
  const arma::uword n = eta.n_elem;
  if (y.n_elem != n || p.n_elem != n || W.n_elem != n || H_diag.n_elem != n)
    Rcpp::stop("krls_logistic_loo_loss_cpp: all inputs must have equal length");
  const double eps = 1e-15;
  double dev = 0.0;
  for (arma::uword i = 0; i < n; ++i) {
    double denom = 1.0 - W(i) * H_diag(i);
    if (denom < 1e-8) denom = 1e-8;
    const double eta_loo = eta(i) - (y(i) - p(i)) * H_diag(i) / denom;
    double e = eta_loo;
    if (e >  30.0) e =  30.0;
    if (e < -30.0) e = -30.0;
    const double p_loo = 1.0 / (1.0 + std::exp(-e));
    const double pi    = std::min(std::max(p_loo, eps), 1.0 - eps);
    dev += y(i) * std::log(pi) + (1.0 - y(i)) * std::log(1.0 - pi);
  }
  return -2.0 * dev;
}
