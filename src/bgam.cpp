// roadrunner -- src/bgam.cpp
//
// Component-wise P-spline gradient boosting engine for bgam().
//
// Exports (in call order during fit):
//   bgam_bspline_basis_cpp  -- B-spline design matrix (one predictor).
//   bgam_diff_penalty_cpp   -- difference penalty matrix P = D'D.
//   bgam_prefactor_cpp      -- Cholesky of A_j = B_j'B_j + lambda_j P_j
//                              (gaussian path only; cached once per fit).
//   bgam_boost_cpp          -- main boosting loop (TBB parallel_reduce over
//                              predictors, serial update).
//   bgam_predict_cpp        -- prediction on new data (link scale).
//
// Determinism invariant:
//   TBB parallel_reduce over j = 0..p-1 uses a BestResult accumulator that
//   always resolves ties by LOWEST j index. The join() method compares the
//   minimum rss of two accumulators and keeps the lower-j entry on a tie.
//   This produces byte-identical j* selection regardless of thread partition
//   or nthreads value. The nthreads argument is forwarded to
//   RcppParallel::parallelReduce via its numThreads parameter; nthreads=1
//   forces a serial scan with the same tie-breaking rule.
//
// No RNG in C++ code. All random operations (CV fold assignment, bagging)
// are performed in R with set.seed().

// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]

#define ARMA_NO_DEBUG

#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

// ============================================================================
// Utility: stable sigmoid
// ============================================================================

static inline double bgam_sigmoid(double x) {
  if (x >= 0.0)
    return 1.0 / (1.0 + std::exp(-x));
  double e = std::exp(x);
  return e / (1.0 + e);
}

// ============================================================================
// Cholesky solve helper: solve L L' x = b (lower triangular L).
// ============================================================================

static inline arma::vec bgam_chol_solve(const arma::mat& L,
                                        const arma::vec& b) {
  arma::vec y = arma::solve(arma::trimatl(L), b);
  return arma::solve(arma::trimatu(L.t()), y);
}

// ============================================================================
// bgam_bspline_basis_cpp
//
// Builds an n x K B-spline design matrix for a single predictor using the
// Cox-de Boor recursive algorithm with the given knot sequence.
//
// Arguments:
//   x      -- n training values for one predictor
//   knots  -- full knot sequence (including boundary repeats), length
//             nknots + 2*degree
//   degree -- B-spline degree (1..5)
//
// Returns a List with:
//   $B   -- n x K arma::mat (K = length(knots) - degree - 1)
//   $K   -- integer, number of B-spline columns
// ============================================================================

// [[Rcpp::export]]
Rcpp::List bgam_bspline_basis_cpp(
    const arma::vec& x,
    const arma::vec& knots,
    int degree
) {
  int n  = (int)x.n_elem;
  int nk = (int)knots.n_elem;
  int K  = nk - degree - 1;

  if (K < 1)
    Rcpp::stop("bgam_bspline_basis_cpp: not enough knots for given degree.");

  arma::mat B(n, K, arma::fill::zeros);

  double knot_lo = knots[0];
  double knot_hi = knots[nk - 1];

  // Index of the last non-degenerate lower span (contains [knot_lo, next unique]).
  // For a point at or below knot_lo, activate span 0 (first basis function = 1
  // after recursion for clamped splines).
  // For a point at or above knot_hi, activate span nk-2-degree (last non-
  // degenerate span, so the last basis function = 1 after recursion).
  // This gives boundary-constant extrapolation (same prediction as nearest
  // boundary point), providing non-zero SE outside the training range.
  int lo_span = 0;
  int hi_span = nk - 2 - degree;  // last non-degenerate span (0-indexed)
  // Validate hi_span; for K >= 1 this is always >= 0.
  if (hi_span < 0) hi_span = 0;

  for (int i = 0; i < n; ++i) {
    double xi = x[i];
    bool at_lo = (xi < knot_lo);
    bool at_hi = (xi > knot_hi);

    // Workspace: nk-1 B-spline values at degree d
    std::vector<double> bval(nk - 1, 0.0);

    if (at_lo) {
      // Extrapolation below lower boundary: activate the first span.
      // After degree recursion this gives the first basis function = 1
      // (boundary-constant extrapolation on the lower side).
      bval[lo_span] = 1.0;
    } else if (at_hi) {
      // Extrapolation above upper boundary: activate the last non-degenerate
      // span. After recursion this gives the last basis function = 1.
      bval[hi_span] = 1.0;
    } else {
      // Interior evaluation: standard degree-0 indicator.
      for (int j = 0; j < nk - 1; ++j) {
        double lo = knots[j];
        double hi = knots[j + 1];
        // Include upper boundary at the last knot span
        if (xi >= lo && (xi < hi || (j == nk - 2 && xi <= hi))) {
          bval[j] = 1.0;
        }
      }
    }

    // Recurse up to degree d (de Boor recursion)
    for (int d = 1; d <= degree; ++d) {
      int new_sz = nk - 1 - d;
      std::vector<double> bval_new(new_sz, 0.0);
      for (int j = 0; j < new_sz; ++j) {
        double left_span  = knots[j + d]     - knots[j];
        double right_span = knots[j + d + 1] - knots[j + 1];
        double t1 = 0.0, t2 = 0.0;
        if (left_span  > 0.0) t1 = (xi - knots[j])          / left_span  * bval[j];
        if (right_span > 0.0) t2 = (knots[j + d + 1] - xi)  / right_span * bval[j + 1];
        bval_new[j] = t1 + t2;
      }
      bval.resize(new_sz);
      bval = bval_new;
    }

    for (int j = 0; j < K; ++j)
      B(i, j) = bval[j];
  }

  return Rcpp::List::create(
    Rcpp::Named("B") = B,
    Rcpp::Named("K") = K
  );
}

// ============================================================================
// bgam_diff_penalty_cpp
//
// Builds the dpen-th order difference penalty matrix P = D'D of size K x K.
// First-order difference matrix D_1: D[i,i]=-1, D[i,i+1]=1, size (K-1) x K.
// Higher orders: D_dpen = D_1 * D_{dpen-1} (recursive composition).
//
// Returns: K x K symmetric matrix P = D'D.
// ============================================================================

// [[Rcpp::export]]
arma::mat bgam_diff_penalty_cpp(int K, int dpen) {
  if (dpen < 1)
    Rcpp::stop("bgam_diff_penalty_cpp: dpen must be >= 1.");
  if (K <= dpen)
    Rcpp::stop("bgam_diff_penalty_cpp: K must be > dpen.");

  int rows = K - 1;
  arma::mat D(rows, K, arma::fill::zeros);
  for (int i = 0; i < rows; ++i) {
    D(i, i)     = -1.0;
    D(i, i + 1) =  1.0;
  }

  for (int d = 2; d <= dpen; ++d) {
    int new_rows = (int)D.n_rows - 1;
    if (new_rows < 1)
      Rcpp::stop("bgam_diff_penalty_cpp: K too small for requested dpen.");
    arma::mat D1(new_rows, (int)D.n_rows, arma::fill::zeros);
    for (int i = 0; i < new_rows; ++i) {
      D1(i, i)     = -1.0;
      D1(i, i + 1) =  1.0;
    }
    D = D1 * D;
  }

  return D.t() * D;
}


// ============================================================================
// bgam_chol_with_jitter
//
// Wrap arma::chol with adaptive ridge retry. Some LAPACK builds (notably
// macOS Accelerate) are stricter on near-singular matrices than Linux/Windows.
// Retry up to 6 times with geometrically growing jitter on the diagonal.
// Returns false only if all retries fail.
// ============================================================================
static inline bool bgam_chol_with_jitter(arma::mat& L, const arma::mat& A) {
  if (arma::chol(L, A, "lower")) return true;
  const arma::uword n = A.n_rows;
  const double base = std::max(1e-12, arma::trace(A) / static_cast<double>(std::max<arma::uword>(1, n)));
  double jitter = base * 1e-8;
  arma::mat I = arma::eye<arma::mat>(n, n);
  for (int k = 0; k < 6; ++k) {
    if (arma::chol(L, A + jitter * I, "lower")) return true;
    jitter *= 100.0;
  }
  return false;
}
// ============================================================================
// bgam_prefactor_cpp
//
// Computes the Cholesky factorisation of A_j = B_j'B_j + lambda * P_j for
// the gaussian path (working weights constant = 1 or fixed prior weights).
// Called once per predictor at fit time; result is cached.
//
// Returns a List with:
//   $chol  -- K x K lower Cholesky factor L_j of A_j
//   $BtB   -- K x K, B_j'B_j
//   $Bt    -- K x n, B_j' (transposed design matrix)
// ============================================================================

// [[Rcpp::export]]
Rcpp::List bgam_prefactor_cpp(
    const arma::mat& B,
    const arma::mat& DtD,
    double lambda,
    Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue
) {
  // When w is non-NULL/non-empty, compute A_j = B_j'WB_j + lambda*DtD
  // (weighted path). When w is NULL (default), unweighted path.
  arma::mat BtB;
  if (w.isNotNull()) {
    arma::vec wv = Rcpp::as<arma::vec>(w.get());
    arma::vec sw = arma::sqrt(wv);
    arma::mat Bw = B.each_col() % sw;
    BtB = Bw.t() * Bw;
  } else {
    BtB = B.t() * B;
  }
  arma::mat A = BtB + lambda * DtD;

  arma::mat L;
  bool ok = bgam_chol_with_jitter(L, A);
  if (!ok)
    Rcpp::stop("bgam_prefactor_cpp: Cholesky failed (A_j not positive "
               "definite). Check that lambda > 0 and B_j has full column "
               "rank.");

  return Rcpp::List::create(
    Rcpp::Named("chol") = L,
    Rcpp::Named("BtB")  = BtB,
    Rcpp::Named("Bt")   = B.t()
  );
}

// ============================================================================
// BestResult struct for parallel_reduce accumulation over predictors.
//
// Each reduce slot holds (min_rss, best_j, best_beta). The join() method:
//   - keeps the entry with strictly lower rss, OR
//   - on equal rss, keeps the entry with the LOWER j index.
// This tie-break rule is applied in BOTH the serial scan (ascending j loop)
// and the parallel join, so the result is byte-identical across nthreads.
// ============================================================================

struct BestResult {
  double    min_rss;
  int       best_j;
  arma::vec best_beta;

  BestResult() : min_rss(std::numeric_limits<double>::infinity()),
                 best_j(-1) {}

  void consider(double rss, int j, const arma::vec& beta) {
    // Accept if strictly better, OR equal rss with lower j (tie-break).
    if (best_j < 0 ||
        rss < min_rss ||
        (rss == min_rss && j < best_j)) {
      min_rss   = rss;
      best_j    = j;
      best_beta = beta;
    }
  }

  // Merge another BestResult into this one using the same tie-break rule.
  void join(const BestResult& other) {
    if (other.best_j < 0) return;
    consider(other.min_rss, other.best_j, other.best_beta);
  }
};

// ============================================================================
// GaussianBoostWorker: RcppParallel::Worker for gaussian boosting step.
// Uses cached Cholesky factors (A_j does not change across steps).
//
// Splitting constructor takes RcppParallel::Split (not tbb::split) as
// required by RcppParallel::parallelReduce's ReducerWrapper.
// ============================================================================

struct GaussianBoostWorker : public RcppParallel::Worker {
  const Rcpp::List& B_list;
  const Rcpp::List& Bt_list;
  const Rcpp::List& chol_list;
  const arma::vec&  u;
  const arma::vec&  w_prior;   // prior weights (unit weights = unweighted)
  BestResult        result;

  GaussianBoostWorker(const Rcpp::List& B_list_,
                      const Rcpp::List& Bt_list_,
                      const Rcpp::List& chol_list_,
                      const arma::vec&  u_,
                      const arma::vec&  w_prior_)
    : B_list(B_list_), Bt_list(Bt_list_), chol_list(chol_list_), u(u_),
      w_prior(w_prior_) {}

  // RcppParallel splitting constructor -- uses RcppParallel::Split
  GaussianBoostWorker(const GaussianBoostWorker& other, RcppParallel::Split)
    : B_list(other.B_list), Bt_list(other.Bt_list),
      chol_list(other.chol_list), u(other.u),
      w_prior(other.w_prior) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t j = begin; j < end; ++j) {
      const arma::mat& Bt_j = Rcpp::as<arma::mat>(Bt_list[j]);
      const arma::mat& B_j  = Rcpp::as<arma::mat>(B_list[j]);
      const arma::mat& L_j  = Rcpp::as<arma::mat>(chol_list[j]);

      // w_prior is always a valid vector (unit weights if no prior weights).
      arma::vec rhs    = Bt_j * (w_prior % u);
      arma::vec beta_j = bgam_chol_solve(L_j, rhs);
      arma::vec f_j    = B_j * beta_j;
      arma::vec diff   = u - f_j;
      double    rss_j  = arma::dot(w_prior % diff, diff);

      result.consider(rss_j, (int)j, beta_j);
    }
  }

  void join(const GaussianBoostWorker& other) {
    result.join(other.result);
  }
};

// ============================================================================
// BinomialBoostWorker: RcppParallel::Worker for binomial boosting step.
// Refactors A_j^(w) = B_j'WB_j + lambda_j P_j at every step.
//
// Splitting constructor takes RcppParallel::Split.
// ============================================================================

struct BinomialBoostWorker : public RcppParallel::Worker {
  const Rcpp::List& B_list;
  const Rcpp::List& DtD_list;
  const arma::vec&  lambda_vec;
  const arma::vec&  u;
  const arma::vec&  w;
  BestResult        result;
  std::vector<int>  chol_failures;

  BinomialBoostWorker(const Rcpp::List& B_list_,
                      const Rcpp::List& DtD_list_,
                      const arma::vec&  lambda_vec_,
                      const arma::vec&  u_,
                      const arma::vec&  w_)
    : B_list(B_list_), DtD_list(DtD_list_),
      lambda_vec(lambda_vec_), u(u_), w(w_) {}

  // RcppParallel splitting constructor
  BinomialBoostWorker(const BinomialBoostWorker& other, RcppParallel::Split)
    : B_list(other.B_list), DtD_list(other.DtD_list),
      lambda_vec(other.lambda_vec), u(other.u), w(other.w) {}

  void operator()(std::size_t begin, std::size_t end) {
    arma::vec sw = arma::sqrt(w);   // sqrt(working weights)
    for (std::size_t j = begin; j < end; ++j) {
      const arma::mat& B_j   = Rcpp::as<arma::mat>(B_list[j]);
      const arma::mat& DtD_j = Rcpp::as<arma::mat>(DtD_list[j]);
      double lambda_j        = lambda_vec[j];

      // A_j^(w) = (sqrt(w) B_j)' (sqrt(w) B_j) + lambda_j P_j
      arma::mat Bw   = B_j.each_col() % sw;
      arma::mat A_jw = Bw.t() * Bw + lambda_j * DtD_j;

      arma::mat L_jw;
      bool ok = bgam_chol_with_jitter(L_jw, A_jw);
      if (!ok) {
        chol_failures.push_back((int)j);
        continue;
      }

      arma::vec rhs    = B_j.t() * (w % u);
      arma::vec beta_j = bgam_chol_solve(L_jw, rhs);
      arma::vec f_j    = B_j * beta_j;
      arma::vec diff   = u - f_j;
      double    wrss_j = arma::dot(w % diff, diff);

      result.consider(wrss_j, (int)j, beta_j);
    }
  }

  void join(const BinomialBoostWorker& other) {
    result.join(other.result);
    for (int jf : other.chol_failures)
      chol_failures.push_back(jf);
  }
};

// ============================================================================
// bgam_boost_cpp
//
// Main component-wise gradient boosting loop.
// family: 0 = gaussian, 1 = binomial.
// nthreads: forwarded to parallelReduce (0 or -1 = TBB default; 1 = serial).
//
// Returns a List with:
//   $selection_path  -- integer vector of length mstop (0-based j* per step)
//   $fitted          -- n-vector F after mstop steps
//   $beta_list       -- list of p accumulated (nu-scaled) beta vectors
//   $loss_path       -- numeric vector (length mstop): MSE or deviance/n
// ============================================================================

// [[Rcpp::export]]
Rcpp::List bgam_boost_cpp(
    const arma::vec& y,
    const arma::vec& F0,
    const Rcpp::List& B_list,
    const Rcpp::List& Bt_list,
    const Rcpp::List& chol_list,
    const Rcpp::List& DtD_list,
    const arma::vec& lambda_vec,
    double nu,
    int mstop,
    int family,
    int nthreads,
    Rcpp::Nullable<Rcpp::NumericVector> w_prior = R_NilValue
) {
  int n = (int)y.n_elem;
  int p = (int)B_list.size();

  if (p < 1)   Rcpp::stop("bgam_boost_cpp: B_list must have >= 1 element.");
  if (mstop < 1) Rcpp::stop("bgam_boost_cpp: mstop must be >= 1.");

  // Accumulated beta vectors (zero-initialised)
  std::vector<arma::vec> acc_beta(p);
  for (int j = 0; j < p; ++j) {
    const arma::mat& B_j = Rcpp::as<arma::mat>(B_list[j]);
    acc_beta[j].zeros(B_j.n_cols);
  }

  arma::vec F = F0;

  // Prior weights for gaussian WLS (unit weights when w_prior is NULL).
  arma::vec wp(n, arma::fill::ones);
  if (w_prior.isNotNull()) {
    wp = Rcpp::as<arma::vec>(w_prior.get());
  }

  // Initialise pseudo-residuals and (for binomial) working weights
  arma::vec u(n), w(n, arma::fill::ones);
  if (family == 0) {
    u = y - F;
  } else {
    for (int i = 0; i < n; ++i) {
      double mu_i = bgam_sigmoid(F[i]);
      mu_i = std::max(1e-6, std::min(1.0 - 1e-6, mu_i));
      u[i] = y[i] - mu_i;
      w[i] = mu_i * (1.0 - mu_i);
    }
  }

  Rcpp::IntegerVector selection_path(mstop);
  Rcpp::NumericVector loss_path(mstop);

  for (int m = 0; m < mstop; ++m) {
    int       jstar = -1;
    arma::vec beta_star;

    if (family == 0) {
      // ---- Gaussian path: use cached Cholesky ---------------------------
      if (nthreads == 1) {
        // Serial path: ascending j scan preserves tie-break by lowest j.
        double min_rss = std::numeric_limits<double>::infinity();
        for (int j = 0; j < p; ++j) {
          const arma::mat& Bt_j = Rcpp::as<arma::mat>(Bt_list[j]);
          const arma::mat& B_j  = Rcpp::as<arma::mat>(B_list[j]);
          const arma::mat& L_j  = Rcpp::as<arma::mat>(chol_list[j]);

          // wp is always valid (unit weights when no prior weights supplied).
          arma::vec rhs    = Bt_j * (wp % u);
          arma::vec beta_j = bgam_chol_solve(L_j, rhs);
          arma::vec f_j    = B_j * beta_j;
          arma::vec diff   = u - f_j;
          double    rss_j  = arma::dot(wp % diff, diff);

          // Strict less-than: ascending scan gives tie-break by lowest j.
          if (rss_j < min_rss) {
            min_rss   = rss_j;
            jstar     = j;
            beta_star = beta_j;
          }
        }
      } else {
        // Parallel path: parallelReduce with join()-level tie-break.
        GaussianBoostWorker worker(B_list, Bt_list, chol_list, u, wp);
        RcppParallel::parallelReduce(0, (std::size_t)p, worker,
                                     /*grainSize=*/1, nthreads);
        jstar     = worker.result.best_j;
        beta_star = worker.result.best_beta;
      }

    } else {
      // ---- Binomial path: refactor A_j^(w) per step --------------------
      if (nthreads == 1) {
        double min_wrss = std::numeric_limits<double>::infinity();
        arma::vec sw = arma::sqrt(w);
        for (int j = 0; j < p; ++j) {
          const arma::mat& B_j   = Rcpp::as<arma::mat>(B_list[j]);
          const arma::mat& DtD_j = Rcpp::as<arma::mat>(DtD_list[j]);
          double lambda_j        = lambda_vec[j];

          arma::mat Bw   = B_j.each_col() % sw;
          arma::mat A_jw = Bw.t() * Bw + lambda_j * DtD_j;
          arma::mat L_jw;
          bool ok = bgam_chol_with_jitter(L_jw, A_jw);
          if (!ok) {
            Rcpp::warning("bgam: base-learner %d singular at iteration %d; "
                          "skipping.", j, m + 1);
            continue;
          }

          arma::vec rhs    = B_j.t() * (w % u);
          arma::vec beta_j = bgam_chol_solve(L_jw, rhs);
          arma::vec f_j    = B_j * beta_j;
          arma::vec diff   = u - f_j;
          double    wrss_j = arma::dot(w % diff, diff);

          if (wrss_j < min_wrss) {
            min_wrss  = wrss_j;
            jstar     = j;
            beta_star = beta_j;
          }
        }
      } else {
        BinomialBoostWorker worker(B_list, DtD_list, lambda_vec, u, w);
        RcppParallel::parallelReduce(0, (std::size_t)p, worker,
                                     /*grainSize=*/1, nthreads);
        if (!worker.chol_failures.empty()) {
          for (int jf : worker.chol_failures)
            Rcpp::warning("bgam: base-learner %d singular at iteration %d; "
                          "skipping.", jf, m + 1);
        }
        jstar     = worker.result.best_j;
        beta_star = worker.result.best_beta;
      }
    }

    if (jstar < 0)
      Rcpp::stop("bgam_boost_cpp: all base-learners were singular at "
                 "iteration %d; cannot proceed.", m + 1);

    // ---- Serial update ---------------------------------------------------
    selection_path[m] = jstar;

    const arma::mat& B_jstar = Rcpp::as<arma::mat>(B_list[jstar]);
    acc_beta[jstar] += nu * beta_star;
    F += nu * (B_jstar * beta_star);

    if (family == 0) {
      u = y - F;
      loss_path[m] = arma::dot(u, u) / n;
    } else {
      double dev = 0.0;
      for (int i = 0; i < n; ++i) {
        double mu_i = bgam_sigmoid(F[i]);
        mu_i = std::max(1e-6, std::min(1.0 - 1e-6, mu_i));
        u[i] = y[i] - mu_i;
        w[i] = mu_i * (1.0 - mu_i);
        // Deviance contribution
        double mu_d = std::max(1e-10, std::min(1.0 - 1e-10, mu_i));
        dev -= 2.0 * (y[i] > 0.5 ? std::log(mu_d) : std::log(1.0 - mu_d));
      }
      loss_path[m] = dev / n;
    }
  }

  Rcpp::List beta_list_out(p);
  for (int j = 0; j < p; ++j)
    beta_list_out[j] = Rcpp::wrap(acc_beta[j]);

  return Rcpp::List::create(
    Rcpp::Named("selection_path") = selection_path,
    Rcpp::Named("fitted")         = Rcpp::wrap(F),
    Rcpp::Named("beta_list")      = beta_list_out,
    Rcpp::Named("loss_path")      = loss_path
  );
}

// ============================================================================
// bgam_predict_cpp
//
// Computes link-scale predictions on new data.
//
// Arguments:
//   B_new_list   -- p matrices (m_new x K_j) for new data
//   beta_list    -- p accumulated (nu-scaled) beta vectors
//   F0_scalar    -- scalar initialisation (mean(y) or logit(mean(y)))
//
// Returns: m_new-vector of fitted values on the link scale.
// F_new = F0_scalar + sum_{j: any(beta_j != 0)} B_new_j * beta_j
// ============================================================================

// [[Rcpp::export]]
arma::vec bgam_predict_cpp(
    const Rcpp::List& B_new_list,
    const Rcpp::List& beta_list,
    double F0_scalar
) {
  int p = (int)B_new_list.size();
  if (p < 1)
    Rcpp::stop("bgam_predict_cpp: B_new_list must not be empty.");

  const arma::mat& B0 = Rcpp::as<arma::mat>(B_new_list[0]);
  int m_new = (int)B0.n_rows;

  arma::vec F(m_new, arma::fill::value(F0_scalar));

  for (int j = 0; j < p; ++j) {
    const arma::vec& beta_j = Rcpp::as<arma::vec>(beta_list[j]);
    if (!arma::any(beta_j != 0.0)) continue;
    const arma::mat& B_j = Rcpp::as<arma::mat>(B_new_list[j]);
    F += B_j * beta_j;
  }

  return F;
}
