// src/plda.cpp
#include <RcppArmadillo.h>
#include <RcppParallel.h>
// [[Rcpp::depends(RcppArmadillo, RcppParallel)]]
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

// Condat (2013) direct 1-D total-variation denoising:
// minimize 0.5 * sum (b_j - u_j)^2 + lam * sum |b_{j+1} - b_j|.
// Reference: L. Condat, "A direct algorithm for 1-D total variation denoising,"
// IEEE Signal Processing Letters, 2013.
static arma::vec tv1d(const arma::vec& u, double lam) {
  const int N = (int) u.n_elem;
  arma::vec x(N);
  if (N == 0) return x;
  if (lam <= 0.0) { x = u; return x; } // lam=0: no fusion penalty, the prox is the identity.
  int k = 0, k0 = 0, kplus = 0, kminus = 0;
  double vmin = u[0] - lam, vmax = u[0] + lam;
  double umin = lam, umax = -lam;
  // Algorithm state invariants (Condat's notation):
  //   [k0, k]  - current unresolved segment; x[k0..k] not yet written.
  //   vmin/vmax - current lower/upper bound on the pending output level.
  //   umin = sum_{j=k0..k}(u[j] - vmin) + lam   (cumulative min-side slack)
  //   umax = sum_{j=k0..k}(u[j] - vmax) - lam   (cumulative max-side slack)
  //   kminus = last k at which vmin was adjusted upward (umin clamped to lam).
  //   kplus  = last k at which vmax was adjusted downward (umax clamped to -lam).
  // Two terminal checks: one at the top of the loop (reset-to-last-element) and
  // one at the bottom (end-of-segment after ++k).
  while (true) {
    if (k == N - 1) {
      x[k] = vmin + umin;
      break;
    }
    if (u[k + 1] + umin < vmin - lam) {
      for (int i = k0; i <= kminus; ++i) x[i] = vmin;
      k = k0 = kminus = kplus = kminus + 1;
      vmin = u[k]; vmax = u[k] + 2 * lam;
      umin = lam; umax = -lam;
    } else if (u[k + 1] + umax > vmax + lam) {
      for (int i = k0; i <= kplus; ++i) x[i] = vmax;
      k = k0 = kminus = kplus = kplus + 1;
      vmin = u[k] - 2 * lam; vmax = u[k];
      umin = lam; umax = -lam;
    } else {
      ++k;
      umin += u[k] - vmin;
      umax += u[k] - vmax;
      if (umin >= lam) { vmin += (umin - lam) / (k - k0 + 1); umin = lam; kminus = k; }
      if (umax <= -lam) { vmax += (umax + lam) / (k - k0 + 1); umax = -lam; kplus = k; }
    }
    if (k == N - 1) {
      if (umin < 0.0) {
        for (int i = k0; i <= kminus; ++i) x[i] = vmin;
        k = k0 = kminus = kminus + 1;
        vmin = u[k]; umin = lam; umax = u[k] + lam - vmax;
      } else if (umax > 0.0) {
        for (int i = k0; i <= kplus; ++i) x[i] = vmax;
        k = k0 = kplus = kplus + 1;
        vmax = u[k]; umax = -lam; umin = u[k] - lam - vmin;
      } else {
        for (int i = k0; i <= N - 1; ++i) x[i] = vmin + umin / (k - k0 + 1);
        break;
      }
    }
  }
  return x;
}

// [[Rcpp::export]]
Rcpp::NumericVector plda_tv1d_cpp(const arma::vec& u, double lam) {
  arma::vec v = tv1d(u, lam);
  return Rcpp::NumericVector(v.begin(), v.end());
}

// Standardize x: subtract global column means, divide by within-class sd.
// Drops constant features by leaving their sd as 1 (column becomes ~0).
struct Standardized { arma::mat xs; arma::rowvec mu; arma::vec sdw; };

static Standardized standardize(const arma::mat& x, const arma::ivec& y, int G) {
  Standardized S;
  S.mu  = arma::mean(x, 0);
  S.sdw = wcsd_impl(x, y, G);
  S.sdw.transform([](double v) { return v > 1e-12 ? v : 1.0; });
  S.xs = x;
  S.xs.each_row() -= S.mu;
  S.xs.each_row() /= S.sdw.t();
  return S;
}

// G x p matrix of class means of the standardized data, and class weights n_g/n.
static void class_means(const arma::mat& xs, const arma::ivec& y, int G,
                        arma::mat& M, arma::vec& w) {
  const arma::uword n = xs.n_rows;
  M.set_size(G, xs.n_cols); w.set_size(G);
  for (int g = 1; g <= G; ++g) {
    arma::uvec idx = arma::find(y == g);
    M.row(g - 1) = arma::mean(xs.rows(idx), 0);
    w(g - 1) = (double) idx.n_elem / (double) n;
  }
}

// Penalized-PCA objective: t(beta) B' P B beta - d * P(beta).  (P here = lambda
// penalty, distinct from the projection matrix Pmat.)  Matches penalizedLDA's
// PenalizedPCACrit; used for the criterion-based convergence test.
static double ppca_crit(const arma::mat& BtPB, const arma::vec& beta, double d,
                        double lambda, double lambda2, int penalty) {
  double quad = arma::as_scalar(beta.t() * BtPB * beta);
  double pen;
  if (penalty == 0) {
    pen = lambda * arma::sum(arma::abs(beta));
  } else {
    pen = lambda * arma::sum(arma::abs(beta))
        + lambda2 * arma::sum(arma::abs(arma::diff(beta)));
  }
  return quad - d * pen;
}

// One penalized discriminant via minorize-maximize, following penalizedLDA's
// PenalizedPCA inner loop.  B is the (G x p) between-class root-scatter matrix
// sqrt.sigma.bet; Pmat is the (G x G) deflation projection (identity for k=1).
// prox: 0 = L1 soft-threshold; 1 = fused (tv1d then soft-threshold).
//
// `ok` is an out-parameter: it is set to false (and the loop bailed) if the MM
// criterion ever decreases — a numerical breakdown. Callers running on the main
// R thread translate `ok == false` into Rcpp::stop(); callers running inside a
// TBB worker (the CV harness) must NOT throw across the thread boundary and so
// inspect `ok` after the parallelFor completes. `Rcpp::stop` is intentionally
// never reached from this function for that reason.
static arma::vec mm_discriminant(const arma::mat& B, const arma::mat& Pmat,
                                 double lambda, double lambda2, int penalty,
                                 int maxit, double tol, bool& ok) {
  // BtP = t(x) %*% P  (p x G); svd gives d and the warm-start direction.
  arma::mat BtP = B.t() * Pmat;
  arma::mat U, V; arma::vec s;
  arma::svd(U, s, V, BtP);
  double d = (s.n_elem > 0) ? s[0] * s[0] : 0.0;
  arma::vec beta = (U.n_cols > 0) ? arma::vec(U.col(0))
                                  : arma::vec(B.n_cols, arma::fill::zeros);
  arma::mat BtPB = BtP * B;                      // p x p effective scatter

  std::vector<double> crits;
  crits.push_back(ppca_crit(BtPB, beta, d, lambda, lambda2, penalty));
  for (int iter = 0; iter < maxit; ++iter) {
    // crits needs >= 4 entries before the relative-change test fires (matches penalizedLDA).
    bool converged = (crits.size() >= 4) &&
      (std::abs(crits.back() - crits[crits.size() - 2]) /
         std::max(0.001, crits.back()) <= tol);
    if (converged || arma::accu(arma::abs(beta)) == 0.0) break;
    arma::vec tmp = BtPB * beta;
    arma::vec prox;
    // penalizedLDA scales the L1 threshold by d/2, where d is the leading squared
    // singular value of the projected between-class root-scatter matrix.
    if (penalty == 0) prox = soft_threshold(tmp, d * lambda / 2.0);
    else              prox = soft_threshold(tv1d(tmp, d * lambda2 / 2.0),
                                            d * lambda / 2.0);
    beta = normalize_l2(prox);
    beta.replace(arma::datum::nan, 0.0);
    crits.push_back(ppca_crit(BtPB, beta, d, lambda, lambda2, penalty));
    if (crits.size() >= 2 &&
        crits.back() < crits[crits.size() - 2] - 1e-6) {
      ok = false;
      return beta;
    }
  }
  return beta;
}

// Engine output of a single pLDA fit. `ok` is false if the MM inner loop hit a
// numerical breakdown; the caller decides whether to throw (main thread) or
// flag the result (TBB worker).
struct PldaFit {
  arma::mat discrim;   // p x K discriminant vectors
  arma::rowvec mu;     // global feature means
  arma::vec sdw;       // within-class feature sds
  arma::mat cmeans;    // G x p standardized class means
  arma::vec cw;        // class weights n_g / n
  bool ok = true;
};

// Core pLDA fit. Pure numeric — calls no Rcpp throwing primitives, so it is
// safe to invoke from inside a TBB worker. Determinism: every step is a fixed,
// data-only computation; running it on N threads (one fit per thread) and on 1
// thread produces byte-identical PldaFit results because the fits are wholly
// independent and never share mutable state.
static PldaFit plda_fit_core(const arma::mat& x, const arma::ivec& y,
                             int G, int K, double lambda, double lambda2,
                             int penalty, int maxit, double tol) {
  PldaFit out;
  Standardized S = standardize(x, y, G);
  arma::mat M; arma::vec w;
  class_means(S.xs, y, G, M, w);

  // Between-class root-scatter matrix B = sqrt.sigma.bet = diag(sqrt(w)) %*% M
  // (G x p); B' B equals the matrix-free Sigma_b = M' diag(w) M.
  arma::mat B = M;
  B.each_col() %= arma::sqrt(w);

  arma::mat discrim(S.xs.n_cols, K, arma::fill::zeros);
  for (int k = 0; k < K; ++k) {
    arma::mat Pmat;
    if (k == 0) {
      Pmat = arma::eye<arma::mat>(G, G);
    } else {
      // Deflate in G-space: P = I - U U', U = left sing. vecs of B %*% betas.
      arma::mat proj = B * discrim.cols(0, k - 1);
      arma::mat U, V; arma::vec s;
      arma::svd(U, s, V, proj);
      arma::uvec keep = arma::find(s > 1e-10);
      Pmat = arma::eye<arma::mat>(G, G);
      if (!keep.is_empty()) {
        arma::mat Uk = U.cols(keep);
        Pmat -= Uk * Uk.t();
      }
    }
    bool ok = true;
    arma::vec beta = mm_discriminant(B, Pmat, lambda, lambda2, penalty,
                                     maxit, tol, ok);
    if (!ok) out.ok = false;
    discrim.col(k) = beta;
  }
  out.discrim = discrim;
  out.mu = S.mu;
  out.sdw = S.sdw;
  out.cmeans = M;
  out.cw = w;
  return out;
}

// [[Rcpp::export]]
List plda_fit_cpp(const arma::mat& x, const arma::ivec& y, int G, int K,
                  double lambda, double lambda2, int penalty,
                  int maxit, double tol) {
  if (K > G - 1) Rcpp::stop("plda_fit_cpp: K must be <= G-1 (G = %d)", G);
  if (K < 1) Rcpp::stop("plda_fit_cpp: K must be >= 1");
  PldaFit f = plda_fit_core(x, y, G, K, lambda, lambda2, penalty, maxit, tol);
  if (!f.ok)
    Rcpp::stop("plda: minorize-maximize criterion decreased — numerical breakdown.");
  return List::create(_["discrim"] = f.discrim, _["mu"] = f.mu,
                      _["sdw"] = f.sdw, _["cmeans"] = f.cmeans, _["cw"] = f.cw);
}

// Project new data onto stored discriminant vectors after standardizing
// with the training mu/sdw. Returns an m x K score matrix.
// [[Rcpp::export]]
arma::mat plda_project_cpp(const arma::mat& xnew, const arma::rowvec& mu,
                           const arma::vec& sdw, const arma::mat& discrim) {
  arma::mat xs = xnew;
  xs.each_row() -= mu;
  xs.each_row() /= sdw.t();
  return xs * discrim;
}

// ---------------------------------------------------------------------------
// Parallel k-fold cross-validation harness (pLDA Task 13).
//
// Runs the full nfold x length(lambda_grid) grid. One TBB task == one fold.
// Each fold's misclassification counts go into a private slot (err_by_fold[f],
// a length(lambda_grid) x K matrix). After the parallelFor, the slots are
// summed in fixed fold order 0,1,...,nfold-1 — so the reduction is byte-
// identical regardless of how many threads ran the folds, satisfying the
// roadrunner determinism invariant.
//
// Determinism notes:
//   * Folds are assigned in R (set.seed(0) + sample(...)) and passed in as a
//     1-based integer vector `folds`; C++ never randomizes.
//   * Each task only reads shared const inputs and writes its own err_by_fold
//     slot — no cross-thread mutable aliasing.
//   * plda_fit_core is pure numeric; one fit per thread is independent of one
//     fit on the main thread, bit-for-bit.
//   * The final err matrix is sum_f err_by_fold[f] accumulated in ascending f.
//
// Classification mirrors the former R loop exactly: for prefix dimension k,
// nearest-centroid in the first k discriminant scores, ties broken "first"
// (the lowest class index), via an argmin over squared distances.
// ---------------------------------------------------------------------------

// Misclassification counts for one fold, all lambdas, all k in 1..K.
// Pure numeric: safe inside a TBB worker.
static void plda_cv_one_fold(const arma::mat& x, const arma::ivec& y,
                             const arma::uvec& tr, const arma::uvec& te,
                             int G, int K, const arma::vec& lambda_grid,
                             double lambda2, int penalty, int maxit, double tol,
                             arma::mat& err_slot, bool& ok_slot) {
  const arma::uword nlam = lambda_grid.n_elem;
  arma::mat x_tr = x.rows(tr);
  arma::ivec y_tr = y.elem(tr);
  arma::mat x_te = x.rows(te);
  arma::ivec y_te = y.elem(te);
  const arma::uword nte = te.n_elem;

  for (arma::uword li = 0; li < nlam; ++li) {
    PldaFit f = plda_fit_core(x_tr, y_tr, G, K, lambda_grid(li), lambda2,
                              penalty, maxit, tol);
    if (!f.ok) ok_slot = false;

    // Test-point scores (nte x K) and class-mean scores (G x K).
    arma::mat sc  = plda_project_cpp(x_te, f.mu, f.sdw, f.discrim);
    arma::mat csc = f.cmeans * f.discrim;

    for (int k = 1; k <= K; ++k) {
      // Squared distance in the first k score dimensions; nearest centroid.
      arma::uword miscls = 0;
      for (arma::uword i = 0; i < nte; ++i) {
        double best = arma::datum::inf;
        int bestg = 1;
        for (int g = 0; g < G; ++g) {
          double d2 = 0.0;
          for (int kk = 0; kk < k; ++kk) {
            double diff = sc(i, kk) - csc(g, kk);
            d2 += diff * diff;
          }
          // ties.method = "first": strict < keeps the lowest class index.
          if (d2 < best) { best = d2; bestg = g + 1; }
        }
        if (bestg != y_te(i)) ++miscls;
      }
      err_slot(li, k - 1) = (double) miscls;
    }
  }
}

struct PldaCvWorker : public RcppParallel::Worker {
  const arma::mat& x;
  const arma::ivec& y;
  const std::vector<arma::uvec>& tr_idx;
  const std::vector<arma::uvec>& te_idx;
  int G, K;
  const arma::vec& lambda_grid;
  double lambda2;
  int penalty, maxit;
  double tol;
  std::vector<arma::mat>& err_by_fold;   // one fixed slot per fold
  std::vector<unsigned char>& ok_by_fold;

  PldaCvWorker(const arma::mat& x_, const arma::ivec& y_,
               const std::vector<arma::uvec>& tr_, const std::vector<arma::uvec>& te_,
               int G_, int K_, const arma::vec& lg_, double lambda2_,
               int penalty_, int maxit_, double tol_,
               std::vector<arma::mat>& err_, std::vector<unsigned char>& ok_)
    : x(x_), y(y_), tr_idx(tr_), te_idx(te_), G(G_), K(K_),
      lambda_grid(lg_), lambda2(lambda2_), penalty(penalty_), maxit(maxit_),
      tol(tol_), err_by_fold(err_), ok_by_fold(ok_) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t f = begin; f < end; ++f) {
      bool ok = true;
      plda_cv_one_fold(x, y, tr_idx[f], te_idx[f], G, K, lambda_grid,
                       lambda2, penalty, maxit, tol, err_by_fold[f], ok);
      ok_by_fold[f] = ok ? 1 : 0;
    }
  }
};

// [[Rcpp::export]]
Rcpp::List plda_cv_inner_cpp(const arma::mat& x, const arma::ivec& y,
                             const arma::ivec& folds, int nfold,
                             int G, int K, const arma::vec& lambda_grid,
                             double lambda2, int penalty,
                             int maxit, double tol, int nthreads) {
  const arma::uword nlam = lambda_grid.n_elem;

  // Pre-resolve per-fold train/test row index vectors (1-based folds in R).
  std::vector<arma::uvec> tr_idx(nfold), te_idx(nfold);
  for (int f = 0; f < nfold; ++f) {
    tr_idx[f] = arma::find(folds != (f + 1));
    te_idx[f] = arma::find(folds == (f + 1));
  }

  // Disjoint output slots: one error matrix and one ok flag per fold.
  std::vector<arma::mat> err_by_fold(nfold);
  for (int f = 0; f < nfold; ++f)
    err_by_fold[f] = arma::mat(nlam, K, arma::fill::zeros);
  std::vector<unsigned char> ok_by_fold(nfold, 1);

  // Clamp worker count: >= 1, no idle workers beyond the number of folds.
  int n_workers = nthreads;
  if (n_workers < 1) n_workers = 1;
  if (n_workers > nfold) n_workers = nfold;

  PldaCvWorker worker(x, y, tr_idx, te_idx, G, K, lambda_grid, lambda2,
                      penalty, maxit, tol, err_by_fold, ok_by_fold);
  RcppParallel::parallelFor(0, nfold, worker, /*grainSize=*/1,
                            /*numThreads=*/n_workers);

  // Serial-order reduction: sum the per-fold slots in ascending fold index so
  // the result is byte-identical for any n_workers.
  arma::mat err(nlam, K, arma::fill::zeros);
  bool all_ok = true;
  for (int f = 0; f < nfold; ++f) {
    err += err_by_fold[f];
    if (!ok_by_fold[f]) all_ok = false;
  }

  Rcpp::NumericMatrix err_out(nlam, K);
  for (arma::uword i = 0; i < nlam; ++i)
    for (int k = 0; k < K; ++k)
      err_out(i, k) = err(i, k);

  return Rcpp::List::create(
    Rcpp::Named("err")          = err_out,
    Rcpp::Named("ok")           = all_ok,
    Rcpp::Named("nthreads_used") = n_workers
  );
}
