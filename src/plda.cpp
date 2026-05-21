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

// Condat (2013) direct 1-D total-variation denoising:
// minimize 0.5 * sum (b_j - u_j)^2 + lam * sum |b_{j+1} - b_j|.
// Reference: L. Condat, "A direct algorithm for 1-D total variation denoising,"
// IEEE Signal Processing Letters, 2013.
static arma::vec tv1d(const arma::vec& u, double lam) {
  const int N = (int) u.n_elem;
  arma::vec x(N);
  if (N == 0) return x;
  if (lam <= 0.0) { x = u; return x; }
  int k = 0, k0 = 0, kplus = 0, kminus = 0;
  double vmin = u[0] - lam, vmax = u[0] + lam;
  double umin = lam, umax = -lam;
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
static arma::vec mm_discriminant(const arma::mat& B, const arma::mat& Pmat,
                                 double lambda, double lambda2, int penalty,
                                 int maxit, double tol) {
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
        crits.back() < crits[crits.size() - 2] - 1e-6)
      Rcpp::stop("plda: minorize-maximize criterion decreased — numerical breakdown.");
  }
  return beta;
}

// [[Rcpp::export]]
List plda_fit_cpp(const arma::mat& x, const arma::ivec& y, int G, int K,
                  double lambda, double lambda2, int penalty,
                  int maxit, double tol) {
  if (K > G - 1) Rcpp::stop("plda_fit_cpp: K must be <= G-1 (G = %d)", G);
  if (K < 1) Rcpp::stop("plda_fit_cpp: K must be >= 1");
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
    arma::vec beta = mm_discriminant(B, Pmat, lambda, lambda2, penalty,
                                     maxit, tol);
    discrim.col(k) = beta;
  }
  return List::create(_["discrim"] = discrim, _["mu"] = S.mu,
                      _["sdw"] = S.sdw, _["cmeans"] = M, _["cw"] = w);
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
