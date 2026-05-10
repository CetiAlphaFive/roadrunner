// ares — Fast Multivariate Adaptive Regression Splines
// Author: Jack Trametta <jtrametta@gmail.com>
// License: MIT
// References: Friedman (1991), "Multivariate Adaptive Regression Splines",
//             Annals of Statistics 19(1):1-67.

// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::plugins(cpp17)]]

#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
#include <numeric>

using namespace Rcpp;

namespace ares {

// ----- Helpers ----------------------------------------------------------------

inline double safe_div(double a, double b, double eps = 1e-300) {
  return std::fabs(b) < eps ? 0.0 : a / b;
}

// Auto minspan / endspan (Friedman 1991 / earth conventions)
static int auto_minspan(int p, int n) {
  if (p <= 0) p = 1;
  if (n <= 0) return 1;
  double v = -std::log2(-(1.0 / (double(p) * double(n))) * std::log(1.0 - 0.05)) / 2.5;
  int ms = int(std::floor(v));
  if (ms < 1) ms = 1;
  return ms;
}

static int auto_endspan(int p) {
  if (p <= 0) p = 1;
  double v = std::floor(3.0 - std::log2(0.05 / double(p)));
  int es = int(v);
  if (es < 1) es = 1;
  return es;
}

// Hinge function h(s, x, t) = max(0, s*(x - t))
inline double hinge(int s, double x, double t) {
  double v = (s == 1) ? (x - t) : (t - x);
  return v > 0.0 ? v : 0.0;
}

// Build basis column for a single term given dirs/cuts row and selected x
static void build_term_column(const double* x, int n, int p,
                              const int* dirs_row, const double* cuts_row,
                              double* out) {
  for (int i = 0; i < n; ++i) out[i] = 1.0;
  for (int j = 0; j < p; ++j) {
    int d = dirs_row[j];
    if (d == 0) continue;
    double t = cuts_row[j];
    if (d == 2) {
      for (int i = 0; i < n; ++i) out[i] *= x[i + size_t(j) * n];
    } else {
      int s = d;
      for (int i = 0; i < n; ++i) {
        double v = (s == 1) ? (x[i + size_t(j) * n] - t) : (t - x[i + size_t(j) * n]);
        out[i] *= (v > 0.0 ? v : 0.0);
      }
    }
  }
}

// OLS least squares via QR (LAPACK dgels-style via Eigen-free hand-rolled QR using
// Householder reflections). For modest M, plain QR is fine and avoids extra deps.
// Returns RSS; fills beta (length M).
static double ols_qr(const double* B, int n, int M, const double* y,
                     std::vector<double>& beta, std::vector<double>& fitted) {
  // Copy B into a column-major working array A(n x M) and y into b
  std::vector<double> A(size_t(n) * M);
  for (int j = 0; j < M; ++j)
    for (int i = 0; i < n; ++i)
      A[i + size_t(j)*n] = B[i + size_t(j)*n];
  std::vector<double> b(y, y + n);

  // Householder QR with column pivoting — for simplicity, no pivoting (rank check via tiny diag)
  std::vector<double> tau(M, 0.0);
  for (int k = 0; k < M; ++k) {
    // Householder vector for column k below diagonal
    double norm_x = 0.0;
    for (int i = k; i < n; ++i) norm_x += A[i + size_t(k)*n] * A[i + size_t(k)*n];
    norm_x = std::sqrt(norm_x);
    if (norm_x == 0.0) { tau[k] = 0.0; continue; }
    double sign = A[k + size_t(k)*n] >= 0 ? 1.0 : -1.0;
    double alpha = -sign * norm_x;
    A[k + size_t(k)*n] -= alpha;
    double v_norm_sq = 0.0;
    for (int i = k; i < n; ++i) v_norm_sq += A[i + size_t(k)*n] * A[i + size_t(k)*n];
    if (v_norm_sq == 0.0) { tau[k] = 0.0; continue; }
    tau[k] = 2.0 / v_norm_sq;
    // Apply Householder to remaining columns
    for (int j = k + 1; j < M; ++j) {
      double dot = 0.0;
      for (int i = k; i < n; ++i) dot += A[i + size_t(k)*n] * A[i + size_t(j)*n];
      dot *= tau[k];
      for (int i = k; i < n; ++i) A[i + size_t(j)*n] -= dot * A[i + size_t(k)*n];
    }
    // Apply to b
    double dot_b = 0.0;
    for (int i = k; i < n; ++i) dot_b += A[i + size_t(k)*n] * b[i];
    dot_b *= tau[k];
    for (int i = k; i < n; ++i) b[i] -= dot_b * A[i + size_t(k)*n];
    // Restore A(k,k) = alpha (but we use the diag of R we just computed)
    A[k + size_t(k)*n] = alpha;
  }
  // Back-substitute: R beta = b[0:M]
  beta.assign(M, 0.0);
  for (int k = M - 1; k >= 0; --k) {
    double sum = b[k];
    for (int j = k + 1; j < M; ++j) sum -= A[k + size_t(j)*n] * beta[j];
    double dkk = A[k + size_t(k)*n];
    if (std::fabs(dkk) < 1e3 * std::numeric_limits<double>::epsilon() * (1.0 + std::fabs(sum))) {
      beta[k] = 0.0;       // rank-deficient; pseudo-zero
    } else {
      beta[k] = sum / dkk;
    }
  }
  // Compute fitted = B beta (n vector) and rss
  fitted.assign(n, 0.0);
  for (int j = 0; j < M; ++j) {
    double bj = beta[j];
    if (bj == 0.0) continue;
    for (int i = 0; i < n; ++i) fitted[i] += bj * B[i + size_t(j)*n];
  }
  double rss = 0.0;
  for (int i = 0; i < n; ++i) {
    double r = y[i] - fitted[i];
    rss += r * r;
  }
  return rss;
}

// ----- Forward-pass parallel worker ------------------------------------------
struct Candidate {
  double rss_red = -1.0;        // RSS reduction (larger is better)
  int parent     = -1;
  int var        = -1;
  double cut     =  0.0;
};

inline bool better(const Candidate& a, const Candidate& b) {
  // larger rss_red wins; tie-break: smaller var, then smaller cut, then smaller parent
  if (a.rss_red != b.rss_red) return a.rss_red > b.rss_red;
  if (a.var     != b.var    ) return a.var     < b.var;
  if (a.cut     != b.cut    ) return a.cut     < b.cut;
  return a.parent < b.parent;
}

// For a single (parent, var) pair, scan all eligible knots and compute the best
// RSS reduction from adding the *pair* of hinges (h+, h-) multiplied by the
// parent column. Uses Friedman fast-LS scoring: each candidate column is
// orthogonalized against the FULL current basis Q (orthonormal columns) and
// then the joint reduction from the pair is computed using the proper
// Gram-Schmidt formula.
struct KnotScanner {
  // inputs
  const double* parent_col;     // n
  const double* xj;             // n
  const double* resid;          // n  (current residual: y - B beta, already orthogonal to span(Q))
  const double* Qmat;           // n*M_q column-major (orthonormal basis of B)
  int n;
  int Mq;
  int minspan;
  int endspan;
  // output
  Candidate best;

  void run(int parent_idx, int var_idx) {
    // Sort indices of rows where parent_col != 0 by xj
    std::vector<int> idx;
    idx.reserve(n);
    for (int i = 0; i < n; ++i) if (parent_col[i] != 0.0) idx.push_back(i);
    int n_eli = (int)idx.size();
    if (n_eli < 2 * endspan + 3) return;
    std::stable_sort(idx.begin(), idx.end(),
              [this](int a, int b){
                if (xj[a] != xj[b]) return xj[a] < xj[b];
                return a < b;  // stable tie-break by row index
              });

    // Candidate knots: positions endspan .. n_eli - endspan - 1, with minspan stride;
    // unique x values only.
    std::vector<int> knot_pos;
    int last_pos = -1;
    for (int k = endspan; k < n_eli - endspan; ++k) {
      if (last_pos >= 0 && (k - last_pos) < minspan) continue;
      if (k > 0 && xj[idx[k]] == xj[idx[k - 1]]) continue;
      knot_pos.push_back(k);
      last_pos = k;
    }
    if (knot_pos.empty()) return;

    // For each candidate knot t, compute RSS reduction from adding pair
    //   c+ = parent_col * max(0, xj - t),  c- = parent_col * max(0, t - xj).
    // Both columns are orthogonalized against Q (full basis). Then pair-joint
    // reduction = (c+_perp . r)^2/||c+_perp||^2 + (c-_perp_after_cp . r)^2/||c-_perp_after_cp||^2.
    Candidate local_best;
    std::vector<double> cp(n, 0.0), cm(n, 0.0);
    std::vector<double> qcp(Mq, 0.0), qcm(Mq, 0.0);
    for (int kp_idx = 0; kp_idx < (int)knot_pos.size(); ++kp_idx) {
      int kp = knot_pos[kp_idx];
      double t = xj[idx[kp]];
      // Build c+ and c- as full n-vectors (zero outside parent's support)
      double cp_dot_r = 0.0, cm_dot_r = 0.0;
      double cp_norm2 = 0.0, cm_norm2 = 0.0;
      for (int i = 0; i < n; ++i) {
        double pi = parent_col[i];
        if (pi == 0.0) { cp[i] = 0.0; cm[i] = 0.0; continue; }
        double xi = xj[i];
        double hp = (xi > t) ? pi * (xi - t) : 0.0;
        double hm = (xi < t) ? pi * (t - xi) : 0.0;
        cp[i] = hp; cm[i] = hm;
        cp_dot_r += hp * resid[i];
        cm_dot_r += hm * resid[i];
        cp_norm2 += hp * hp;
        cm_norm2 += hm * hm;
      }
      if (cp_norm2 < 1e-300 && cm_norm2 < 1e-300) continue;
      // Compute Q' c+ and Q' c-
      for (int q = 0; q < Mq; ++q) {
        const double* Qcol = Qmat + size_t(q) * n;
        double dp = 0.0, dm = 0.0;
        for (int i = 0; i < n; ++i) { dp += Qcol[i] * cp[i]; dm += Qcol[i] * cm[i]; }
        qcp[q] = dp; qcm[q] = dm;
      }
      // ||c+_perp||^2 = ||c+||^2 - sum_q (Q'c+)_q^2
      double cp_perp_norm2 = cp_norm2;
      double cm_perp_norm2 = cm_norm2;
      for (int q = 0; q < Mq; ++q) {
        cp_perp_norm2 -= qcp[q] * qcp[q];
        cm_perp_norm2 -= qcm[q] * qcm[q];
      }
      // c+_perp . c-_perp
      double cpcm = 0.0;
      for (int i = 0; i < n; ++i) cpcm += cp[i] * cm[i];
      // (note: cp .* cm == 0 because hinge supports are disjoint, so cpcm == 0
      //  in exact arithmetic; we still compute the perp inner product:
      //  c+_perp . c-_perp = c+ . c- - sum_q (Q'c+)_q (Q'c-)_q)
      double cp_perp_dot_cm_perp = cpcm;
      for (int q = 0; q < Mq; ++q) cp_perp_dot_cm_perp -= qcp[q] * qcm[q];

      // Joint reduction via 2-step Gram-Schmidt:
      // step 1: e1 = c+_perp / ||c+_perp||;  reduction1 = (cp_dot_r)^2 / cp_perp_norm2
      // step 2: e2 = c-_perp - (e1 . c-_perp) e1;  ||e2||^2 = cm_perp_norm2 - (cp_perp_dot_cm_perp)^2 / cp_perp_norm2
      // reduction2 = (e2 . r)^2 / ||e2||^2 where e2 . r = cm_dot_r - (cp_perp_dot_cm_perp / cp_perp_norm2) * cp_dot_r
      double red_p = 0.0, red_m = 0.0;
      const double tiny = 1e-12;
      if (cp_perp_norm2 > tiny) {
        red_p = cp_dot_r * cp_dot_r / cp_perp_norm2;
        // e2 norm and projection
        double alpha = cp_perp_dot_cm_perp / cp_perp_norm2;
        double e2_norm2 = cm_perp_norm2 - alpha * cp_perp_dot_cm_perp;
        double e2_dot_r = cm_dot_r - alpha * cp_dot_r;
        if (e2_norm2 > tiny) red_m = e2_dot_r * e2_dot_r / e2_norm2;
      } else if (cm_perp_norm2 > tiny) {
        // Only "-" side is non-degenerate
        red_m = cm_dot_r * cm_dot_r / cm_perp_norm2;
      }
      double red = red_p + red_m;
      if (red > local_best.rss_red ||
          (red == local_best.rss_red && var_idx < local_best.var) ||
          (red == local_best.rss_red && var_idx == local_best.var && t < local_best.cut)) {
        local_best.rss_red = red;
        local_best.parent  = parent_idx;
        local_best.var     = var_idx;
        local_best.cut     = t;
      }
    }
    if (better(local_best, best)) best = local_best;
  }
};

struct ForwardWorker : public RcppParallel::Worker {
  const double* x;             // n*p flat column-major
  const double* B_cols;         // n*M_active flat column-major (parent columns, in same order as parent_global_idx)
  const double* resid;          // n
  const double* Qmat;           // n*Mq orthonormal basis of current B
  int Mq;
  const std::vector<std::pair<int,int>>* pairs;  // (parent_local_idx, var_idx)
  const std::vector<int>* parent_global_idx;
  int n, p;
  int minspan, endspan;
  std::vector<Candidate>* out;  // size = pairs->size()

  ForwardWorker(const double* x_, const double* B_, const double* r_,
                const double* Q_, int Mq_,
                const std::vector<std::pair<int,int>>& pairs_,
                const std::vector<int>& parent_global_,
                int n_, int p_, int ms_, int es_,
                std::vector<Candidate>& out_)
    : x(x_), B_cols(B_), resid(r_), Qmat(Q_), Mq(Mq_),
      pairs(&pairs_), parent_global_idx(&parent_global_),
      n(n_), p(p_), minspan(ms_), endspan(es_), out(&out_) {}

  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t i = begin; i < end; ++i) {
      int parent_local = (*pairs)[i].first;
      int var_idx     = (*pairs)[i].second;
      KnotScanner sc{
        B_cols + size_t(parent_local) * n,
        x      + size_t(var_idx)     * n,
        resid, Qmat, n, Mq, minspan, endspan, Candidate{}
      };
      sc.run((*parent_global_idx)[parent_local], var_idx);
      (*out)[i] = sc.best;
    }
  }
};

// ----- Main fit --------------------------------------------------------------

static void recompute_residual(const double* B, int n, int M,
                               const double* y,
                               std::vector<double>& beta,
                               std::vector<double>& fitted,
                               std::vector<double>& resid,
                               double& rss) {
  rss = ols_qr(B, n, M, y, beta, fitted);
  resid.resize(n);
  for (int i = 0; i < n; ++i) resid[i] = y[i] - fitted[i];
}

// Compute orthonormal Q (n x M) of B (n x M) via Modified Gram-Schmidt with
// pivoting-free dropping of near-zero columns. Returns the actual rank Mq <= M.
// Q is laid out column-major in q_out (size n*M, only first Mq columns valid).
static int build_Q(const double* B, int n, int M, double* q_out) {
  int Mq = 0;
  std::vector<double> v(n);
  const double tol = 1e3 * std::numeric_limits<double>::epsilon();
  for (int j = 0; j < M; ++j) {
    for (int i = 0; i < n; ++i) v[i] = B[i + size_t(j) * n];
    // Project out previously accepted q-columns
    for (int k = 0; k < Mq; ++k) {
      const double* qk = q_out + size_t(k) * n;
      double dot = 0.0;
      for (int i = 0; i < n; ++i) dot += qk[i] * v[i];
      for (int i = 0; i < n; ++i) v[i] -= dot * qk[i];
    }
    double norm2 = 0.0;
    for (int i = 0; i < n; ++i) norm2 += v[i] * v[i];
    double norm = std::sqrt(norm2);
    // Reference scale: max abs of column entries
    double col_scale = 0.0;
    for (int i = 0; i < n; ++i) col_scale = std::max(col_scale, std::fabs(B[i + size_t(j)*n]));
    if (norm < tol * (col_scale + 1.0)) continue;  // rank-deficient — skip
    double inv = 1.0 / norm;
    double* qj = q_out + size_t(Mq) * n;
    for (int i = 0; i < n; ++i) qj[i] = v[i] * inv;
    ++Mq;
  }
  return Mq;
}

// Count how many variables a term uses (interaction depth)
static int term_degree(const int* dirs_row, int p) {
  int d = 0;
  for (int j = 0; j < p; ++j) if (dirs_row[j] != 0) ++d;
  return d;
}

// Check if a term uses variable j
static bool term_uses_var(const int* dirs_row, int j) {
  return dirs_row[j] != 0;
}

} // namespace ares

// =============================================================================
// Public R-callable entries
// =============================================================================

// [[Rcpp::export]]
List mars_fit_cpp(const NumericMatrix& x_in,
                  const NumericVector& y_in,
                  int degree, int nk, double penalty, double thresh,
                  int minspan_in, int endspan_in, int nprune,
                  int pmethod, int trace, int nthreads) {
  using namespace ares;
  int n = x_in.nrow();
  int p = x_in.ncol();
  if (n < 3) stop("ares: need at least 3 observations.");
  if (p < 1) stop("ares: x must have at least one column.");
  if (y_in.size() != n) stop("ares: length(y) must equal nrow(x).");
  if (nk < 3) nk = 3;
  if (degree < 1) stop("ares: degree must be >= 1.");

  // Flatten to std::vector<double> for thread-safety
  std::vector<double> X(size_t(n) * p);
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < n; ++i)
      X[i + size_t(j) * n] = x_in(i, j);
  std::vector<double> Y(y_in.begin(), y_in.end());

  int ms = minspan_in > 0 ? minspan_in : auto_minspan(p, n);
  int es = endspan_in > 0 ? endspan_in : auto_endspan(p);

  // Cap nk at n - 1 for safety
  int nk_cap = std::min(nk, n - 1);
  if (nk_cap < 3) nk_cap = std::min(3, n - 1);

  // Allocate basis: B(n x nk_cap), dirs(nk_cap x p), cuts(nk_cap x p)
  std::vector<double> B(size_t(n) * nk_cap, 0.0);
  std::vector<int> dirs(size_t(nk_cap) * p, 0);
  std::vector<double> cuts(size_t(nk_cap) * p, 0.0);

  // Term 0: intercept
  for (int i = 0; i < n; ++i) B[i + 0 * size_t(n)] = 1.0;
  int M = 1;

  // Initial fit (intercept only)
  std::vector<double> beta, fitted, resid;
  double rss, rss0;
  recompute_residual(B.data(), n, M, Y.data(), beta, fitted, resid, rss);
  rss0 = rss;
  if (rss0 <= 0.0) rss0 = 1.0; // pathological all-zero y

  // Thread count is set from the R-side wrapper via RcppParallel::setThreadOptions.
  // Here we just observe the requested count to choose the serial vs parallel path.
  (void)nthreads;

  // --- Forward pass ---
  while (M < nk_cap - 1) {
    // Build candidate (parent, var) pairs
    std::vector<std::pair<int,int>> pairs;
    std::vector<int> parent_global_idx; // for each local parent index, its M-index
    std::vector<int> parent_locals;     // dedup locals
    // We need to iterate over current M parent terms; for each parent, exclude
    // variables already used in that parent (so depth doesn't exceed `degree`).
    pairs.reserve(size_t(M) * p);
    int parent_local_counter = 0;
    std::vector<int> parent_local_for_global(M, -1);
    // Pre-stack parent columns
    std::vector<double> B_cols_flat;
    B_cols_flat.reserve(size_t(M) * n);
    for (int parent_idx = 0; parent_idx < M; ++parent_idx) {
      int* dr = &dirs[size_t(parent_idx) * p];
      int td = term_degree(dr, p);
      if (td >= degree) continue;
      // candidate vars: those not already used in this parent term
      bool any_var = false;
      int local = parent_local_counter;
      for (int j = 0; j < p; ++j) {
        if (term_uses_var(dr, j)) continue;
        pairs.emplace_back(local, j);
        any_var = true;
      }
      if (any_var) {
        // commit this parent's column
        B_cols_flat.insert(B_cols_flat.end(),
                           B.begin() + size_t(parent_idx) * n,
                           B.begin() + size_t(parent_idx + 1) * n);
        parent_global_idx.push_back(parent_idx);
        ++parent_local_counter;
      }
    }
    if (pairs.empty()) break;

    // Build orthonormal Q of current basis B[:, 0:M] (read-only by workers)
    std::vector<double> Q(size_t(n) * M);
    int Mq = ares::build_Q(B.data(), n, M, Q.data());

    std::vector<Candidate> out(pairs.size());
    ForwardWorker w(X.data(), B_cols_flat.data(), resid.data(),
                    Q.data(), Mq,
                    pairs, parent_global_idx, n, p, ms, es, out);
    if (nthreads <= 1) {
      // serial
      w(0, pairs.size());
    } else {
      RcppParallel::parallelFor(0, pairs.size(), w);
    }
    // reduce
    Candidate best;
    for (auto& c : out) if (better(c, best)) best = c;
    if (best.rss_red <= 0.0 || best.parent < 0) break;

    // Add the hinge pair (-1 first then +1) — earth's order is "-" then "+"
    int parent = best.parent;
    int var = best.var;
    double cut = best.cut;
    int* dr_p = &dirs[size_t(parent) * p];
    double* cr_p = &cuts[size_t(parent) * p];
    int added = 0;
    for (int s : {1, -1}) {
      if (M >= nk_cap) break;
      int* dr_new = &dirs[size_t(M) * p];
      double* cr_new = &cuts[size_t(M) * p];
      // copy parent's dirs/cuts rows
      for (int j = 0; j < p; ++j) { dr_new[j] = dr_p[j]; cr_new[j] = cr_p[j]; }
      dr_new[var] = s;
      cr_new[var] = cut;
      // compute new column = parent_col * max(0, s*(x_var - cut))
      double* Bcol = &B[size_t(M) * n];
      const double* parent_col = &B[size_t(parent) * n];
      const double* xj = &X[size_t(var) * n];
      for (int i = 0; i < n; ++i) {
        double v = (s == 1) ? (xj[i] - cut) : (cut - xj[i]);
        Bcol[i] = (v > 0.0) ? parent_col[i] * v : 0.0;
      }
      ++M;
      ++added;
    }
    if (added == 0) break;
    // Recompute rss with new basis
    double new_rss;
    recompute_residual(B.data(), n, M, Y.data(), beta, fitted, resid, new_rss);
    double rel_imp = (rss - new_rss) / rss0;
    if (trace > 0) {
      Rcpp::Rcout << "forward step: M=" << M << " rss=" << new_rss
                  << " rel_imp=" << rel_imp << " (parent=" << parent
                  << ", var=" << var << ", cut=" << cut << ")\n";
    }
    if (rel_imp < thresh) {
      rss = new_rss;
      break;
    }
    rss = new_rss;
  }

  // --- Backward pass (pmethod == 0) ---
  std::vector<int> selected_idx;
  std::vector<double> rss_per_subset(M, NA_REAL);
  std::vector<double> gcv_per_subset(M, NA_REAL);
  double final_rss = rss, final_gcv;
  std::vector<double> final_beta;

  auto compute_gcv = [&](double rss_val, int n_terms) {
    double C = double(n_terms) + penalty * (double(n_terms) - 1.0) / 2.0;
    double denom = 1.0 - C / double(n);
    if (denom <= 0.0) return std::numeric_limits<double>::infinity();
    return rss_val / (double(n) * denom * denom);
  };

  if (pmethod == 0 && M > 1) {
    // Backward subset-selection: at each step drop the term whose removal causes
    // the smallest RSS increase, recompute RSS via OLS on retained columns.
    std::vector<int> cur(M);
    std::iota(cur.begin(), cur.end(), 0);
    // Initial fit (all M)
    double cur_rss = ols_qr(B.data(), n, M, Y.data(), beta, fitted);
    double cur_gcv = compute_gcv(cur_rss, M);
    rss_per_subset[M - 1] = cur_rss;
    gcv_per_subset[M - 1] = cur_gcv;
    int best_size = M;
    double best_gcv = cur_gcv;
    std::vector<int> best_set = cur;
    std::vector<double> best_beta = beta;

    while ((int)cur.size() > 1) {
      // Try removing each non-intercept term
      double best_drop_rss = std::numeric_limits<double>::infinity();
      int best_drop_idx = -1;
      std::vector<double> best_drop_beta;
      for (int ki = 1; ki < (int)cur.size(); ++ki) {
        std::vector<int> cand(cur);
        cand.erase(cand.begin() + ki);
        // build candidate B
        std::vector<double> Bc(size_t(n) * cand.size());
        for (int j = 0; j < (int)cand.size(); ++j)
          for (int i = 0; i < n; ++i)
            Bc[i + size_t(j) * n] = B[i + size_t(cand[j]) * n];
        std::vector<double> bb, ff;
        double cand_rss = ols_qr(Bc.data(), n, (int)cand.size(), Y.data(), bb, ff);
        if (cand_rss < best_drop_rss) {
          best_drop_rss = cand_rss;
          best_drop_idx = ki;
          best_drop_beta = bb;
        }
      }
      if (best_drop_idx < 0) break;
      cur.erase(cur.begin() + best_drop_idx);
      cur_rss = best_drop_rss;
      cur_gcv = compute_gcv(cur_rss, (int)cur.size());
      rss_per_subset[cur.size() - 1] = cur_rss;
      gcv_per_subset[cur.size() - 1] = cur_gcv;
      if (cur_gcv < best_gcv) {
        best_gcv = cur_gcv;
        best_size = (int)cur.size();
        best_set = cur;
        best_beta = best_drop_beta;
      }
    }
    selected_idx = best_set;
    final_beta = best_beta;
    final_rss = rss_per_subset[best_size - 1];
    final_gcv = best_gcv;
  } else {
    // pmethod == "none": keep all M terms
    selected_idx.resize(M);
    std::iota(selected_idx.begin(), selected_idx.end(), 0);
    double cur_rss = ols_qr(B.data(), n, M, Y.data(), beta, fitted);
    final_beta = beta;
    final_rss = cur_rss;
    final_gcv = compute_gcv(cur_rss, M);
    rss_per_subset[M - 1] = cur_rss;
    gcv_per_subset[M - 1] = final_gcv;
  }

  // Build outputs
  int L = (int)selected_idx.size();
  NumericMatrix bx(n, L);
  for (int j = 0; j < L; ++j)
    for (int i = 0; i < n; ++i)
      bx(i, j) = B[i + size_t(selected_idx[j]) * n];
  NumericVector coefs(final_beta.begin(), final_beta.end());

  // dirs and cuts: keep ALL forward-pass rows (size M x p)
  IntegerMatrix dirs_out(M, p);
  NumericMatrix cuts_out(M, p);
  for (int t = 0; t < M; ++t)
    for (int j = 0; j < p; ++j) {
      dirs_out(t, j) = dirs[size_t(t) * p + j];
      cuts_out(t, j) = cuts[size_t(t) * p + j];
    }

  IntegerVector sel(L);
  for (int j = 0; j < L; ++j) sel[j] = selected_idx[j] + 1; // 1-indexed for R

  NumericVector rss_ps(rss_per_subset.begin(), rss_per_subset.end());
  NumericVector gcv_ps(gcv_per_subset.begin(), gcv_per_subset.end());

  return List::create(
    _["coefficients"]   = coefs,
    _["bx"]             = bx,
    _["dirs"]           = dirs_out,
    _["cuts"]           = cuts_out,
    _["selected.terms"] = sel,
    _["rss"]            = final_rss,
    _["gcv"]            = final_gcv,
    _["rss.per.subset"] = rss_ps,
    _["gcv.per.subset"] = gcv_ps,
    _["nk"]             = nk_cap,
    _["thresh"]         = thresh,
    _["penalty"]        = penalty,
    _["minspan"]        = ms,
    _["endspan"]        = es,
    _["degree"]         = degree,
    _["nthreads"]       = nthreads
  );
}

// [[Rcpp::export]]
NumericMatrix mars_basis_cpp(const NumericMatrix& xnew,
                             const IntegerMatrix& dirs,
                             const NumericMatrix& cuts,
                             const IntegerVector& selected) {
  using namespace ares;
  int n = xnew.nrow();
  int p = xnew.ncol();
  int L = selected.size();
  if (dirs.nrow() != cuts.nrow() || dirs.ncol() != cuts.ncol())
    stop("ares: dirs and cuts must have the same dimensions.");
  if (dirs.ncol() != p)
    stop("ares: ncol(xnew) does not match ncol(dirs)/ncol(cuts).");
  std::vector<double> X(size_t(n) * p);
  for (int j = 0; j < p; ++j)
    for (int i = 0; i < n; ++i)
      X[i + size_t(j) * n] = xnew(i, j);
  std::vector<int> dr(size_t(dirs.nrow()) * p);
  std::vector<double> ct(size_t(cuts.nrow()) * p);
  for (int t = 0; t < dirs.nrow(); ++t)
    for (int j = 0; j < p; ++j) {
      dr[size_t(t) * p + j] = dirs(t, j);
      ct[size_t(t) * p + j] = cuts(t, j);
    }
  NumericMatrix bx(n, L);
  std::vector<double> col(n);
  for (int k = 0; k < L; ++k) {
    int term_id = selected[k] - 1;
    if (term_id < 0 || term_id >= dirs.nrow())
      stop("ares: selected.terms index out of range.");
    build_term_column(X.data(), n, p,
                      &dr[size_t(term_id) * p], &ct[size_t(term_id) * p],
                      col.data());
    for (int i = 0; i < n; ++i) bx(i, k) = col[i];
  }
  return bx;
}
