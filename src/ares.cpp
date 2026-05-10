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

// AVX2 intrinsics for the KnotScanner q-inner loops. Guarded so the
// package still compiles on architectures without AVX2; the scalar
// fallback path stays intact.
#if defined(__AVX2__)
#include <immintrin.h>
#define ARES_HAVE_AVX2 1
#else
#define ARES_HAVE_AVX2 0
#endif

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
  bool is_boundary = false;     // best knot was at the leftmost or rightmost
                                // eligible position — Auto.linpreds substitutes
                                // a linear basis here.
};

inline bool better(const Candidate& a, const Candidate& b) {
  // Tie-break: largest rss_red wins; on tie, smallest var, then smallest cut,
  // then smallest parent. Empirically aligns with earth on Friedman-1
  // benchmark seeds; deviates from Friedman 1991 paper's literal first-found
  // ordering, which would put parent ahead of var, but earth's behaviour wins
  // here because earth is the parity target.
  if (a.rss_red != b.rss_red) return a.rss_red > b.rss_red;
  if (a.var     != b.var    ) return a.var     < b.var;
  if (a.cut     != b.cut    ) return a.cut     < b.cut;
  return a.parent < b.parent;
}

// For a single (parent, var) pair, scan all eligible knots and compute the best
// RSS reduction from adding the *pair* of hinges (h+, h-) multiplied by the
// parent column.
//
// Friedman fast-LS scoring (1991, §3.5): instead of recomputing
// (Q' h+, Q' h-, ||h+||^2, ||h-||^2, h+ . r, h- . r) at O(n·Mq) per knot, we
// maintain prefix sums over the eligible rows sorted ascending by x_j and
// derive each per-knot quantity as a closed-form combination of those sums in
// O(Mq) per knot. Total cost per (parent, var) pair drops from O(K·n·Mq) to
// O((n + K)·Mq), where K is the number of candidate knots.
//
// Identities (with low-side L = rows with x_i < t, high-side H = rows with x_i > t):
//   h+(i)·r_i summed over H = (sum_{H} p_i x_i r_i) - t · (sum_{H} p_i r_i)
//   h-(i)·r_i summed over L = t · (sum_{L} p_i r_i) - (sum_{L} p_i x_i r_i)
//   ||h+||^2 = sum_{H} p_i^2 x_i^2 - 2t · sum_{H} p_i^2 x_i + t^2 · sum_{H} p_i^2
//   ||h-||^2 = t^2 · sum_{L} p_i^2 - 2t · sum_{L} p_i^2 x_i + sum_{L} p_i^2 x_i^2
//   (Q' h+)_q = sum_{H} Q_{i,q} p_i x_i - t · sum_{H} Q_{i,q} p_i
//   (Q' h-)_q = t · sum_{L} Q_{i,q} p_i - sum_{L} Q_{i,q} p_i x_i
// All "sum_{H}" expressions are reconstructed as Total - sum_{L} - row_k.
struct KnotScanner {
  // inputs
  const double* parent_col;     // n
  const double* xj;             // n
  const double* resid;          // n  (current residual: y - B beta, already orthogonal to span(Q))
  const double* QT;             // n × QT_stride row-major-over-q layout: QT[r*QT_stride + q] is Q[r,q].
                                // Per-row Mq slice is contiguous within the QT_stride block →
                                // unit-stride q-loop, auto-SIMD.
  const int* var_sort;          // length-n permutation of rows sorted by xj (ascending).
                                // Built once per variable per forward step (amortises across all
                                // parents sharing this variable). Filter to parent's support is O(n).
  int n;
  int Mq;                       // number of valid Q-cols to scan (cols 0..Mq-1)
  int QT_stride;                // byte stride between row r and r+1 in QT (≥ Mq)
  int minspan_user;             // user override; <=0 means auto via Friedman eq. 43 with N_m = n_eli
  int endspan;
  int p_pred;                   // total #predictors (used in auto_minspan formula)
  // output
  Candidate best;

  void run(int parent_idx, int var_idx) {
    // Filter the precomputed `var_sort` (rows sorted by xj) down to those
    // also in this parent's support (parent_col[r] != 0). Linear pass; the
    // O(n_eli log n_eli) sort that this used to do per-pair is now amortised
    // across every parent that shares this variable.
    std::vector<int> idx;
    idx.reserve(n);
    for (int i = 0; i < n; ++i) {
      int r = var_sort[i];
      if (parent_col[r] != 0.0) idx.push_back(r);
    }
    int n_eli = (int)idx.size();
    if (n_eli < 2 * endspan + 3) return;

    // Use raw n in auto-minspan (earth's behavior) when minspan_user <= 0.
    // Friedman 1991 eq. 43 uses N_m (parent's nonzero count); earth does not,
    // and the parity target is earth, so we match earth.
    int minspan = (minspan_user > 0) ? minspan_user : auto_minspan(p_pred, n);

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

    // Boundary-knot values: leftmost = xj at knot_pos.front(); rightmost at back().
    // Used at end of scan to flag is_boundary for Auto.linpreds.
    double leftmost_t  = xj[idx[knot_pos.front()]];
    double rightmost_t = xj[idx[knot_pos.back()]];

    // ---- Precompute totals over the eligible (parent_col != 0) rows --------
    // T_*: scalar sums; T_Q*: per-Q-column sums (length Mq each).
    double T_pr = 0.0, T_pxr = 0.0;
    double T_pp = 0.0, T_ppx = 0.0, T_ppxx = 0.0;
    std::vector<double> T_Qp(Mq, 0.0), T_Qpx(Mq, 0.0);
    for (int i = 0; i < n_eli; ++i) {
      int r = idx[i];
      double p  = parent_col[r];
      double xv = xj[r];
      double rv = resid[r];
      double pp = p * p;
      T_pr   += p * rv;
      T_pxr  += p * xv * rv;
      T_pp   += pp;
      T_ppx  += pp * xv;
      T_ppxx += pp * xv * xv;
      const double* QTr = QT + size_t(r) * QT_stride;  // unit-stride Mq slice
      double pxv = p * xv;
      int q = 0;
#if ARES_HAVE_AVX2
      __m256d vp = _mm256_set1_pd(p);
      __m256d vpxv = _mm256_set1_pd(pxv);
      for (; q + 4 <= Mq; q += 4) {
        __m256d Qq = _mm256_loadu_pd(QTr + q);
        __m256d Tp = _mm256_loadu_pd(T_Qp.data() + q);
        __m256d Tpx = _mm256_loadu_pd(T_Qpx.data() + q);
        Tp  = _mm256_add_pd(Tp,  _mm256_mul_pd(Qq, vp));
        Tpx = _mm256_add_pd(Tpx, _mm256_mul_pd(Qq, vpxv));
        _mm256_storeu_pd(T_Qp.data() + q, Tp);
        _mm256_storeu_pd(T_Qpx.data() + q, Tpx);
      }
#endif
      for (; q < Mq; ++q) {
        double Qiq = QTr[q];
        T_Qp[q]  += Qiq * p;
        T_Qpx[q] += Qiq * pxv;
      }
    }

    // ---- Sweep ascending in x_j, scoring knots in order --------------------
    // S_* = prefix sums over [0..k-1] (the "L" / low-side support at knot k).
    Candidate local_best;
    double S_pr = 0.0, S_pxr = 0.0;
    double S_pp = 0.0, S_ppx = 0.0, S_ppxx = 0.0;
    std::vector<double> S_Qp(Mq, 0.0), S_Qpx(Mq, 0.0);
    std::vector<double> qcp(Mq), qcm(Mq);

    const double tiny = 1e-12;
    size_t kp_pos = 0;  // next index into knot_pos to score

    for (int k = 0; k < n_eli; ++k) {
      // If row k is a scoring knot, score using prefix sums *before* incorporating row k.
      if (kp_pos < knot_pos.size() && knot_pos[kp_pos] == k) {
        int rk      = idx[k];
        double pk   = parent_col[rk];
        double xk   = xj[rk];
        double rk_v = resid[rk];
        double t    = xk;
        double pk2  = pk * pk;

        // High-side sums = Total - low-side - row k.
        double H_pr   = T_pr   - S_pr   - pk * rk_v;
        double H_pxr  = T_pxr  - S_pxr  - pk * xk * rk_v;
        double H_pp   = T_pp   - S_pp   - pk2;
        double H_ppx  = T_ppx  - S_ppx  - pk2 * xk;
        double H_ppxx = T_ppxx - S_ppxx - pk2 * xk * xk;

        // h+ . r and h- . r from precomputed sums.
        double cp_dot_r = H_pxr - t * H_pr;
        double cm_dot_r = t * S_pr - S_pxr;
        // ||h+||^2 and ||h-||^2.
        double cp_norm2 = H_ppxx - 2.0 * t * H_ppx + t * t * H_pp;
        double cm_norm2 = t * t * S_pp - 2.0 * t * S_ppx + S_ppxx;

        if (cp_norm2 >= 1e-300 || cm_norm2 >= 1e-300) {
          // (Q' h+)_q = H_Qpx_q - t * H_Qp_q;  (Q' h-)_q = t * S_Qp_q - S_Qpx_q
          double cp_perp_norm2 = cp_norm2;
          double cm_perp_norm2 = cm_norm2;
          double cp_perp_dot_cm_perp = 0.0;  // h+·h- = 0 (disjoint supports)
          const double* QTrk = QT + size_t(rk) * QT_stride;  // contiguous Mq slice
          double pkxk = pk * xk;
          int q = 0;
#if ARES_HAVE_AVX2
          __m256d vpk   = _mm256_set1_pd(pk);
          __m256d vpkxk = _mm256_set1_pd(pkxk);
          __m256d vt    = _mm256_set1_pd(t);
          __m256d acc_pp = _mm256_setzero_pd();
          __m256d acc_mm = _mm256_setzero_pd();
          __m256d acc_pm = _mm256_setzero_pd();
          for (; q + 4 <= Mq; q += 4) {
            __m256d Qkq = _mm256_loadu_pd(QTrk + q);
            __m256d Tp  = _mm256_loadu_pd(T_Qp.data() + q);
            __m256d Tpx = _mm256_loadu_pd(T_Qpx.data() + q);
            __m256d Sp  = _mm256_loadu_pd(S_Qp.data() + q);
            __m256d Spx = _mm256_loadu_pd(S_Qpx.data() + q);
            __m256d HQp  = _mm256_sub_pd(_mm256_sub_pd(Tp, Sp),  _mm256_mul_pd(Qkq, vpk));
            __m256d HQpx = _mm256_sub_pd(_mm256_sub_pd(Tpx, Spx), _mm256_mul_pd(Qkq, vpkxk));
            __m256d dp = _mm256_sub_pd(HQpx, _mm256_mul_pd(vt, HQp));
            __m256d dm = _mm256_sub_pd(_mm256_mul_pd(vt, Sp), Spx);
            _mm256_storeu_pd(qcp.data() + q, dp);
            _mm256_storeu_pd(qcm.data() + q, dm);
            acc_pp = _mm256_add_pd(acc_pp, _mm256_mul_pd(dp, dp));
            acc_mm = _mm256_add_pd(acc_mm, _mm256_mul_pd(dm, dm));
            acc_pm = _mm256_add_pd(acc_pm, _mm256_mul_pd(dp, dm));
          }
          alignas(32) double rbuf[4];
          _mm256_store_pd(rbuf, acc_pp);
          cp_perp_norm2       -= rbuf[0] + rbuf[1] + rbuf[2] + rbuf[3];
          _mm256_store_pd(rbuf, acc_mm);
          cm_perp_norm2       -= rbuf[0] + rbuf[1] + rbuf[2] + rbuf[3];
          _mm256_store_pd(rbuf, acc_pm);
          cp_perp_dot_cm_perp -= rbuf[0] + rbuf[1] + rbuf[2] + rbuf[3];
#endif
          for (; q < Mq; ++q) {
            double Qkq = QTrk[q];
            double H_Qp_q  = T_Qp[q]  - S_Qp[q]  - Qkq * pk;
            double H_Qpx_q = T_Qpx[q] - S_Qpx[q] - Qkq * pkxk;
            double dp = H_Qpx_q - t * H_Qp_q;
            double dm = t * S_Qp[q] - S_Qpx[q];
            qcp[q] = dp;
            qcm[q] = dm;
            cp_perp_norm2       -= dp * dp;
            cm_perp_norm2       -= dm * dm;
            cp_perp_dot_cm_perp -= dp * dm;
          }

          // Joint reduction via 2-step Gram-Schmidt (same as before).
          double red_p = 0.0, red_m = 0.0;
          if (cp_perp_norm2 > tiny) {
            red_p = cp_dot_r * cp_dot_r / cp_perp_norm2;
            double alpha = cp_perp_dot_cm_perp / cp_perp_norm2;
            double e2_norm2 = cm_perp_norm2 - alpha * cp_perp_dot_cm_perp;
            double e2_dot_r = cm_dot_r - alpha * cp_dot_r;
            if (e2_norm2 > tiny) red_m = e2_dot_r * e2_dot_r / e2_norm2;
          } else if (cm_perp_norm2 > tiny) {
            red_m = cm_dot_r * cm_dot_r / cm_perp_norm2;
          }
          double red = red_p + red_m;
          // Tie-break: prefer larger t on equal RSS reduction (paper §3.9
          // descending-knot sweep, first-encountered wins). var_idx is fixed
          // within this scanner so the var clause is dead code; kept for
          // symmetry with `better()` if scanner is ever re-used cross-pair.
          if (red > local_best.rss_red ||
              (red == local_best.rss_red && t < local_best.cut)) {
            local_best.rss_red = red;
            local_best.parent  = parent_idx;
            local_best.var     = var_idx;
            local_best.cut     = t;
          }
        }
        ++kp_pos;
      }

      // Accumulate row k into low-side prefix sums for the NEXT knot.
      int r    = idx[k];
      double p = parent_col[r];
      double xv = xj[r];
      double rv = resid[r];
      double pp = p * p;
      S_pr   += p * rv;
      S_pxr  += p * xv * rv;
      S_pp   += pp;
      S_ppx  += pp * xv;
      S_ppxx += pp * xv * xv;
      const double* QTr2 = QT + size_t(r) * QT_stride;
      double pxv2 = p * xv;
      int q2 = 0;
#if ARES_HAVE_AVX2
      __m256d vp2 = _mm256_set1_pd(p);
      __m256d vpxv2 = _mm256_set1_pd(pxv2);
      for (; q2 + 4 <= Mq; q2 += 4) {
        __m256d Qq = _mm256_loadu_pd(QTr2 + q2);
        __m256d Sp = _mm256_loadu_pd(S_Qp.data() + q2);
        __m256d Spx = _mm256_loadu_pd(S_Qpx.data() + q2);
        Sp  = _mm256_add_pd(Sp,  _mm256_mul_pd(Qq, vp2));
        Spx = _mm256_add_pd(Spx, _mm256_mul_pd(Qq, vpxv2));
        _mm256_storeu_pd(S_Qp.data() + q2, Sp);
        _mm256_storeu_pd(S_Qpx.data() + q2, Spx);
      }
#endif
      for (; q2 < Mq; ++q2) {
        double Qiq = QTr2[q2];
        S_Qp[q2]  += Qiq * p;
        S_Qpx[q2] += Qiq * pxv2;
      }
    }
    if (local_best.parent >= 0) {
      local_best.is_boundary = (local_best.cut == leftmost_t ||
                                local_best.cut == rightmost_t);
    }
    if (better(local_best, best)) best = local_best;
  }
};

struct ForwardWorker : public RcppParallel::Worker {
  const double* x;             // n*p flat column-major
  const double* B_cols;         // n*M_active flat column-major (parent columns, in same order as parent_global_idx)
  const double* resid;          // n
  const double* QT;             // n × QT_stride, row-major-over-q (per-row Mq slice contiguous)
  int Mq;
  int QT_stride;
  const std::vector<std::pair<int,int>>* pairs;  // (parent_local_idx, var_idx)
  const std::vector<int>* parent_global_idx;
  const std::vector<int>* pair_endspans;          // per-pair endspan (Adjust.endspan applied)
  const int* var_sort_flat;     // p * n ints, [v*n + i] = rows sorted by x[:, v] ascending
  int n, p;
  int minspan_user;             // minspan_user <=0 ⇒ auto per parent
  std::vector<Candidate>* out;  // size = pairs->size()

  ForwardWorker(const double* x_, const double* B_, const double* r_,
                const double* QT_, int Mq_, int QT_stride_,
                const std::vector<std::pair<int,int>>& pairs_,
                const std::vector<int>& parent_global_,
                const std::vector<int>& pair_endspans_,
                const int* var_sort_flat_,
                int n_, int p_, int ms_,
                std::vector<Candidate>& out_)
    : x(x_), B_cols(B_), resid(r_), QT(QT_), Mq(Mq_), QT_stride(QT_stride_),
      pairs(&pairs_), parent_global_idx(&parent_global_),
      pair_endspans(&pair_endspans_),
      var_sort_flat(var_sort_flat_),
      n(n_), p(p_), minspan_user(ms_), out(&out_) {}

  void operator()(std::size_t begin, std::size_t end) override {
    for (std::size_t i = begin; i < end; ++i) {
      int parent_local = (*pairs)[i].first;
      int var_idx     = (*pairs)[i].second;
      int es_pair     = (*pair_endspans)[i];
      KnotScanner sc{
        B_cols + size_t(parent_local) * n,
        x      + size_t(var_idx)     * n,
        resid, QT,
        var_sort_flat + size_t(var_idx) * n,
        n, Mq, QT_stride, minspan_user, es_pair, p, Candidate{}
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

// ============================================================================
// QR-with-downdate machinery for the backward pass.
//
// Backward subset selection drops one term per step; previously every step
// rebuilt a full Householder QR for each of M-1 candidate removals at O(n·M²)
// per trial. Profiling at v0.1 showed this loop was 63 % of wall-clock at
// n=5000. Replace with: factor B once, then maintain R (M×M upper triangular)
// across removals via Givens-rotation column downdate at O(M²) per step;
// pick the best column to drop in closed form using
//   ΔRSS_j = β_j² / [(R'R)^{-1}]_{j,j}
// (the standard "increase from constrained-to-zero" formula). One full
// Householder QR up front: O(n·M²); subsequent M backward steps: O(M³) each
// for the diag-of-inverse + downdate; total O(n·M² + M⁴).
// ============================================================================

// In-place Householder QR of (B|y). On exit:
//   - B's upper triangle (rows 0..M-1, cols 0..M-1) holds R; lower part holds
//     Householder reflectors (we don't need them again).
//   - y[0..M-1] holds Q'y (the regression "RHS" Qty).
//   - rss = ||y||^2 - ||Qty||^2 = sum of squared residual on rows M..n-1.
// R_out is written column-major with stride R_stride (>= M); the M×M
// upper-tri block lives at R_out[i + j*R_stride] for i, j in [0, M).
static void householder_qr_R(double* B, int n, int M, double* y,
                             double* R_out, int R_stride,
                             double* Qty, double& rss) {
  for (int k = 0; k < M; ++k) {
    double norm_x = 0.0;
    for (int i = k; i < n; ++i) norm_x += B[i + size_t(k)*n] * B[i + size_t(k)*n];
    norm_x = std::sqrt(norm_x);
    if (norm_x == 0.0) continue;
    double sign = B[k + size_t(k)*n] >= 0 ? 1.0 : -1.0;
    double alpha = -sign * norm_x;
    B[k + size_t(k)*n] -= alpha;
    double v_norm_sq = 0.0;
    for (int i = k; i < n; ++i) v_norm_sq += B[i + size_t(k)*n] * B[i + size_t(k)*n];
    if (v_norm_sq == 0.0) continue;
    double tau = 2.0 / v_norm_sq;
    for (int j = k + 1; j < M; ++j) {
      double dot = 0.0;
      for (int i = k; i < n; ++i) dot += B[i + size_t(k)*n] * B[i + size_t(j)*n];
      dot *= tau;
      for (int i = k; i < n; ++i) B[i + size_t(j)*n] -= dot * B[i + size_t(k)*n];
    }
    double dot_y = 0.0;
    for (int i = k; i < n; ++i) dot_y += B[i + size_t(k)*n] * y[i];
    dot_y *= tau;
    for (int i = k; i < n; ++i) y[i] -= dot_y * B[i + size_t(k)*n];
    B[k + size_t(k)*n] = alpha;
  }
  for (int j = 0; j < M; ++j) {
    for (int i = 0; i <= j; ++i) R_out[i + size_t(j)*R_stride] = B[i + size_t(j)*n];
    for (int i = j + 1; i < M; ++i) R_out[i + size_t(j)*R_stride] = 0.0;
  }
  for (int i = 0; i < M; ++i) Qty[i] = y[i];
  rss = 0.0;
  for (int i = M; i < n; ++i) rss += y[i] * y[i];
}

// Solve R x = b where R is upper triangular (size M_active in an M_alloc-stride
// column-major buffer). Pseudo-inverts near-zero diagonal pivots.
static void solve_R_back(const double* R, int M_alloc, int M_active,
                         const double* b, double* x) {
  const double tol_eps = 1e3 * std::numeric_limits<double>::epsilon();
  for (int k = M_active - 1; k >= 0; --k) {
    double sum = b[k];
    for (int j = k + 1; j < M_active; ++j) sum -= R[k + size_t(j)*M_alloc] * x[j];
    double dkk = R[k + size_t(k)*M_alloc];
    x[k] = (std::fabs(dkk) < tol_eps * (1.0 + std::fabs(sum))) ? 0.0 : sum / dkk;
  }
}

// diag[j] = [(R'R)^{-1}]_{j,j} = sum_{k>=j} W[j,k]^2 where W = R^{-1}
// (upper triangular). O(M³) per call.
static void chol_diag_inverse(const double* R, int M_alloc, int M_active,
                              double* diag) {
  const double tol_eps = 1e3 * std::numeric_limits<double>::epsilon();
  std::vector<double> W(size_t(M_active) * M_active, 0.0);
  for (int j = 0; j < M_active; ++j) {
    double dkk = R[j + size_t(j)*M_alloc];
    W[j + size_t(j)*M_active] = (std::fabs(dkk) < tol_eps) ? 0.0 : 1.0 / dkk;
    for (int k = j - 1; k >= 0; --k) {
      double sum = 0.0;
      for (int i = k + 1; i <= j; ++i)
        sum -= R[k + size_t(i)*M_alloc] * W[i + size_t(j)*M_active];
      double rkk = R[k + size_t(k)*M_alloc];
      W[k + size_t(j)*M_active] = (std::fabs(rkk) < tol_eps) ? 0.0 : sum / rkk;
    }
  }
  for (int j = 0; j < M_active; ++j) {
    double s = 0.0;
    for (int k = j; k < M_active; ++k) {
      double w = W[j + size_t(k)*M_active];
      s += w * w;
    }
    diag[j] = s;
  }
}

// Remove column `k` from the M_active-active R (col-major in M_alloc-stride
// storage), and apply the resulting row-Givens rotations to b (length M_active).
// On exit, R is (M_active-1)×(M_active-1) upper triangular in cols 0..M_active-2
// and rows 0..M_active-2. b[0..M_active-2] is the new Q'y; b[M_active-1] is
// the orthogonal residue, whose square is exactly the RSS increase from
// dropping column k. Caller must decrement M_active.
static double qr_downdate_col(double* R, int M_alloc, int M_active, int k,
                              double* b) {
  // 1) Shift columns k+1..M_active-1 left by one (overwriting column k).
  for (int j = k; j < M_active - 1; ++j) {
    for (int i = 0; i <= j + 1; ++i)
      R[i + size_t(j)*M_alloc] = R[i + size_t(j + 1)*M_alloc];
  }
  // 2) Re-triangulate via Givens rotations on rows (c, c+1), c = k..M_active-2.
  for (int c = k; c < M_active - 1; ++c) {
    double a = R[c     + size_t(c)*M_alloc];
    double v = R[c + 1 + size_t(c)*M_alloc];
    if (v == 0.0) continue;
    double r = std::hypot(a, v);
    if (r == 0.0) continue;
    double cs = a / r;
    double sn = v / r;
    R[c     + size_t(c)*M_alloc] = r;
    R[c + 1 + size_t(c)*M_alloc] = 0.0;
    for (int j = c + 1; j < M_active - 1; ++j) {
      double rcj  = R[c     + size_t(j)*M_alloc];
      double rc1j = R[c + 1 + size_t(j)*M_alloc];
      R[c     + size_t(j)*M_alloc] =  cs * rcj + sn * rc1j;
      R[c + 1 + size_t(j)*M_alloc] = -sn * rcj + cs * rc1j;
    }
    double bc  = b[c];
    double bc1 = b[c + 1];
    b[c]     =  cs * bc + sn * bc1;
    b[c + 1] = -sn * bc + cs * bc1;
  }
  // RSS increase = (b[M_active-1])^2 (orthogonal residue from the eliminated row)
  return b[M_active - 1] * b[M_active - 1];
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
                  int minspan_in, int endspan_in,
                  int adjust_endspan, int auto_linpreds,
                  int nprune,
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

  // ms <= 0 means "auto per parent"; KnotScanner derives auto_minspan from
  // the actual N_m at scan time. ms > 0 forces a fixed user value.
  int ms = (minspan_in > 0) ? minspan_in : 0;
  int es = (endspan_in > 0) ? endspan_in : auto_endspan(p);
  // For reporting in the returned list, expose the auto-derived value at the
  // top-level n (the value used for full-support parents like the intercept).
  int ms_reported = (minspan_in > 0) ? minspan_in : auto_minspan(p, n);

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

  // Incremental QR maintenance (replaces per-step build_Q + ols_qr).
  //   Q1: n × nk_cap, orthonormal columns of B[:, 0:M_active] (col-major).
  //   QT: n × nk_cap, row-major-over-q transpose of Q1 (per-row Mq slice
  //       contiguous → unit-stride q-loop in KnotScanner, auto-vectorisable).
  //   R:  nk_cap × nk_cap upper triangular, the QR factor of B[:, 0:M_active].
  //   Qty: length nk_cap, Qty[q] = Q1[:, q]' · y (the OLS RHS in QR form).
  //   resid = y - Q1[:, 0:M_active] · Qty[0:M_active]. Per added column this
  //   updates as resid -= Q1[:, M] * Qty[M]; an O(n) update vs the old
  //   O(n·M²) recompute_residual call.
  std::vector<double> Q1(size_t(n) * nk_cap, 0.0);
  std::vector<double> QT(size_t(n) * nk_cap, 0.0);
  std::vector<double> R(size_t(nk_cap) * nk_cap, 0.0);
  std::vector<double> Qty(nk_cap, 0.0);
  std::vector<double> resid(n);

  // Initialise QR for M=1 (intercept). Column is all-ones; ||1|| = sqrt(n);
  // Q1[:, 0] = 1/sqrt(n); R[0,0] = sqrt(n); Qty[0] = sum(y)/sqrt(n).
  {
    double sn = std::sqrt(double(n));
    double inv_sn = 1.0 / sn;
    double sum_y = 0.0;
    for (int i = 0; i < n; ++i) sum_y += Y[i];
    R[0] = sn;
    Qty[0] = sum_y * inv_sn;
    for (int i = 0; i < n; ++i) {
      Q1[i] = inv_sn;
      QT[size_t(i) * nk_cap + 0] = inv_sn;
    }
    double fitted0 = Qty[0] * inv_sn;
    for (int i = 0; i < n; ++i) resid[i] = Y[i] - fitted0;
  }
  double sum_y2 = 0.0;
  for (int i = 0; i < n; ++i) sum_y2 += Y[i] * Y[i];
  double rss = sum_y2 - Qty[0] * Qty[0];
  double rss0 = rss;
  if (rss0 <= 0.0) rss0 = 1.0; // pathological all-zero y

  // Lambda: append one column c (length n) to the maintained QR. Returns
  // true if the column was accepted (full-rank) and Q1/R/Qty/resid/rss
  // were updated; false if c is in span(existing basis) within tolerance,
  // in which case the maintained state is left unchanged. The caller must
  // then NOT increment M (so the rejected column never enters B).
  // Keeping rank-deficient columns out of B means the maintained R has no
  // zero diagonals — the Givens downdate path in backward works directly
  // off this R, removing the need for the per-step Householder commit.
  const double tol_eps = 1e3 * std::numeric_limits<double>::epsilon();
  auto qr_append_col = [&](const double* c) -> bool {
    int M_now = M;
    if (M_now >= nk_cap) return false;
    std::vector<double> u(M_now, 0.0);
    for (int q = 0; q < M_now; ++q) {
      const double* Qcol = Q1.data() + size_t(q) * n;
      double s = 0.0;
      for (int i = 0; i < n; ++i) s += Qcol[i] * c[i];
      u[q] = s;
    }
    double c_norm2 = 0.0;
    for (int i = 0; i < n; ++i) c_norm2 += c[i] * c[i];
    double u_norm2 = 0.0;
    for (int q = 0; q < M_now; ++q) u_norm2 += u[q] * u[q];
    double d2 = c_norm2 - u_norm2;
    double col_scale = 0.0;
    for (int i = 0; i < n; ++i) {
      double a = std::fabs(c[i]);
      if (a > col_scale) col_scale = a;
    }
    if (d2 < tol_eps * tol_eps * (col_scale + 1.0) * (col_scale + 1.0)) {
      return false;  // rank-deficient: reject; caller will skip this slot
    }
    double d = std::sqrt(d2);
    double inv_d = 1.0 / d;
    double* Qcol_new = Q1.data() + size_t(M_now) * n;
    for (int i = 0; i < n; ++i) {
      double cp = c[i];
      for (int q = 0; q < M_now; ++q) cp -= Q1[size_t(q) * n + i] * u[q];
      Qcol_new[i] = cp * inv_d;
      QT[size_t(i) * nk_cap + M_now] = Qcol_new[i];
    }
    for (int q = 0; q < M_now; ++q) R[q + size_t(M_now) * nk_cap] = u[q];
    R[M_now + size_t(M_now) * nk_cap] = d;
    double cy = 0.0;
    for (int i = 0; i < n; ++i) cy += c[i] * Y[i];
    double uQ = 0.0;
    for (int q = 0; q < M_now; ++q) uQ += u[q] * Qty[q];
    double qty_new = (cy - uQ) * inv_d;
    Qty[M_now] = qty_new;
    for (int i = 0; i < n; ++i) resid[i] -= Qcol_new[i] * qty_new;
    rss -= qty_new * qty_new;
    return true;
  };

  // Thread count is set from the R-side wrapper via RcppParallel::setThreadOptions.
  (void)nthreads;

  // --- Forward pass ---
  while (M < nk_cap - 1) {
    // Build candidate (parent, var) pairs
    std::vector<std::pair<int,int>> pairs;
    std::vector<int> pair_endspans;     // per-pair endspan (Adjust.endspan applied)
    std::vector<int> parent_global_idx; // for each local parent index, its M-index
    // We need to iterate over current M parent terms; for each parent, exclude
    // variables already used in that parent (so depth doesn't exceed `degree`).
    pairs.reserve(size_t(M) * p);
    pair_endspans.reserve(size_t(M) * p);
    int parent_local_counter = 0;
    // Pre-stack parent columns
    std::vector<double> B_cols_flat;
    B_cols_flat.reserve(size_t(M) * n);
    for (int parent_idx = 0; parent_idx < M; ++parent_idx) {
      int* dr = &dirs[size_t(parent_idx) * p];
      int td = term_degree(dr, p);
      if (td >= degree) continue;
      // Adjust.endspan: when adding a hinge would create a deg>=2 term, scale
      // up the endspan to keep candidate knots away from variable boundaries.
      // Earth default `Adjust.endspan = 2`. Setting it to 1 disables.
      int eff_es = (td >= 1) ? adjust_endspan * es : es;
      // candidate vars: those not already used in this parent term
      bool any_var = false;
      int local = parent_local_counter;
      for (int j = 0; j < p; ++j) {
        if (term_uses_var(dr, j)) continue;
        pairs.emplace_back(local, j);
        pair_endspans.push_back(eff_es);
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

    // QT, resid, rss are maintained incrementally across forward steps via
    // qr_append_col below (replacing the per-step build_Q + ols_qr cost).
    // Mq == M because we always append a column on each emission (zeroing
    // it instead of dropping when rank-deficient — keeps M_q in sync).
    int Mq = M;

    // Per-var sorted index: rows sorted ascending by x[:, v] (stable on row
    // index for ties). Built once per forward step and shared by all
    // (parent, var) pairs that pick this v. KnotScanner filters this in
    // O(n) instead of running its own O(n log n) sort per pair.
    // The p sorts are independent → parallelise.
    std::vector<int> var_sort_flat(size_t(p) * n);
    {
      struct VarSortWorker : public RcppParallel::Worker {
        const double* X_ptr;
        int n;
        std::vector<int>* dst;
        VarSortWorker(const double* X_, int n_, std::vector<int>& d)
          : X_ptr(X_), n(n_), dst(&d) {}
        void operator()(std::size_t begin, std::size_t end) override {
          for (std::size_t v = begin; v < end; ++v) {
            int* d = dst->data() + size_t(v) * n;
            const double* xv = X_ptr + size_t(v) * n;
            for (int i = 0; i < n; ++i) d[i] = i;
            std::stable_sort(d, d + n,
                             [xv](int a, int b){
                               if (xv[a] != xv[b]) return xv[a] < xv[b];
                               return a < b;
                             });
          }
        }
      };
      VarSortWorker vw(X.data(), n, var_sort_flat);
      if (nthreads <= 1 || p < 2) {
        vw(0, (size_t)p);
      } else {
        RcppParallel::parallelFor(0, (size_t)p, vw);
      }
    }

    std::vector<Candidate> out(pairs.size());
    ForwardWorker w(X.data(), B_cols_flat.data(), resid.data(),
                    QT.data(), Mq, nk_cap,
                    pairs, parent_global_idx, pair_endspans,
                    var_sort_flat.data(),
                    n, p, ms, out);
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

    // Emission. Default is the hinge pair (s=+1 then s=-1). When
    // Auto.linpreds is enabled and the best knot was at the variable's
    // boundary (leftmost or rightmost eligible position in the parent's
    // support), substitute a single linear term (dirs=2) instead — the
    // hinge would degenerate to a half-line that reproduces a linear
    // function over the parent's support.
    int parent = best.parent;
    int var = best.var;
    double cut = best.cut;
    int* dr_p = &dirs[size_t(parent) * p];
    double* cr_p = &cuts[size_t(parent) * p];
    int added = 0;
    double rss_before = rss;
    if (auto_linpreds && best.is_boundary) {
      if (M < nk_cap) {
        double* Bcol = &B[size_t(M) * n];
        const double* parent_col = &B[size_t(parent) * n];
        const double* xj = &X[size_t(var) * n];
        for (int i = 0; i < n; ++i) Bcol[i] = parent_col[i] * xj[i];
        if (qr_append_col(Bcol)) {
          int* dr_new = &dirs[size_t(M) * p];
          double* cr_new = &cuts[size_t(M) * p];
          for (int j = 0; j < p; ++j) { dr_new[j] = dr_p[j]; cr_new[j] = cr_p[j]; }
          dr_new[var] = 2;
          cr_new[var] = 0.0;
          ++M;
          ++added;
        }
      }
    } else {
      for (int s : {1, -1}) {
        if (M >= nk_cap) break;
        double* Bcol = &B[size_t(M) * n];
        const double* parent_col = &B[size_t(parent) * n];
        const double* xj = &X[size_t(var) * n];
        for (int i = 0; i < n; ++i) {
          double v = (s == 1) ? (xj[i] - cut) : (cut - xj[i]);
          Bcol[i] = (v > 0.0) ? parent_col[i] * v : 0.0;
        }
        if (qr_append_col(Bcol)) {
          int* dr_new = &dirs[size_t(M) * p];
          double* cr_new = &cuts[size_t(M) * p];
          for (int j = 0; j < p; ++j) { dr_new[j] = dr_p[j]; cr_new[j] = cr_p[j]; }
          dr_new[var] = s;
          cr_new[var] = cut;
          ++M;
          ++added;
        }
      }
    }
    if (added == 0) break;
    // Class (e) numerical guard: with the QR maintained incrementally, rss
    // is the exact OLS minimum so a blowup would only come from a
    // catastrophically ill-conditioned new column. The qr_append_col path
    // already zeroes the column on rank deficiency (rss unchanged), so the
    // guard is now mostly a backstop.
    if (!std::isfinite(rss) || rss > 1e6 * rss0 + 1.0) {
      // Roll back: M -= added; rss not easily restorable, but we're stopping
      // anyway. Caller sees the rolled-back M.
      M -= added;
      if (trace > 0) {
        Rcpp::Rcout << "forward step: REJECTED ill-conditioned pair "
                    << "(parent=" << parent << ", var=" << var
                    << ", cut=" << cut << ") — rss=" << rss
                    << "; stopping forward pass at M=" << M << "\n";
      }
      break;
    }
    double rel_imp = (rss_before - rss) / rss0;
    if (trace > 0) {
      Rcpp::Rcout << "forward step: M=" << M << " rss=" << rss
                  << " rel_imp=" << rel_imp << " (parent=" << parent
                  << ", var=" << var << ", cut=" << cut << ")\n";
    }
    if (rel_imp < thresh) break;
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
    // Backward subset-selection via QR-with-downdate.
    //   1. One initial Householder QR of B[:, 0:M] gives R (M_alloc×M_alloc),
    //      Qty (length M_alloc) and rss_full = ||y||^2 - ||Qty||^2.
    //   2. Per backward step:
    //      - For each candidate drop ki ∈ [1, M_active-1], run a Givens
    //        downdate on a *copy* of (R, Qty); the orthogonal-residue squared
    //        b[M_active-1]^2 of the rotated b is the exact RSS increase.
    //        Trials are independent and parallelisable.
    //      - Pick smallest RSS = cur_rss + rss_inc.
    //      - Re-do the chosen downdate on the actual (R, Qty) to commit.
    //      - cur.erase(ki*); M_active--.
    // Per-step cost: O(M^2) per trial × M-1 trials = O(M^3); plus O(M^2)
    // commit. Total backward: O(M^4). At n=5000, M=21 this is ~200K ops vs
    // ~10^9 ops for the per-trial Householder rebuild approach — ~5000x.
    std::vector<int> cur(M);
    std::iota(cur.begin(), cur.end(), 0);

    // Initial Householder QR for backward. Conceptually we could reuse the
    // forward-maintained (R, Qty) directly, but empirically the
    // outer-state has slightly different rounding and forward and
    // backward end up disagreeing on which knot to drop on tightly tied
    // designs. Cost was 0.8% in v0.4 — keep the fresh init for safety.
    int M_alloc = M;
    std::vector<double> R_b(size_t(M_alloc) * M_alloc, 0.0);
    std::vector<double> Qty_b(M_alloc, 0.0);
    double cur_rss;
    {
      std::vector<double> Bw(B.begin(), B.begin() + size_t(M) * n);
      std::vector<double> yw(Y.begin(), Y.end());
      ares::householder_qr_R(Bw.data(), n, M, yw.data(),
                             R_b.data(), M_alloc, Qty_b.data(), cur_rss);
    }
    int M_active = M;

    std::vector<double> beta_active(M_active, 0.0);
    ares::solve_R_back(R_b.data(), M_alloc, M_active, Qty_b.data(), beta_active.data());

    double cur_gcv = compute_gcv(cur_rss, M_active);
    rss_per_subset[M_active - 1] = cur_rss;
    gcv_per_subset[M_active - 1] = cur_gcv;
    int best_size = M_active;
    double best_gcv = cur_gcv;
    std::vector<int> best_set = cur;
    std::vector<double> best_beta = beta_active;

    struct DowndateTrialWorker : public RcppParallel::Worker {
      const double* R_ptr;
      const double* Qty_ptr;
      int M_alloc;
      int M_active;
      std::vector<double>* rss_inc_out;  // length M_active - 1
      DowndateTrialWorker(const double* R_, const double* Q_, int Ma, int Mc,
                          std::vector<double>& out)
        : R_ptr(R_), Qty_ptr(Q_), M_alloc(Ma), M_active(Mc), rss_inc_out(&out) {}
      void operator()(std::size_t begin, std::size_t end) override {
        // Each trial works on its own copy of R and Qty (so ki may run in any order).
        std::vector<double> R_trial(size_t(M_alloc) * M_active);
        std::vector<double> Qty_trial(M_active);
        for (std::size_t t = begin; t < end; ++t) {
          int ki = (int)t + 1;
          if (ki >= M_active) continue;
          // Copy active sub-region of R (cols 0..M_active-1) and Qty.
          for (int j = 0; j < M_active; ++j)
            for (int i = 0; i <= j && i < M_active; ++i)
              R_trial[i + size_t(j) * M_alloc] = R_ptr[i + size_t(j) * M_alloc];
          for (int i = 0; i < M_active; ++i) Qty_trial[i] = Qty_ptr[i];
          double inc = ares::qr_downdate_col(R_trial.data(), M_alloc, M_active,
                                             ki, Qty_trial.data());
          (*rss_inc_out)[t] = inc;
        }
      }
    };

    while (M_active > 1) {
      std::vector<double> rss_inc(M_active - 1, std::numeric_limits<double>::infinity());
      DowndateTrialWorker w(R_b.data(), Qty_b.data(), M_alloc, M_active, rss_inc);
      if (nthreads <= 1 || M_active - 1 < 4) {
        w(0, (size_t)(M_active - 1));
      } else {
        RcppParallel::parallelFor(0, (size_t)(M_active - 1), w);
      }
      double best_inc = std::numeric_limits<double>::infinity();
      int best_t = -1;
      for (int t = 0; t < M_active - 1; ++t) {
        if (rss_inc[t] < best_inc) { best_inc = rss_inc[t]; best_t = t; }
      }
      if (best_t < 0) break;
      int drop_ki = best_t + 1;

      // Commit: redo the chosen Givens downdate on the actual (R_b, Qty_b).
      // Periodically (every 4 steps) refresh from a fresh Householder QR
      // of the surviving B columns to eat any cumulative rotation drift.
      // Drops the per-step O(n·M²) Householder cost on most steps.
      double commit_inc = ares::qr_downdate_col(R_b.data(), M_alloc, M_active,
                                                drop_ki, Qty_b.data());
      cur.erase(cur.begin() + drop_ki);
      cur_rss += commit_inc;
      --M_active;

      static const int REFRESH_EVERY = 4;
      bool need_refresh = (M_active >= 2) && ((M - M_active) % REFRESH_EVERY == 0);
      if (need_refresh) {
        std::vector<double> Bw_re(size_t(n) * M_active);
        for (int j = 0; j < M_active; ++j)
          for (int i = 0; i < n; ++i)
            Bw_re[i + size_t(j) * n] = B[i + size_t(cur[j]) * n];
        std::vector<double> yw_re(Y.begin(), Y.end());
        std::fill(R_b.begin(), R_b.end(), 0.0);
        std::fill(Qty_b.begin(), Qty_b.end(), 0.0);
        ares::householder_qr_R(Bw_re.data(), n, M_active, yw_re.data(),
                               R_b.data(), M_alloc, Qty_b.data(), cur_rss);
      }

      beta_active.assign(M_active, 0.0);
      ares::solve_R_back(R_b.data(), M_alloc, M_active, Qty_b.data(), beta_active.data());

      cur_gcv = compute_gcv(cur_rss, M_active);
      rss_per_subset[M_active - 1] = cur_rss;
      gcv_per_subset[M_active - 1] = cur_gcv;
      if (cur_gcv < best_gcv) {
        best_gcv = cur_gcv;
        best_size = M_active;
        best_set = cur;
        best_beta = beta_active;
      }
    }
    selected_idx = best_set;
    final_beta = best_beta;
    final_rss = rss_per_subset[best_size - 1];
    final_gcv = best_gcv;
  } else {
    // pmethod == "none": keep all M terms. The maintained outer (R, Qty)
    // is the QR of B[:, 0:M] and `rss` is the residual norm² already; just
    // back-substitute for β.
    selected_idx.resize(M);
    std::iota(selected_idx.begin(), selected_idx.end(), 0);
    std::vector<double> beta_full(M, 0.0);
    ares::solve_R_back(R.data(), nk_cap, M, Qty.data(), beta_full.data());
    final_beta = beta_full;
    final_rss = rss;
    final_gcv = compute_gcv(rss, M);
    rss_per_subset[M - 1] = rss;
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
    _["minspan"]        = ms_reported,
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
