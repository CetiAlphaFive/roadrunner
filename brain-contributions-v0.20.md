# Brain Contributions — ares v0.0.0.9001 through v0.0.0.9020

## Proposed Entries

---

### Entry 1: Friedman Fast-LS Scoring via Prefix Sums for MARS Forward Pass

# Friedman Fast-LS Scoring via Prefix Sums for MARS Forward Pass

<!-- brain-entry -->
<!-- domain: builder -->
<!-- subdomain: math-methods -->
<!-- tags: MARS, hinge regression, prefix sums, knot scoring, fast-LS, computational geometry, forward pass, O(n) amortization -->
<!-- contributor: @jackrametta -->
<!-- contributed: 2026-05-10 -->

## Summary

The MARS forward-pass inner loop scores each candidate knot by computing
projections of hinge-pair columns onto the current basis. The naive approach
costs O(n * Mq) per knot; using prefix sums over rows sorted by the candidate
predictor reduces the per-pair cost from O(K * n * Mq) to O((n + K) * Mq),
where K is the number of candidate knots per (parent, variable) pair.

## Knowledge

For each (parent, variable) pair, rows are sorted ascending by x_j. Five scalar
prefix-sum accumulators and two length-Mq vector accumulators are maintained:

    Scalar: S_pr, S_pxr, S_pp, S_ppx, S_ppxx
    Vector: S_Qp[Mq], S_Qpx[Mq]

where p = parent column value, x = predictor value, r = residual, Q = columns
of the current orthonormal basis factor. Totals over all eligible rows (T_*)
are precomputed in one O(n * Mq) pass.

At knot position k (threshold t = x_k), high-side sums are recovered as:
    H_* = T_* - S_* - (contribution of row k)

The key per-knot quantities follow in O(Mq):

    h+ . r = H_pxr - t * H_pr
    h- . r = t * S_pr - S_pxr
    ||h+||^2 = H_ppxx - 2t * H_ppx + t^2 * H_pp
    ||h-||^2 = t^2 * S_pp - 2t * S_ppx + S_ppxx
    (Q' h+)_q = H_Qpx_q - t * H_Qp_q
    (Q' h-)_q = t * S_Qp_q - S_Qpx_q

RSS reduction from adding the hinge pair is then computed via 2-step
Gram-Schmidt: add h+ first, then project h- orthogonal to h+ (already
orthogonal to span(Q)). The reduction formula is:

    red_p = (h+ . r)^2 / ||h+_perp||^2
    alpha = (h+_perp . h-_perp) / ||h+_perp||^2
    e2    = h-_perp - alpha * h+_perp
    red_m = (e2 . r)^2 / ||e2||^2
    total_red = red_p + red_m

Note that h+ and h- have disjoint support (h+(i) > 0 iff x_i > t, h-(i) > 0
iff x_i < t), so (Q' h+)_q and (Q' h-)_q capture the bulk of the
cross-term through their Q-projections, not through direct inner products.

## When to Use

- Implementing the MARS forward-pass inner loop (Friedman 1991, §3.5) in any
  language. The approach is algebraically general and applies to weighted MARS
  as well (replace p = parent_col with p = parent_col * sqrt(weight)).
- Whenever K candidate knots are large relative to n (e.g. n = 1000, K ~ 200)
  and the naive O(K * n * Mq) loop is the bottleneck.
- The prefix-sum approach naturally auto-vectorises over the Mq dimension
  (inner q-loop is unit-stride and independent), making it AVX2/NEON friendly.

## Example

```cpp
// Pseudo-code for prefix-sum MARS knot scoring
// Precompute totals over eligible rows in one O(n*Mq) pass
double T_pr = 0, T_pxr = 0, T_pp = 0, T_ppx = 0, T_ppxx = 0;
std::vector<double> T_Qp(Mq, 0), T_Qpx(Mq, 0);
for (int i = 0; i < n_eli; ++i) {
  int r = sorted_idx[i];
  // accumulate T_*
}

// Sweep ascending x_j; accumulate low-side S_* before scoring each knot
double S_pr = 0, S_pxr = 0, S_pp = 0, S_ppx = 0, S_ppxx = 0;
std::vector<double> S_Qp(Mq, 0), S_Qpx(Mq, 0);
for (int k = 0; k < n_eli; ++k) {
  if (is_candidate_knot(k)) {
    double t = x[sorted_idx[k]];
    // Recover high-side as Total - Low - row_k
    double H_pr = T_pr - S_pr - p[k]*r[k];
    // ... compute rss_reduction from H_* and S_* in O(Mq) ...
  }
  // accumulate row k into S_*
}
```

## Pitfalls

- The precomputed `var_sort` (sort once per variable per forward step) should
  be hoisted entirely outside the forward loop if x does not change (it
  depends only on the design matrix, not on the residual). Failure to hoist
  produces O(nk_cap) redundant sorts — a measured 28% wall-clock fraction on
  small designs in profiling.
- Prefix sums are built over eligible rows only (parent support), in ascending
  x_j order. Mixing in ineligible rows or wrong sort order silently corrupts
  every knot score.
- The high-side recovery (T - S - row_k) requires that row k itself is
  excluded from both the low-side prefix sum S and the high-side sum. Score
  the knot before accumulating row k.
- On near-rank-deficient designs (||h+_perp||^2 near zero), guard with a
  threshold (e.g. 1e-12) before dividing.

---

### Entry 2: Incremental QR Maintenance for MARS Forward Pass

# Incremental QR Maintenance for MARS Forward Pass

<!-- brain-entry -->
<!-- domain: builder -->
<!-- subdomain: math-methods -->
<!-- tags: QR factorization, incremental update, Householder, rank-deficiency, forward pass, linear algebra, column append -->
<!-- contributor: @jackrametta -->
<!-- contributed: 2026-05-10 -->

## Summary

In a greedy forward-selection algorithm (MARS, stepwise regression, etc.),
recomputing the full QR factorization of the basis matrix at each step costs
O(n * M^2) per step. Maintaining (Q1, R, Qty) incrementally via column-append
reduces the per-step marginal cost to O(n * M), with a full QR refresh only
needed at specific downstream transition points (e.g., backward pruning pass).

## Knowledge

Maintain the thin QR factorization: B = Q1 * R, where Q1 is n x M orthonormal,
R is M x M upper-triangular, and Qty = Q1' * y.

When a new column c is appended to B:

    u        = Q1' c           (O(n*M) matrix-vector product)
    c_perp   = c - Q1 * u      (O(n*M))
    d        = ||c_perp||

Rank guard: if d < eps_rank * ||c|| (e.g. eps_rank = 1e3 * machine_epsilon):
  - Mark column as rank-deficient; zero out the new Q1 column, set R[M,M] = 0.
  - This keeps M_q (number of valid columns) in step with M without gaps.

Otherwise:
    Q1[:, M] = c_perp / d
    R[0:M, M] = u;  R[M, M] = d
    Qty[M]   = Q1[:, M]' * y  (= (c_perp' * y) / d, closed-form from above)
    resid    -= Q1[:, M] * Qty[M]
    rss      -= Qty[M]^2

For the hot knot-scanning inner loop, store Q1 in TRANSPOSED layout (QT where
QT[row * stride + q] = Q1[row, q]) so the q-loop is unit-stride and
auto-vectorisable by the compiler or by hand with SIMD intrinsics.

The backward pruning pass must NOT inherit the forward-maintained QR because:
1. Rank-deficient columns produce zero diagonals in R that Givens downdate
   cannot handle cleanly.
2. Backward drops need the full M x M R, not the M_q x M_q reduced form.
Perform a fresh Householder QR on the surviving forward-pass basis before
entering backward. The cost is O(n * M^2), which profiling typically places
at under 1% of total wall-clock.

## When to Use

- Any greedy column-addition algorithm where the same orthonormal basis Q is
  needed for multiple subsequent operations (residual projection, Gram-Schmidt
  scoring) and the basis grows one column at a time.
- MARS forward pass, forward-stepwise regression, Matching Pursuit, OMP.
- When the basis dimension M is small (M <= 200) relative to n (n >= 500):
  the incremental approach cuts multi-thread scaling bottlenecks because the
  serial O(n * M^2) Gram-Schmidt rebuild is replaced with O(n * M) per step.

## Example

```cpp
// Incremental Q1 append
// Q1: n x current_M (stored column-major)
// R:  current_M x current_M (upper triangular)
// Qty: length current_M

void append_column(const double* c, int n, int M,
                   double* Q1, double* R, double* Qty,
                   double* resid, double& rss) {
  std::vector<double> u(M);
  for (int q = 0; q < M; ++q) {
    double dot = 0;
    for (int i = 0; i < n; ++i) dot += Q1[i + q*n] * c[i];
    u[q] = dot;
  }
  std::vector<double> c_perp(c, c + n);
  for (int q = 0; q < M; ++q)
    for (int i = 0; i < n; ++i)
      c_perp[i] -= u[q] * Q1[i + q*n];
  double d = 0;
  for (int i = 0; i < n; ++i) d += c_perp[i] * c_perp[i];
  d = std::sqrt(d);
  if (d < eps_rank) { /* rank-deficient: zero column */ return; }
  for (int i = 0; i < n; ++i) Q1[i + M*n] = c_perp[i] / d;
  for (int q = 0; q < M; ++q) R[q + M*(M+1)] = u[q];
  R[M + M*(M+1)] = d;
  double qty_M = 0;
  for (int i = 0; i < n; ++i) qty_M += Q1[i + M*n] * resid[i] + 0; // resid is y - Q1*Qty
  // Actually: qty_M = Q1[:,M]' * y (precomputed from c_perp' * y before normalising)
  Qty[M] = qty_M;
  for (int i = 0; i < n; ++i) resid[i] -= Q1[i + M*n] * qty_M;
  rss -= qty_M * qty_M;
}
```

## Pitfalls

- The maintained RSS can drift below zero on near-singular designs if
  rank-deficient columns accumulate Qty updates. Guard: clamp rss = max(rss, 0).
- The transposed Q layout (QT) must use a stride >= nk_cap (the maximum
  possible M), not the current M, so the active region lives contiguously
  without reallocation or repacking as M grows.
- Incremental Q tracks the TRUE OLS minimum RSS at each M, whereas a fresh
  Householder QR on the same (potentially rank-clamped) basis tracks the
  pseudo-minimum. On near-rank-deficient designs, incremental Q gives lower RSS
  than competitors using full refactorization — this is correct behaviour, but
  will cause divergences in downstream parity tests against implementations that
  clamp rank-deficient columns to zero.

---

### Entry 3: Givens Downdate for Backward Subset Selection with Periodic Householder Refresh

# Givens Downdate for Backward Subset Selection with Periodic Householder Refresh

<!-- brain-entry -->
<!-- domain: builder -->
<!-- subdomain: math-methods -->
<!-- tags: Givens rotations, QR downdate, backward selection, GCV, numerical stability, periodic refresh -->
<!-- contributor: @jackrametta -->
<!-- contributed: 2026-05-10 -->

## Summary

In GCV-based backward subset selection (MARS backward pass, stepwise deletion),
each candidate term removal requires recomputing the OLS RSS after deletion.
Givens downdating reduces the per-candidate cost from O(n * M^2) (full QR) to
O(M^2). Applying Givens downdates cumulatively across all M backward steps
introduces numerical drift; refreshing with a full Householder QR every 4 steps
eliminates drift at low amortised cost.

## Knowledge

**Per-candidate Givens downdate:**

Given the current R (M x M upper triangular) and Qty (length M), simulate
column k being removed from the basis:
1. Copy R and Qty to R_trial, Qty_trial.
2. Delete row k and column k from R_trial (renumber remaining rows/cols).
3. Restore upper-triangular form of R_trial using (M - k - 1) Givens rotations
   on rows k through M-1, applying the same rotations to Qty_trial.
4. RSS_without_k = ||y||^2 - ||Qty_trial||^2

This is O(M^2) per candidate, versus O(n * M^2) for a full QR. With M candidates
to evaluate, the step cost drops from O(n * M^3) to O(M^3).

**Periodic Householder refresh:**

After committing a drop (deleting the chosen column k), apply Givens downdate
to the live R and Qty. Every r steps (empirically r = 4 works well), perform
a fresh Householder QR on the surviving basis columns. This refresh:
- Costs O(n * M^2) once per r steps, amortised to O(n * M^2 / r) per step.
- Eliminates cumulative rotation drift that would otherwise cause the backward
  pass to pick suboptimal drops (producing 2-3% RSS errors relative to ground truth).

**The cadence r = 4 is a sweet spot**: r = 1 (refresh every step) reproduces
exactly the O(n * M^3) cost of full QR. r = infinity (no refresh) improves
1-thread wall-clock ~6% but causes multi-thread regression and some cells to
cross back above earth's cost. r = 4 cuts per-step cost while keeping drift
below FP test tolerance across a wide range of designs.

## When to Use

- Any backward subset-selection loop where candidate-deletion RSS needs to be
  evaluated cheaply (MARS, forward-stepwise, subset selection in linear models).
- When M is large enough that O(n * M^2) per candidate is the bottleneck but
  you cannot afford cumulative drift from pure Givens chains.
- NOT suitable when R has zero diagonal entries (rank-deficient columns). The
  backward pass must begin from a fresh Householder QR on the survived
  forward-pass terms, rejecting rank-deficient columns before entering the
  downdate loop.

## Example

```cpp
// Pseudo-code: backward step with Givens downdate + periodic refresh
// Inputs: R (M x M), Qty (M), y_norm2 = ||y||^2
// Maintained state: refresh_counter (reset to 0 after each Householder refresh)

double best_rss_increase = std::numeric_limits<double>::infinity();
int    best_drop = -1;

for (int k = 1; k < M; ++k) {  // k=0 is intercept, never dropped
  std::vector<double> Rt = R;   // copy
  std::vector<double> Qt = Qty;
  givens_delete_col(Rt, Qt, k, M);
  double rss_k = y_norm2 - dot(Qt, Qt);
  double increase = rss_k - current_rss;
  if (increase < best_rss_increase) {
    best_rss_increase = increase;
    best_drop = k;
  }
}
// commit drop
commit_givens_delete(R, Qty, best_drop, M);
--M;
++refresh_counter;
if (refresh_counter >= 4) {
  householder_qr_refresh(B_current, n, M, R, Qty);
  refresh_counter = 0;
}
```

## Pitfalls

- Zero diagonals in R (from rank-deficient column appends in the forward pass)
  cause Givens downdate to silently produce wrong RSS values. Always start the
  backward loop from a fresh QR that rejects rank-deficient columns.
- Closed-form Cholesky downdate (DELTA_RSS_j = beta_j^2 / [(R'R)^{-1}]_jj)
  is an attractive O(M) alternative but is unstable on rank-deficient designs.
  Tested and rejected; Givens is more robust.
- The refresh cadence (r = 4) is empirical. On strongly interaction-heavy
  designs or high M (>100), r = 2 or r = 3 may be needed; validate by comparing
  with a ground-truth full-QR backward pass on a representative held-out cell.

---

### Entry 4: AVX2 SIMD Cross-Lane Reduction for MARS Inner Loop

# AVX2 SIMD Cross-Lane Reduction for MARS Inner Loop

<!-- brain-entry -->
<!-- domain: builder -->
<!-- subdomain: math-methods -->
<!-- tags: AVX2, SIMD, vectorisation, C++, inner loop, accumulator, cross-lane reduction, FMA-safe -->
<!-- contributor: @jackrametta -->
<!-- contributed: 2026-05-10 -->

## Summary

The MARS knot-scoring q-loop (accumulating dot products over the current basis
dimension Mq) is amenable to AVX2 manual vectorisation with 4-lane double
accumulators. The pattern processes 4 q-indices per iteration, reduces across
lanes once per knot (not once per lane update), and is safe with FMA enabled
(`-ffp-contract=on`) because it uses `_mm256_mul_pd` + `_mm256_add_pd` rather
than fused intrinsics.

## Knowledge

The loop computes three accumulated sums over q in [0, Mq):

    acc_pp += dp[q]^2
    acc_mm += dm[q]^2
    acc_pm += dp[q] * dm[q]

where dp[q] and dm[q] are scalar expressions depending on precomputed prefix
sums and the current Q row.

AVX2 pattern:

```cpp
__m256d acc_pp = _mm256_setzero_pd();
__m256d acc_mm = _mm256_setzero_pd();
__m256d acc_pm = _mm256_setzero_pd();
int q = 0;
for (; q + 4 <= Mq; q += 4) {
  __m256d dp = /* compute 4-lane dp */;
  __m256d dm = /* compute 4-lane dm */;
  acc_pp = _mm256_add_pd(acc_pp, _mm256_mul_pd(dp, dp));
  acc_mm = _mm256_add_pd(acc_mm, _mm256_mul_pd(dm, dm));
  acc_pm = _mm256_add_pd(acc_pm, _mm256_mul_pd(dp, dm));
}
// Tail scalar fallback for Mq % 4 != 0
for (; q < Mq; ++q) { /* scalar accumulate */ }
// Cross-lane horizontal reduction
alignas(32) double rbuf[4];
_mm256_store_pd(rbuf, acc_pp);
double cp_correction = rbuf[0] + rbuf[1] + rbuf[2] + rbuf[3];
_mm256_store_pd(rbuf, acc_mm);
double cm_correction = rbuf[0] + rbuf[1] + rbuf[2] + rbuf[3];
_mm256_store_pd(rbuf, acc_pm);
double pm_correction = rbuf[0] + rbuf[1] + rbuf[2] + rbuf[3];
// Apply corrections to pre-computed scalar norms
cp_perp_norm2 -= cp_correction;
cm_perp_norm2 -= cm_correction;
cp_perp_dot_cm_perp -= pm_correction;
```

Guard with `#if defined(__AVX2__)` so the package still compiles on
non-AVX2 targets; the scalar fallback path (identical arithmetic) follows
the SIMD block.

**FMA safety:** Using `_mm256_mul_pd + _mm256_add_pd` (two separate instructions)
is not the same as `_mm256_fmadd_pd`. With `-ffp-contract=on`, the compiler may
fuse nearby scalar multiplications into FMAs in the scalar path, but the manual
SIMD path with explicit mul+add is not subject to that contraction. This means
the SIMD and scalar paths can give slightly different per-knot reductions in
principle; in practice, on the designs tested, this did not break run-to-run
determinism.

The hot precompute loop (accumulating T_Qp and T_Qpx totals over all rows) uses
the same AVX2 pattern and benefits proportionally.

## When to Use

- Inner loops accumulating dot products or squared norms over a moderate-length
  vector (Mq in the range 5-200) that runs inside an outer loop executed
  thousands of times per fit (once per candidate knot, per (parent, variable) pair).
- C++ code compiled with `-march=native` on x86-64 hardware; guard with
  `__AVX2__` or a CMake/Makevars feature test.
- When FMA can be left enabled (`-ffp-contract=on`): explicit mul+add in the
  SIMD path avoids the scalar-path FMA contraction, preserving determinism.

## Example

See Knowledge section above for the full pattern.

## Pitfalls

- Do NOT store intermediate per-lane values to memory inside the SIMD loop
  (e.g. storing qcp[q] or qcm[q]). If those stores are dead code (values never
  read after the loop), the compiler may not eliminate them and the extra writes
  add cache pressure that outweighs SIMD gains. Profile before adding any
  per-lane writes.
- The scalar fallback tail (`for (; q < Mq; ++q)`) must use the SAME arithmetic
  as the scalar path that would have run without SIMD; otherwise the SIMD and
  non-SIMD totals diverge, breaking determinism checks.
- Horizontal reduction (summing 4 lanes) via `_mm256_store_pd + scalar sum` is
  portable. `_mm256_hadd_pd` requires two instructions and is not faster on
  modern microarchitectures for a one-time 4-element sum.
- On ARM (NEON/SVE), the same pattern applies with `float64x2_t`; the algebra
  is identical, only the intrinsic names differ.

---

### Entry 5: Fast-MARS Priority Cache — Winner-from-Rescored-Set Rule

# Fast-MARS Priority Cache — Winner-from-Rescored-Set Rule

<!-- brain-entry -->
<!-- domain: builder -->
<!-- subdomain: math-methods -->
<!-- tags: MARS, priority cache, fast-MARS, Friedman 1993, stale scores, forward pass, greedy search -->
<!-- contributor: @jackrametta -->
<!-- contributed: 2026-05-10 -->

## Summary

The Friedman (1993) fast-MARS priority cache reduces forward-pass work by
re-scoring only a subset of (parent, variable) pairs each step. The critical
correctness rule is: pairs from the stale cache must NEVER be picked as the
step winner directly from their cached score. The winner must always come from
the rescored (freshly evaluated) set. Violating this rule collapses model
quality (R^2 dropping from ~0.99 to ~0.35 was observed in testing).

## Knowledge

**Cache structure:**
Each (parent, variable) pair maintains a `Candidate` struct with:
- The best RSS reduction observed at its last rescore.
- An age counter (steps since last rescore).
- An age-discounted score: score * discount^age (Friedman uses discount = 1/(1 + beta*age)).

**Per-step protocol:**
1. Always rescore pairs whose parent was added in the immediately preceding step
   (they are "fresh" — their parent column changed).
2. Among remaining stale pairs, rank by age-discounted cached score; rescore
   the top fast_k pairs (fast_k = 10 is a good default).
3. Pick the winner as argmax(RSS reduction) over the RESCORED set only.
4. Update cache entries for all rescored pairs.

**Why the winner-from-rescored-set rule is non-negotiable:**
Stale scores reflect the RSS reduction from a previous step, before new basis
columns were added. Adding a hinge pair whose support was already well-explained
by the new columns would produce near-zero actual reduction but a high stale
score. Picking from stale scores therefore selects mis-scored pairs, causing the
forward pass to terminate prematurely with a near-trivial model.

**fast_k tuning:**
- fast_k = 0: disable cache; every pair is rescored every step. Slowest, highest
  accuracy.
- fast_k = 5: most aggressive cache. ~3-5% accuracy loss on tightly-tied designs.
- fast_k = 10 (default): well-balanced; measured median accuracy loss under 0.5%
  across a variety of DGPs.
- fast_k = 25: almost as fast as fast_k = 10 on large p, negligible further loss.

**Autotune application:** When sweeping fast_k in a grid, use the rule "pick
the smallest non-zero fast_k whose CV-MSE is within 1% of the best-cell CV-MSE;
fall back to fast_k = 0 only if no positive fast_k qualifies." This yields the
cheapest setting that does not measurably hurt accuracy on the actual data.

## When to Use

- Greedy forward-selection algorithms where pair scoring is expensive and the
  number of pairs is large (p * M pairs grow quadratically in p and the step count M).
- Effective when many pairs share a variable but not a parent (most pairs are
  ineligible at any given step), so stale scores are a reasonable proxy for
  relative importance.
- Particularly effective at high p (many predictors) or high nk (many forward steps).

## Example

```r
# R interface (genericized)
fit <- forward_selector(x, y, fast_k = 10, fast_beta = 1.0)

# The fast_k = 0 path matches the exact (slow) forward pass for correctness
# verification:
fit_exact <- forward_selector(x, y, fast_k = 0)
stopifnot(abs(fit$rss - fit_exact$rss) / fit_exact$rss < 0.01)
```

## Pitfalls

- Do not pick the winner from cached/stale scores. The rescore set (fresh pairs
  + top-fast_k stale pairs) must always contain the eventual winner.
- A single parent being added in step t makes ALL its child pairs fresh for
  step t+1. If many parents share a variable, the fresh set can be large;
  budget for this in time estimates.
- On low-n designs (n < 200) with few predictors, the cache provides little
  benefit because the number of pairs is small and rescore cost is negligible.
  The overhead of age tracking can hurt; consider gating the cache on
  n_pairs > some threshold (e.g., 30).

---

### Entry 6: Successive Halving for Hyperparameter Grid Search with Fold-1 Screening

# Successive Halving for Hyperparameter Grid Search with Fold-1 Screening

<!-- brain-entry -->
<!-- domain: builder -->
<!-- subdomain: math-methods -->
<!-- tags: hyperparameter tuning, successive halving, cross-validation, early stopping, grid search, model selection -->
<!-- contributor: @jackrametta -->
<!-- contributed: 2026-05-10 -->

## Summary

When running a multi-fold CV grid search over a large hyperparameter space,
evaluating every cell on all folds wastes most of the budget on clearly inferior
configurations. A simple successive-halving rule — eliminate cells whose fold-1
MSE exceeds 1.5x the running minimum after fold 1 — cuts wall-clock by 40-60%
on designs where one model class dominates, with negligible accuracy loss.

## Knowledge

**Protocol:**

1. Score all cells on fold 1 of the inner CV.
2. Compute: threshold = 1.5 * min(fold-1 MSE across all cells).
3. Eliminate any cell with fold-1 MSE > threshold. Mark it as eliminated; record
   its partial score (fold-1 mean) and skip the remaining folds.
4. Score surviving cells on the remaining (nfold * ncross - 1) fold-evaluations.
5. Final ranking uses only full-grid cells (non-eliminated) for the winner.

**Elimination factor 1.5:**
- Cells within 50% of the leader on fold 1 always survive.
- This is deliberately loose. A 1.1x threshold would miss some good cells due
  to fold-1 sampling noise; 1.5x is empirically safe across a range of DGPs.
- On strongly interaction-heavy targets: degree-1 cells are typically eliminated
  after fold 1, saving (nfold * ncross - 1) evaluations per eliminated cell.
- On linear/additive targets: degree-2 and degree-3 cells may survive if the
  fold-1 MSE gap is small; they proceed to full scoring.

**Integration with parsimony tie-break:**
If two cells have the same final CV-MSE, prefer the smaller degree, then smaller
nk, then smaller penalty. This ensures successsive halving doesn't bias the
winner toward complex models that happen to have lower fold-1 MSE by noise.

**Cost:**
With G cells and nfold * ncross folds, naive cost is G * nfold * ncross
evaluations. With halving eliminating fraction f of cells after fold 1, cost
drops to G + (1-f) * G * (nfold * ncross - 1). At f = 0.5, this is a ~2x
speedup at large nfold * ncross.

## When to Use

- Inner-loop CV grid search for any model where one model class (e.g. linear
  vs. nonlinear, simple vs. complex) is consistently better across datasets.
- Particularly effective when the grid has a "nuisance" dimension (e.g. degree=1
  vs. degree=2 in regression when the true DGP is highly nonlinear): fold-1
  is enough to separate the classes.
- NOT reliable when fold-1 MSE is very noisy (e.g. nfold = 3, n = 100, p = 50).
  In those settings, use a larger minimum_survivors count or raise the threshold.

## Example

```r
# Pseudo-code: successive halving in CV grid search
fold1_mse <- score_all_cells_on_fold1(grid, data)
threshold <- 1.5 * min(fold1_mse)
eliminated <- fold1_mse > threshold

for (fold in 2:(nfold * ncross)) {
  active_cells <- which(!eliminated)
  fold_mse[active_cells, fold] <- score_cells_on_fold(active_cells, data, fold)
}

# Compute mean over all scored folds (ignoring eliminated cells' partial scores)
final_mse <- rowMeans(fold_mse[!eliminated, ], na.rm = FALSE)
winner <- which.min(final_mse)  # then apply parsimony tie-break
```

## Pitfalls

- Eliminated cells must not participate in the final ranking. Their fold-1
  score is biased (only one fold sampled) and will tend to underestimate true MSE.
- Report the `eliminated` flag in the output grid for diagnostics. Users should
  be able to see which cells were screened out and on what basis.
- The halving threshold is DATA-DEPENDENT. A design that looks linear on a 15%
  subsample may look nonlinear on the full data. If using a subsample pre-probe
  (warm-start), re-run the full grid if the subsample result is not decisive.

---

### Entry 7: Subsample Warm-Start to Short-Circuit Hyperparameter Grid Search

# Subsample Warm-Start to Short-Circuit Hyperparameter Grid Search

<!-- brain-entry -->
<!-- domain: builder -->
<!-- subdomain: math-methods -->
<!-- tags: hyperparameter tuning, warm start, subsample, early exit, grid search, MARS, autotune -->
<!-- contributor: @jackrametta -->
<!-- contributed: 2026-05-10 -->

## Summary

A cheap subsample pre-fit can identify an unambiguous winner among model
families before running the expensive full-data grid search. If the best model
family wins by a wide enough margin on the subsample (here: >= 5-10% CV-MSE gap
between the best per-family cell and the next-best family's best cell), adopt
the subsample winner and skip the full grid entirely.

## Knowledge

**Protocol:**

1. Draw a subsample of min(0.15 * n, 200) rows (stratified by y-quantile).
2. Run a coarse grid search on the subsample using a fast setting (fewer folds,
   smaller nk) — this probe is sub-second even for large n.
3. Compute: for each model family (e.g. degree), find its best-cell CV-MSE.
4. Decisiveness check: if min(best_per_family) < (1 - gap_threshold) * second_best_per_family,
   the subsample has identified an unambiguous winner family.
5. If decisive: return the subsample winner's hyperparameters and fit the full
   data once with those params. Skip the full grid.
6. If not decisive: run the full-data grid as normal. The probe adds ~0.4s of
   overhead but does not degrade the full-grid outcome.

**Gap threshold:**
- 5% gap (0.95x ratio): aggressive short-circuiting, slight risk of wrong family.
- 10% gap: safe in practice across tested DGPs.
- The threshold should be tuned per problem domain; for regression on smooth
  structured DGPs a 5-10% gap is reliable.

**Subsample size:**
- 15% of n, capped at 200. Below 200 rows the fit is noisy but fast; above 200
  the probe overhead starts to matter. On n = 600 with a 200-row subsample,
  the probe fires (or doesn't) in ~0.1s.

**Effect:**
- On a clearly-linear DGP (n=600, p=6): warm-start fires and reduces total
  autotune time from ~15s to ~0.5s (30x speedup).
- On a strongly-interaction DGP: probe is not decisive, full grid runs; overhead
  is ~0.4s out of a total ~2.5s run.

## When to Use

- Any multi-cell grid search where the number of cells is large (>20) and the
  dominant distinction is between model families that are reliably separated
  on partial data.
- Particularly effective when degree (linear vs. nonlinear) is the dominant
  hyperparameter dimension: a 200-row subsample is typically enough to detect
  strong nonlinearity.
- NOT reliable for identifying the best penalty or nk within a family — those
  dimensions vary too subtly for a 200-row sample to resolve reliably. Use
  the warm-start only to select the family; sweep fine-grained params on the
  full data.

## Example

```r
# Pseudo-code: subsample warm-start in autotune
if (n >= 200) {
  n_sub <- min(round(0.15 * n), 200)
  idx_sub <- stratified_sample(y, n_sub)
  probe_result <- run_grid(x[idx_sub, ], y[idx_sub], speed = "fast", nfold = 3)
  best_per_family <- aggregate_best_by_family(probe_result$grid)
  if (is_decisive(best_per_family, gap_threshold = 0.10)) {
    winner_params <- extract_winner(probe_result)
    return(fit_full(x, y, params = winner_params))
  }
}
# Fall through to full grid search
full_result <- run_grid(x, y, speed = "balanced", nfold = 5)
```

## Pitfalls

- The warm-start result must be discarded (not used as a partial fold in the
  full-grid search) to avoid contaminating the CV-MSE estimates.
- Report a flag (e.g. `$warmstart_fired = TRUE`) in the output so the user can
  see that the full grid was skipped. Users relying on autotune for reproducible
  comparisons should be able to disable warm-start.
- Seed the subsample draw separately from the main CV fold partition so the two
  streams are reproducible AND distinct. Using `seed_cv + constant` for the
  subsample RNG works well.

---

### Entry 8: Shared Forward Pass Across Hyperparameter Grid Cells

# Shared Forward Pass Across Hyperparameter Grid Cells

<!-- brain-entry -->
<!-- domain: builder -->
<!-- subdomain: math-methods -->
<!-- tags: MARS, shared forward pass, hyperparameter tuning, grid search, backward replay, model selection, GCV -->
<!-- contributor: @jackrametta -->
<!-- contributed: 2026-05-10 -->

## Summary

In MARS (and analogous greedy forward-selection algorithms), the forward pass
output (selected basis structure: term directions and knot positions) depends
only on (degree, nk, fast_k) and is INDEPENDENT of the penalty parameter that
governs backward pruning. Cells that share (degree, nk, fast_k) can therefore
run the expensive forward pass once and replay the cheaper backward pass with
each cell's distinct penalty, yielding a 1.5-2.3x speedup in hyperparameter
grid search.

## Knowledge

**Factorisation:**

In GCV-based MARS:
- Forward pass: greedy hinge selection up to nk terms. Depends on (degree, nk,
  fast_k). DOES NOT depend on penalty.
- Backward pass: sequential subset selection minimising GCV. Depends on penalty
  (which enters the GCV denominator). Also depends on the forward-pass output
  (dirs, cuts). Cost: O(n * M^2 + M^4) for Householder QR init + Givens downdate loop.

Group cells by (degree, nk, fast_k). For each group and each CV fold:
1. Run forward pass once on the fold's training data. Save (dirs, cuts).
2. Build the full basis matrix B from (dirs, cuts) on the training data.
3. For each cell in the group, replay the backward pass with that cell's penalty.
   Each backward replay costs O(n * M^2 + M^4) — the same as the backward
   block of a full fit, but with no forward cost.

The backward replay function must:
- Accept (x, y, dirs, cuts, penalty) as inputs.
- Rebuild B from (dirs, cuts) on the given x.
- Run the same Householder QR init + Givens downdate + periodic refresh as the
  full backward pass.
- Return (coefficients, rss, gcv, rss_per_subset, gcv_per_subset).

For CV-based subset selection (pmethod = "cv"), the backward replay can
optionally return per-size path data so the caller can score every candidate
subset size on the holdout in O(M) per size.

## When to Use

- Grid search over a hyperparameter space with a penalty dimension when other
  hyperparameter dimensions are fixed or few in number.
- The forward pass dominates total grid-search cost (typical when n is large
  and M is moderate).
- In general: any model where the "structure selection" step (expensive, runs
  once) can be separated from the "shrinkage/complexity control" step (cheap,
  runs per cell).

## Example

```cpp
// Pseudo-code: backward-only C++ function for grid replay
// mars_backward_only(x, y, dirs, cuts, penalty, nprune, nthreads):
//   1. Build B from dirs/cuts on x  -- O(n * M)
//   2. Householder QR of B          -- O(n * M^2)
//   3. Givens downdate backward loop -- O(M^3)
//   4. GCV selection, return coefficients and metadata
```

```r
# R side: group cells and dispatch
groups <- split(grid, list(grid$degree, grid$nk, grid$fast_k))
for (grp in groups) {
  for (fold in 1:nfold) {
    dirs_cuts <- run_forward(x_train[[fold]], y_train[[fold]],
                             degree = grp$degree[1], nk = grp$nk[1],
                             fast_k = grp$fast_k[1])
    for (cell in grp) {
      fold_mse[cell, fold] <- backward_only_and_score(
        x_train[[fold]], y_train[[fold]], x_val[[fold]], y_val[[fold]],
        dirs_cuts$dirs, dirs_cuts$cuts, penalty = cell$penalty)
    }
  }
}
```

## Pitfalls

- Numerical identity with the full-fit backward: the backward replay MUST use
  exactly the same algorithm (same QR implementation, same Givens downdate, same
  periodic refresh cadence) as the full-fit backward block, or CV-MSE estimates
  will diverge between the shared-forward and per-cell-full-fit paths. This is
  easy to verify: at the same penalty, backward_only(dirs_cuts) should return
  the same GCV as mars_fit(x, y, penalty).
- The basis matrix B rebuilt from (dirs, cuts) inside the backward replay must
  match the B that the forward pass would have built on the same data. Verify
  with a unit test before using the replay path for model selection.
- If nk varies across grid cells (e.g. nk in {10, 20, 40}), cells with different
  nk cannot share a forward pass even if they share (degree, fast_k). Group only
  by the full tuple (degree, nk, fast_k).

---

### Entry 9: Bootstrap Bagging for MARS Uncertainty Quantification

# Bootstrap Bagging for MARS Uncertainty Quantification

<!-- brain-entry -->
<!-- domain: builder -->
<!-- subdomain: math-methods -->
<!-- tags: bagging, bootstrap, uncertainty quantification, MARS, ensemble, prediction intervals, R -->
<!-- contributor: @jackrametta -->
<!-- contributed: 2026-05-10 -->

## Summary

MARS (like most tree/greedy methods) provides point estimates with no native
uncertainty measure. Row-bootstrap bagging over B replicates is a simple and
effective add-on: average replicate predictions for bias reduction, and report
the per-row standard deviation across replicates as a rough prediction
uncertainty band. The critical design point is that each replicate uses the
SAME hyperparameters as the central fit (not re-tuned), so autotune cost is
paid only once.

## Knowledge

**Protocol:**

1. Fit the central model (with autotune if desired) to the full training data.
   Record the chosen hyperparameters (degree, penalty, nk, fast_k).
2. For b in 1..n_boot:
   a. Draw a row-bootstrap sample (with replacement, n rows).
   b. Fit the model using the central fit's hyperparameters (no re-tuning).
   c. Store the fitted object in a list.
3. Bagged prediction at new x:
   - Compute predictions from the central fit and all n_boot replicates.
   - Return mean as the point estimate.
   - Return SD across (n_boot + 1) predictions as `se` (accessible via
     `se.fit = TRUE` in predict).

**Composition with autotune:**
- Autotune runs ONCE on the central fit.
- Each bootstrap replicate refits with the tuned params. Cost is
  (n_boot + 1) * single_fit_cost, not n_boot * autotune_cost.
- This makes autotune + bagging practical: for n_boot = 50 and a ~0.5s fit,
  total cost is ~25s regardless of autotune complexity.

**RNG separation:**
- When a reproducibility seed is set, use two separate streams: one for CV
  fold partitioning (seed S), one for bootstrap draws (seed S + constant).
  This ensures (a) both streams are reproducible at fixed seed and (b) changing
  the fold partition does not perturb the bootstrap sequence or vice versa.
- Save and restore the caller's `.Random.seed` around all internal RNG use so
  the user's global RNG stream is unaffected.

## When to Use

- When a MARS fit will be used for decision-making that requires uncertainty
  bounds, and you cannot afford the complexity of a full Bayesian approach.
- When the reference implementation (e.g., `earth`) provides no bagging and you
  need a drop-in uncertainty extension.
- When autotune is already paying the hyperparameter selection cost: adding
  n_boot = 50 replicates at single-fit cost is cheap relative to the autotune
  grid scan.
- NOT a substitute for calibrated conformal prediction intervals; the SD across
  bootstrap replicates captures model variance but not aleatoric uncertainty.

## Example

```r
fit <- model_fit(x, y, autotune = TRUE, n_boot = 50, seed = 42)
preds <- predict(fit, newdata = x_test, se.fit = TRUE)
# preds: vector of point estimates
# attr(preds, "sd"): per-row uncertainty (bootstrap SD)
upper_90 <- preds + 1.645 * attr(preds, "sd")
lower_90 <- preds - 1.645 * attr(preds, "sd")
```

## Pitfalls

- Each replicate is fit on n rows sampled WITH replacement. Out-of-bag rows
  are available for validation but are NOT used internally in this pattern;
  the OOB predictions are a bonus diagnostic, not the primary output.
- Bootstrap SD under-covers when: (a) the model is high-variance (overfitting),
  (b) n is small (n < 200), or (c) the true DGP has heavy-tailed residuals. In
  those settings the SD is a lower bound on true uncertainty.
- Do NOT re-run autotune for each replicate. Using the central fit's tuned
  params for all replicates is correct and efficient. The replicate models may
  select slightly different basis structures (because the bootstrap sample is
  different), but the hyperparameters (degree, penalty, nk) are fixed.

---

### Entry 10: CV Scoring Strategy — Aggregate Only Over Universally-Reached Subset Sizes

# CV Scoring Strategy — Aggregate Only Over Universally-Reached Subset Sizes

<!-- brain-entry -->
<!-- domain: tester -->
<!-- subdomain: math-methods -->
<!-- tags: cross-validation, backward selection, MARS, subset size, CV-MSE aggregation, bias, fold coverage -->
<!-- contributor: @jackrametta -->
<!-- contributed: 2026-05-10 -->

## Summary

In K-fold CV for models with an internal sequential pruning step (MARS backward
pass, stepwise regression), different folds may produce different maximum model
sizes because the forward pass terminates differently on each training subset.
Averaging MSE only over sizes reached by ALL folds — not just some folds —
removes a subtle downward bias: folds that terminate at smaller sizes also tend
to have lower MSE at those sizes (the same regularisation that stops early also
reduces variance).

## Knowledge

**The bias:**

Suppose fold A terminates at M_max = 10 and fold B terminates at M_max = 8.
If you average MSE at size 9 over folds that provide it (fold A only), you
underestimate MSE at size 9 because fold A's data tend to support larger models
(less regularisation needed) and therefore have lower holdout MSE at any given
size. The average across all folds (which includes fold B at smaller sizes) gives
a fairer picture.

**The fix:**

1. Collect per-fold, per-size MSE values: fold_mse[fold, size].
2. Compute n_folds_with_size[size] = number of folds that reached at least that size.
3. For CV-MSE aggregation at each size, use only sizes where n_folds_with_size[size] == nfold.
4. If no size is reached by ALL folds (e.g. very small n), fall back to any-fold
   sizes and report a warning.

**1-SE rule integration:**

After selecting size_star via argmin of mean CV-MSE, compute the 1-SE rule as:
    se_star = cv_se[size_star]
    1se_size = min(sizes) where mean_cv_mse[size] <= mean_cv_mse[size_star] + se_star

The 1-SE size is always <= size_star (never larger), since we're looking for
the smallest model within 1 SE of the minimum.

**Return diagnostics:**

Always return:
- cv_mse: length-max_size vector of mean CV-MSE by size (NA for sizes not
  reached by all folds).
- cv_n: count of folds contributing to each size (0 for universally-unreached).
- cv_se: per-size SE (sqrt(var/nfold) over the fold-level MSE values).
- cv_mse_mat: (nfold x max_size) matrix of per-fold per-size MSE. Essential for
  diagnostics and custom selection rules.

## When to Use

- Any greedy model where the number of selected terms can vary across folds
  (MARS, stepwise regression, lasso with warm-start stopping).
- When the user wants to apply a custom selection rule (not just argmin or
  1-SE), the full cv_mse_mat is needed.
- When nfold is small (3-5) or n is small (< 500), fold-level size variation
  is common; the universal-coverage restriction matters most here.

## Example

```r
# Pseudo-code: CV-MSE aggregation with universal-size restriction
fold_mse <- matrix(NA, nrow = nfold, ncol = max_size)
for (fold in 1:nfold) {
  fold_mse[fold, 1:size_reached[fold]] <- compute_holdout_mse(fold)
}
n_folds_with_size <- colSums(!is.na(fold_mse))
universal_sizes <- which(n_folds_with_size == nfold)
mean_cv_mse <- colMeans(fold_mse[, universal_sizes, drop = FALSE])
size_star <- universal_sizes[which.min(mean_cv_mse)]
```

## Pitfalls

- If the universal-size restriction leaves only 1 size (e.g., all folds
  terminate at size 1 due to very small n), revert to any-fold averaging and
  issue a warning; there is no useful CV signal in this case.
- Per-fold MSE at a given size is computed as the OLS prediction error on the
  holdout using the model fit up to that size on the training data (NOT the
  GCV-pruned model). This requires storing per-size path data (coefficients at
  each backward-pass size), which adds memory but enables cheap holdout scoring
  via a single matrix multiply per size.
- Report `size.argmin` even when `cv.1se = TRUE` so the user can see both the
  minimum-MSE size and the parsimonious 1-SE size.

---

## Rejected Entries (not proposed)

| Candidate | Reason for Rejection |
| --- | --- |
| Auto.linpreds heuristic (v0.3) | Attempted from scratch without source; known to be incorrect relative to the reference implementation (boundary-knot detection too eager). Not validated. |
| Class (e) RSS blowup guard (v0.2) | Implementation-specific sanity check (reject column if post-OLS RSS grows > 6 OOM vs intercept RSS). Useful but too narrow to generalize without the larger forward-pass context. Not novel enough. |
| Per-pair sequential array packing | Explicitly failed: pack overhead (~1 double per row) exceeded savings on sparse-parent-support cells. Entry 1 (prefix sums) supersedes this approach. |
| VarBatchScanner (batch all parents sharing a variable in one row sweep) | Tested and failed: per-row branch overhead exceeded memory savings. Documented in NEWS.md v0.12. |
| Closed-form Cholesky downdate without pivoting | Unstable on rank-deficient designs. Tested and reverted. Entry 3 (Givens) is the correct approach. |
| Backward init from forward-maintained (R, Qty) | Tested and reverted: model rebuild cost ate the QR-init savings; tightly-tied drops differ slightly, causing downstream MSE regression. The warning in Entry 2 already documents this. |
| Skip periodic Householder refresh in backward | Tested: 6% 1-thread improvement but multi-thread regression; 2 cells crossed back above reference tool. Documented as the justification for the refresh cadence in Entry 3. |
| Small-n fast path (force serial for n < 1000) | Marginal effect dominated by scheduler jitter. Not generalizable. |

---

## Duplicate Check Results

| Proposed Entry | Nearest Existing Entry | Overlap Assessment |
| --- | --- | --- |
| Entry 1: Friedman Fast-LS via Prefix Sums | None found in brain seedbank (index empty) | New |
| Entry 2: Incremental QR Maintenance | None found | New |
| Entry 3: Givens Downdate + Periodic Refresh | None found | New |
| Entry 4: AVX2 SIMD Cross-Lane Reduction | None found | New |
| Entry 5: Fast-MARS Priority Cache | None found | New |
| Entry 6: Successive Halving for Grid Search | None found | New |
| Entry 7: Subsample Warm-Start | None found | New |
| Entry 8: Shared Forward Pass Across Grid Cells | None found | New |
| Entry 9: Bootstrap Bagging for MARS UQ | None found | New |
| Entry 10: CV Aggregation Over Universal Sizes | None found | New |

Note: the brain seedbank index is empty (no prior contributions uploaded).
The five entries proposed in the earlier v0.0.0.9000 session
(`brain-contributions.md`) were noted in CLAUDE.md as saved locally but not yet
uploaded. Entries 1-10 above cover the v0.0.0.9001 through v0.0.0.9020 work
and are distinct from those earlier five entries (FMA determinism,
hand-rolled QR rationale, RcppParallel pattern, two-tier numerical-parity
acceptance, CRAN-only-ERROR triage).

---

## Privacy Scrub Verification

For each proposed entry:
- [x] No GitHub usernames, repo names, or org names
- [x] No file paths, directory structures, or package names (package name
      appears only where it is the object of discussion, not as a path or
      identifier; genericized in all examples)
- [x] No issue/PR numbers, commit SHAs, or branch names
- [x] No GitHub URLs or email addresses
- [x] All code references use generic placeholder names in examples
- [x] No dataset names, column names, or data file paths
- [x] Benchmark numbers (ratios, wall-clock) are retained as they are
      methodologically informative (not identifying); they describe the
      technique's empirical effect, not a specific user's system.
