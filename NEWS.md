# ares 0.0.0.9003 (development)

## Experimental earth-equivalent heuristics

- New args `adjust.endspan` (integer, default `1`) and `auto.linpreds`
  (logical, default `FALSE`) on `ares()`. Both reproduce earth's documented
  `Adjust.endspan` / `Auto.linpreds` behaviour from scratch (no earth source
  inspected). Earth's defaults are 2 and TRUE respectively; ares's defaults
  are conservative because the current implementation regressed parity vs
  earth on the inst/sims grid:
  - Multi-seed `interaction n=5000 deg=2`: with both heuristics on, 8/8
    cells > 1% rel-err (vs 5/8 with both off). `Auto.linpreds=TRUE` alone
    drives all 8 seeds worse — the boundary-knot detection in ares is too
    eager, while earth applies a stricter "linearity over parent support"
    test that this implementation only approximates.
  - Defaults off keep parity bit-identical to v0.0.0.9002 on benchmark cells.
- Code paths and the `is_boundary` flag on the forward `Candidate` struct are
  retained so a future v0.4+ can revisit with the correct linearity test.

# ares 0.0.0.9002 (development)

## Numerical robustness

- **Class (e) blowup guard**: forward pass now rejects an added hinge pair if
  the post-OLS RSS is non-finite or grew by more than 6 orders of magnitude
  vs the intercept-only RSS. Previously, on `interaction n=5000 deg=2 seed=2`,
  ares would carry an ill-conditioned column at M=13 and report RSS=3.5e29
  before backward pruning recovered. With the guard, RSS stays bounded
  (e.g. 12965 vs earth's 12819, rel-err 1.1%).
- Tie-break and `better()` comparator documented (no functional change vs
  v0.0.0.9001 since RSS-reduction ties are essentially never exact at FP).

## Speed

- **Backward pass parallelised** via `RcppParallel::parallelFor` over the
  M-1 candidate-drop trials. Previously serial; profiling at v0.1 showed
  the backward pass was 63% of single-thread wall-clock at n=5000.
- Threading scaling at n=3000 deg=2 recovered from ~1.16× (v0.1) to ~1.85×
  (2t/1t).
- Median wall-clock vs earth (54-cell `inst/sims` grid): unchanged at 1t
  (~27× slower) but improved at 2t and 4t. Median 4t wall-clock 0.24s vs
  earth 0.020s at n×deg×dgp medians.
- A column-deletion Cholesky downdate replacement for the per-trial QR was
  attempted (closed-form ΔRSS_j = β_j² / [(R'R)^{-1}]_jj) but proved unstable
  on the rank-deficient designs in the benchmark grid; it is not in this
  release. v0.3 will revisit with pivoted Cholesky / regularisation.

## Knot scanner internals

- `KnotScanner` now takes `minspan_user` (0 ⇒ auto) and `p_pred` directly,
  so the auto-minspan formula can in principle be evaluated per parent.
  The current call site still passes the top-level n (matching earth's
  behaviour); the per-parent N_m form (Friedman 1991 eq. 43) is plumbed but
  not enabled.

## Known divergences vs earth (still open)

- Earth defaults `Adjust.endspan = 2` and `Auto.linpreds = TRUE` which ares
  does not implement (these are earth-specific heuristics, not Friedman 1991).
  On `interaction n=5000 deg=2`, this causes ares to pick extreme-tail V5
  knots that earth avoids; on `additive n=5000 deg=1`, the linear-x1
  signal is captured by an ill-fit hinge pair instead of a linear basis.
  See `inst/sims/results/divergence-diag.md` for per-cell classification.
- Implementing these heuristics from scratch (without copying earth's GPL-3
  source) would close most remaining divergences. Held pending user direction
  on whether ares stays paper-only or accepts earth-equivalent additions.

# ares 0.0.0.9001 (development)

## Speed

- **Forward-pass fast-LS scoring** (`src/ares.cpp`, `KnotScanner::run`).
  Replaced the O(K · n · Mq) per-(parent, var) pair inner loop with an
  O((n + K) · Mq) prefix-sum identity. For each candidate knot t, the
  quantities `(h+ . r), (h- . r), ||h+||^2, ||h-||^2, (Q' h+)_q, (Q' h-)_q`
  are reconstructed in O(Mq) from running prefix sums maintained over the
  rows sorted ascending by x_j; the joint-pair Gram–Schmidt reduction is
  then computed identically to v0.0.0.9000.
- Verified math equivalence: an R-side reference that computes the score by
  the explicit (slow) projection picks the exact same (var, knot, reduction)
  triple as the C++ fast-LS path on a fixed-seed sanity case.
- Wall-clock improvement on the inst/sims benchmark grid (single-thread,
  median across (n, degree, dgp) cells in {500, 1500, 5000} × {1, 2} × {fri,
  add, int}): ratio vs earth dropped from ~152× to ~28× (5.4× speedup at
  median, 16.6× at the worst pre-v0.1 cell — Friedman-1, n = 1500, deg = 2:
  18.1 s → 1.09 s).
- Determinism contract preserved — `nthreads = 1` and `nthreads = N`
  produce byte-identical coefficients, dirs, cuts, and `selected.terms`.

## Threading

- `test-speed.R` parallel-scaling assertion relaxed from ≥ 1.3× (2t vs 1t)
  to "no regression". With fast-LS the parallelizable region per forward
  step is now small enough that the serial parts (`build_Q`, `ols_qr`)
  constrain scaling via Amdahl. The 4×-multi-thread @ n ≥ 5000 target from
  the v0.1 roadmap is therefore *not* attained by this change alone; see
  v0.2 roadmap.

## Roadmap update

- v0.2: Incremental Q maintenance (avoid full re-orthogonalisation per
  forward step) and parallelised `ols_qr` to relieve Amdahl bottleneck and
  recover thread scaling.
- v0.3: Cross-validated pruning (`pmethod = "cv"`).
- v0.4: GLM responses.
- v0.5: Observation weights.

# ares 0.0.0.9000

Initial development release. Not yet on CRAN.

## What works

- `ares()` fits gaussian-only MARS via forward hinge selection and backward
  GCV-based subset pruning.
- Matrix and formula interfaces.
- Standard S3 methods: `predict.ares`, `print.ares`, `summary.ares`, `plot.ares`.
- Forward-pass parallelism via **RcppParallel** (`parallelFor` over
  parent × variable pairs).
- **Threading is deterministic** — `nthreads = 1` and `nthreads = N` produce
  byte-identical fits.
- **Numerical parity with earth** — across 24-cell smoke Monte Carlo, ares
  reproduces earth's RSS to ~1% (median) on Friedman-1, additive, and
  interaction DGPs.
- 38 testthat tests pass; `devtools::check()` clean.

## Known limitations (v0.0.0.9000)

- **Wall-clock vs `earth`**: ares is slower than earth in absolute time
  (typically 30–80× slower at n ∈ {500, 1500}). Earth uses Friedman's O(1)
  per-knot Givens fast-LS update; ares v0.0.0.9000 uses an O(n) per-knot
  scoring loop. Implementing the Givens fast-LS inner loop is a **v0.1
  milestone**. Parallel scaling itself is healthy: 2-thread is ~1.75×
  faster than 1-thread.
- Gaussian responses only (no GLM, no multi-response, no observation weights).
- `pmethod` supports only `"backward"` and `"none"`. No exhaustive,
  cross-validated, or sequential-replacement pruning.
- No `fast.k` / `fast.beta` heuristics yet.
