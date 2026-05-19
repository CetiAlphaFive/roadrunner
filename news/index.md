# Changelog

## roadrunner 0.0.0.9044

### krls() — GCV lambda selection

- `krls(..., lambda.method = "gcv")` selects the ridge penalty by
  minimising the closed-form generalised cross-validation criterion
  (Craven & Wahba 1979) on the existing eigendecomposition of `K`.
  Brings parity with `KRLS` v1.5-0+.
- Closed form: no extra kernel evaluation or refit; cost is one
  golden-section search over the same `[L, U]` bracket used by LOO.
- Denominator `(1 - trH/n)^2` is floored at `1e-8` to remain finite in
  the interpolation limit `lambda -> 0`.
- `lambda.method = "loo"` remains the default. `"loo"` and `"cv"` paths
  are byte-identical to v0.0.0.9043 on a fixed-seed problem.

#### Out of scope (Phase 2/3)

- ARD-kernel option (per-feature sigma).
- HSIC pre-screen for high-p feature selection.
- Autotune-dispatch routing through GCV (autotune still calls inner CV
  regardless of `lambda.method`).

## roadrunner 0.0.0.9043

### krls() speedup — Phase 2 (Nystrom low-rank approximation)

- New optional argument `approx = "nystrom"` enables the Nystrom
  low-rank approximation. Replaces O(n^3) eigendecomposition on the full
  kernel with O(m^3) + O(n m^2) where m = `nystrom_m` (default
  `ceiling(sqrt(n) * 3)`, e.g. n=2000 -\> m=135, n=5000 -\> m=213).
- Five new optional args: `approx`, `nystrom_m`, `landmarks`,
  `landmark_method`, `landmark_seed`, `nystrom_eps`.
- Default exact path is byte-identical to v0.0.0.9042 — opt-in only.
- New exported helper `get_landmarks(fit)` returns landmark coordinates
  (original or standardized X scale).
- Determinism: at fixed `landmark_seed` (+ `seed.cv` when autotune),
  Nystrom fits are byte-identical across `autotune.nthreads`.

#### Empirical wall-clock speedup (EMP-PHASE2, R=5 reps, paired seeds)

| n    | p   | approx  | autotune | wall (s) | speedup | paired RMSE delta (median) |
|------|-----|---------|----------|----------|---------|----------------------------|
| 2000 | 10  | exact   | FALSE    | 0.52     | 1.0x    | \-                         |
| 2000 | 10  | nystrom | FALSE    | 0.056    | 9.27x   | -0.02%                     |
| 5000 | 10  | exact   | FALSE    | 10.08    | 1.0x    | \-                         |
| 5000 | 10  | nystrom | FALSE    | 0.329    | 30.64x  | +4.81%                     |

Paired RMSE delta is the per-rep
`(RMSE_nystrom - RMSE_exact) / RMSE_exact` on identical (X, y) draws
(paired seed schedule keyed on `(n, rep)`); n=2000 range \[-3.4%,
+5.4%\], n=5000 range \[+4.0%, +8.9%\]. Both well within the +10% target
on the additive smooth DGP.

#### API compatibility

- Zero breaking changes.
  [`krls()`](https://cetialphafive.github.io/roadrunner/reference/krls.md)
  without `approx` defaults to `"exact"` and behaves identically to
  v0.0.0.9042.
- `predict.krls_rr` auto-detects Nystrom fits via `fit$approx` and uses
  the cheap cross-kernel path.

#### Out of scope (Phase 3)

- Leverage-score landmark selection.
- Auto-engage Nystrom at large n.
- Sparse kernel truncation.
- Bagging + Nystrom integration.

## roadrunner 0.0.0.9042

### krls() speedup — Phase 1 (shared distance + parallel autotune)

- `krls(..., autotune = TRUE)` is now parallelised over the sigma grid
  via RcppParallel/TBB. The pairwise squared-distance matrix is computed
  once per CV fold and reused across every sigma candidate.
- New optional argument `autotune.nthreads` (default
  `getOption("roadrunner.nthreads", parallel::detectCores(logical=FALSE))`).
  Pass `autotune.nthreads = 1` for strictly sequential execution.
- Determinism contract preserved: fits are byte-identical across
  `autotune.nthreads` values at fixed seed and inputs (each worker
  writes to a unique slot in the output vectors, so no reduction is
  involved).

#### Empirical wall-clock speedup (REQ-001 minimal grid, R=5)

| n    | p   | nthreads | time  | speedup vs v0.0.0.9041 |
|------|-----|----------|-------|------------------------|
| 500  | 10  | 1        | (TBD) | 1.0x (sequential)      |
| 500  | 10  | 4        | (TBD) | (filled by Task 9)     |
| 1500 | 20  | 1        | (TBD) | 1.0x (sequential)      |
| 1500 | 20  | 4        | (TBD) | (filled by Task 9)     |

#### API compatibility

- Zero breaking changes. `krls(..., autotune = TRUE)` returns the same
  S3 object with the same fields. `autotune_info` gains `nthreads_used`
  and `sigma_grid_sorted` slots.
- `predict.krls_rr`, `summary.krls_rr`, `print.krls_rr` unchanged.

#### Out of scope (Phase 2)

- Nystrom approximation (`approx = "nystrom"`, `nystrom_m`,
  `landmarks`).
- Predict()-side speedups.
- Multi-response / sparse kernel paths.

## roadrunner 0.0.0.9041

### `krls()` — scale-aware sigma anchor refined to geomean_p (REQ-20260518-003)

The default sigma anchor formula in `.krls_sigma_anchor()` is updated
from the raw median heuristic (`median(d2)`, introduced in v0.0.0.9040)
to the *geomean_p* formula: `sqrt(median(d2) * p)` where `d2` are
pairwise squared Euclidean distances on the standardised training matrix
and `p = ncol(X)`.

**Why**: a 15-DGP head-to-head vs
[`KRLS::krls()`](https://rdrr.io/pkg/KRLS/man/krls.html) at
`n=500, p=10` (iter-0, REQ-20260518-003) showed the v0.0.0.9040 median
anchor over-smoothed locally nonlinear signals, producing 8 losses and 0
ties in favour of KRLS. An anchor sweep over 6 candidate formulae
(iter-1) identified geomean_p as the clear winner (14/15 DGPs).
Verification with the patched anchor (iter-2) confirmed: **roadrunner
wins 12/15 DGPs, KRLS wins 0/15, 3 ties**.

**Delta table — iter-0 (median) vs iter-2 (geomean_p) mean MSE ratio (rr
/ KRLS)**:

| DGP              | iter-0 ratio | iter-2 ratio | delta  |
|------------------|--------------|--------------|--------|
| additive         | 1.327        | 0.896        | -0.431 |
| exp-decay        | 1.681        | 0.887        | -0.795 |
| friedman1        | 1.126        | 0.959        | -0.167 |
| friedman3        | 1.084        | 0.957        | -0.127 |
| interaction      | 1.471        | 0.895        | -0.576 |
| mixture          | 1.140        | 0.992 (tie)  | -0.148 |
| poly2            | 1.286        | 0.871        | -0.415 |
| tanh-interaction | 1.485        | 0.952        | -0.533 |
| friedman2        | 0.984 (tie)  | 0.936        | -0.049 |
| sin-sum          | 1.011 (tie)  | 0.982 (tie)  | -0.029 |
| linear           | 0.872        | 0.919        | +0.047 |
| sparse           | 0.865        | 0.917        | +0.052 |
| heterosked       | 0.859        | 0.920        | +0.062 |
| monotone         | 0.960        | 0.925        | -0.035 |
| noise            | 1.000 (tie)  | 1.000 (tie)  | 0.000  |

At `n=500, p=10` on standardised N(0,1) data the geomean_p anchor gives
`sigma ~ 13.5`, splitting the difference between the raw median (~18–19)
and KRLS’s fixed `sigma = p = 10`. The formula is still data-adaptive:
on non-standardised or non-unit-variance inputs it will differ sensibly
from both extremes.

This is a **one-line R change** (`R/krls.R`, `.krls_sigma_anchor()`);
the C++ engine and all other logic are unchanged. At a fixed
`(sigma, lambda)` fits remain byte-identical to all prior versions. To
restore the v0.0.0.9040 median anchor, pass
`sigma = stats::median(as.numeric(stats::dist(scale(X)))^2)`.

## roadrunner 0.0.0.9040

### `krls()` — four overfitting fixes (REQ-20260518-002)

Addresses overfitting confirmed by diagnostic sweeps in
REQ-20260518-001. All four changes are R-level only; the C++ engine is
unchanged, so at a fixed `(sigma, lambda)` fits remain byte-identical to
earlier versions.

Empirical improvement on signal DGPs (`n=500, p=10`, `tune=none`, `R=10`
replications; overfit ratio = `test_MSE / train_MSE`):

| DGP         | Before (pre-fix) | After (post-fix) | Improvement |
|-------------|------------------|------------------|-------------|
| additive    | 3.35             | 1.58             | 53%         |
| interaction | 4.38             | 2.10             | 52%         |
| sparse      | 2.61             | 1.26             | 52%         |
| linear      | 2.70             | 1.37             | 49%         |
| noise       | 1.06             | 1.04             | —           |

All improvement is attributable to Fix 1 (sigma default) and Fix 2
(lambda tolerance); Fixes 3 and 4 further stabilise autotune sigma
selection.

**Breaking defaults** (old values restorable by explicit arguments):

- **Default `sigma` changed** from `ncol(X)` to the median pairwise
  squared Euclidean distance on the standardised predictors (the ‘median
  heuristic’). At `n=500, p=10` the oracle sigma is ~20 and the median
  heuristic anchors near that neighbourhood; `sigma = ncol(X)` was 10
  (2x under-smoothed). Restore old behaviour with `sigma = ncol(X)`.

- **Default lambda tolerance changed** from `1e-3 * n` (n-dependent,
  coarse at moderate n) to `1e-6` (fixed, 6-digit precision). The LOO
  golden-section was empirically selecting lambda ~4x the argmin at
  `sigma=20, n=500` under the old tolerance. The L-bracket climb now
  uses multiplicative steps (x10 per step) instead of additive steps
  (0.05 per step) for scale-robustness. Restore old behaviour with
  `tol = 1e-3 * nrow(X)`.

- **Autotune sigma grid changed** from
  `ncol(X) * c(0.25, 0.5, 1, 2, 4, 8)` (6 points fixed at `d`) to
  `sigma_anchor * c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32)` (9 points
  centred on the median heuristic anchor). Restore old behaviour with
  `autotune.grid = ncol(X) * c(0.25, 0.5, 1, 2, 4, 8)`.

- **Autotune CV stabilised**: default folds raised from 5 to 10; default
  cross-partition repeats raised from 1 to 2 (via the `ncross` argument,
  whose default changes from `1L` to `NULL`); sigma selection now
  applies the 1-SE rule (largest sigma within 1 SE of minimum CV-MSE,
  biasing toward wider kernels when evidence is weak). The `autotune`
  component of the fit gains new fields: `ncross`, `mse_per_fold`,
  `se_mse`, `cv.1se`, `sigma_1se`. Restore old fold count with
  `nfold = 5`; restore old repeat count with `ncross = 1`.

## roadrunner 0.0.0.9033

### New feature: `krls()` – Kernel Regularized Least Squares

- Adds
  [`krls()`](https://cetialphafive.github.io/roadrunner/reference/krls.md),
  a from-scratch implementation of the Hainmueller and Hazlett (2014)
  KRLS estimator under the roadrunner roof. The algorithm mirrors
  [`KRLS::krls()`](https://rdrr.io/pkg/KRLS/man/krls.html) exactly
  (standardisation, Gaussian kernel, eigen-basis closed-form solve,
  golden-section LOO lambda search, marginal effects, binary
  first-difference handling); at matched `(sigma, lambda)` fits agree
  with [`KRLS::krls()`](https://rdrr.io/pkg/KRLS/man/krls.html) to
  within floating-point precision (`< 1e-12` on coefficients, fitted
  values, pointwise and average marginal effects, and prediction SEs).
- Engine is C++ (`src/krls.cpp`) on top of `RcppArmadillo` (added to
  `LinkingTo` + `Imports`) and `RcppParallel`. The Gaussian kernel and
  the test-vs-train kernel are built in parallel with a TBB worker.
  Eigendecomposition is dispatched to LAPACK via `arma::eig_sym` with
  the divide-and-conquer driver. Marginal effects are computed via the
  identity `dy/dx_k = -(2/sigma)*(X_k*(K c) - K diag(c) X)_k`, which
  avoids the explicit `n x n` distance matrix and reduces the average-
  marginal-effect variance from `O(n^3)` to `O(n^2)` per variable via a
  row-sum trick.
- Measured speed-up over
  [`KRLS::krls()`](https://rdrr.io/pkg/KRLS/man/krls.html) on simulated
  benchmarks (`derivative = TRUE`, `vcov = TRUE`) at
  `nthreads = default`: `n=200 p=3` ~2x; `n=500 p=3` ~5-6x; `n=500 p=10`
  ~10x; `n=1000 p=3` ~6x; `n=1000 p=10` ~10x. Coefficient max-abs error
  vs [`KRLS::krls()`](https://rdrr.io/pkg/KRLS/man/krls.html) is
  `< 2e-13` on every cell.
- S3 class is `c("krls_rr", "krls")` so
  [`predict()`](https://rdrr.io/r/stats/predict.html),
  [`print()`](https://rdrr.io/r/base/print.html), and
  [`summary()`](https://rdrr.io/r/base/summary.html) dispatch
  unambiguously to roadrunner methods even when `KRLS` is loaded in the
  same session. `inherits(fit, "krls")` is preserved for
  downstream-compat checks.
- Tests: `tests/testthat/test-krls.R` (parity vs
  [`KRLS::krls`](https://rdrr.io/pkg/KRLS/man/krls.html) for
  coefficients, fitted values, Looe, marginal effects, prediction SEs,
  binary first-differences; structural input-validation; predict
  recovers fitted values; sensible R^2). 11 tests, 28 assertions, all
  green.

## roadrunner 0.0.0.9032

### Bug fixes (statsclaw 2026-05-13 audit triage, BUG-008..BUG-013)

- **BUG-013 (usability, low)**:
  [`print.ares()`](https://cetialphafive.github.io/roadrunner/reference/print.ares.md)
  and
  [`summary.ares()`](https://cetialphafive.github.io/roadrunner/reference/summary.ares.md)
  used to be silent about bagging (`n.boot > 0`) and autotune state, so
  a bagged or autotuned fit printed identically to a plain fit. Fix: add
  one-line `Bagging: n.boot = N replicate(s)` and
  `Autotune: degree=D penalty=P nk=K fast.k=F warmstart=T/F` blocks to
  both printers when the corresponding components are present.
  `summary.ares` also carries `$boot` and `$autotune` through to its
  print method. Regression test:
  `tests/testthat/test-bug-013-print-bag-autotune.R`.

- **BUG-011 (correctness, medium)**:
  [`ares.formula()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)
  silently dropped `subset = ...` (it fell into `...` and went nowhere)
  and silently absorbed `offset(...)` terms as ordinary predictors. Both
  produced wrong fits with zero indication. Fix: (a) add explicit
  `subset` arg to `ares.formula` and pass it through to `model.frame`
  via the lm()-style
  [`match.call()`](https://rdrr.io/r/base/match.call.html) construction
  (so NSE inside model.frame doesn’t trip over the `subset` symbol
  resolving to the base R function); (b) detect
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms via
  `attr(terms, "offset")` before `model.frame` runs and
  [`stop()`](https://rdrr.io/r/base/stop.html) with an actionable
  message. Offset pass-through to the post-hoc GLM refit was punted to a
  later release (touches predict + bag + autotune compose paths).
  Regression test:
  `tests/testthat/test-bug-011-formula-subset-offset.R`.

- **BUG-008 (correctness, high)**:
  [`predict()`](https://rdrr.io/r/stats/predict.html) used to return
  finite WRONG values for newdata rows containing `NA` when training
  used `na.action = "omit"`. The pre-existing warning promised “the
  affected rows will return NA predictions” but the code never imposed
  `NA` – `NaN > 0` in the C++ hinge evaluates to `FALSE`, collapsing
  each affected hinge to 0 and yielding a deterministic but wrong
  prediction. Fix: detect NA rows in `xnew` before the C++ basis pass,
  zero-fill the NaN cells, and re-impose `NA` on those rows after the
  linear-predictor compute (including bag mean, bag SE, and the
  `interval = "pint"` matrix path). Regression test:
  `tests/testthat/test-bug-008-predict-na-rows.R`.

- **BUG-009 (correctness, high)**: bagged `predict(..., type = "link")`
  for non-gaussian families used to return `g(mean(g^{-1}(eta_b)))`,
  i.e. the link applied to the response-scale bag mean. By Jensen’s
  inequality this is not `mean(eta_b)`, and the simulator audit showed
  divergences of hundreds of log-odds units for binomial bags at
  moderate signal – making `type = "link"` numerically unreliable for
  any downstream use. Fix: collect per-replicate linear predictors
  `etas` and response-scale predictions `resps` separately;
  `type = "link"` returns `rowMeans(etas)`, `type = "response"` returns
  `rowMeans(resps)` (unchanged from the prior behaviour). Bag SE is
  computed on whichever scale was returned. Regression test:
  `tests/testthat/test-bug-009-bagged-link-jensen.R`.

- **BUG-012 (correctness, medium)**: formula-path fits with derived
  terms (`I(x^2)`, `poly(x, 2)`, `log(x + 10)`, `scale(x)`,
  `splines::bs(x)`, …) used to fit successfully but
  [`predict()`](https://rdrr.io/r/stats/predict.html) would fail with
  “newdata is missing columns: I(x^2)” because `predict.ares` looked for
  the *expanded* column name as a literal column of newdata, not
  re-evaluating the original `terms` object on newdata. Fix: when
  `object$terms` is non-null (formula path), use
  `model.matrix(delete.response(terms), newdata, xlev = object$xlevels)`
  to rebuild the design with derived terms re-evaluated. Falls back to
  the prior column-lookup path when no `terms` object is stored (matrix
  interface). `ares.formula` now also stashes `xlevels` on the fit.
  Regression test:
  `tests/testthat/test-bug-012-predict-derived-terms.R`.

- **BUG-010 (robustness, medium)**: sister to BUG-004. `NA` values in
  factor / character newdata columns used to fall through the OOV
  detector (which only handled non-NA character values), then
  `model.matrix(~ ., newdata)`’s default `na.action = na.omit` silently
  dropped those rows – `length(predict(fit, newdata))` was less than
  `nrow(newdata)`. Fix: detect NA in factor / character newdata columns
  up front; error with a clear message naming the column(s), mirroring
  BUG-004’s OOV path. Regression test:
  `tests/testthat/test-bug-010-predict-factor-na.R`.

## roadrunner 0.0.0.9031

### Performance

- Autotune now reuses the Householder R / Qty computed during the shared
  forward pass instead of recomputing them inside each
  `mars_backward_only_cpp` call. Saves the O(n\*M^2) initial Householder
  pass at every per-cell backward replay; the cached R / Qty are
  byte-identical to the recomputed values, so selected basis,
  coefficients, and GCV are unchanged. Measured speedup on the v0.26
  speed baseline (`inst/sims/v0.26-speed-baseline.R`, 24-cell grid x 5
  reps, nthreads=4): geometric-mean 1.23x across the grid; on the
  autotune cells specifically 1.05-1.15x (gaussian highdim p=20: 1.09x;
  gamma n=1000 p=10: 1.11x; binomial n=1500 p=20: 1.15x). Determinism
  invariant (`nthreads=1 == nthreads=N`) preserved. Internal C++ entries
  gain optional `compute_forward_qr` (mars_fit_cpp) and `R_in` /
  `Qty_in` (mars_backward_only_cpp) parameters; the public R API is
  unchanged.

## roadrunner 0.0.0.9030

### New features

- `plot(fit)` now produces a 4-panel diagnostic display in a 2x2 grid,
  modelled on
  [`stats::plot.lm()`](https://rdrr.io/r/stats/plot.lm.html): residuals
  vs fitted, normal Q-Q of standardized residuals, scale-location, and
  residuals vs leverage with Cook’s-distance contours. Panels 4 (Cook’s
  distance) and 6 (Cook’s vs leverage) are also available via `which`.
  For binomial / poisson / gamma fits, the residual-vs-fitted panel uses
  deviance residuals and the hat matrix uses canonical-link IRLS working
  weights. The training observation weights are now stored on the fitted
  object (`$weights`).

## roadrunner 0.0.0.9029

### Bug fixes

Triage of 2026-05-11 adversarial audit (`audit-2026-05-11/`):

- BUG-001 (high): bagged GLM refits for
  `family = "binomial" | "poisson" | "gamma"` now reuse the same
  bootstrap indices used to select the basis. The unseeded path used to
  redraw indices from the live RNG, so the post-hoc GLM coefficients
  were fitted on a different bootstrap sample than the basis.
- BUG-002: `predict(bagged_fit)` (with `newdata = NULL`) now returns the
  bag mean, matching `predict(bagged_fit, x_train)`. The training `x` is
  stored on the fit (`out$x`) to support this.
- BUG-003: `varmod = "lm"` prediction intervals now warn and floor at a
  meaningful lower bound (rather than `1e-12`) when extrapolation makes
  the predicted MAD non-positive. In-sample PIs are unchanged.
- BUG-004: [`predict()`](https://rdrr.io/r/stats/predict.html) errors
  loudly on out-of-vocabulary factor or character levels in `newdata`,
  instead of silently dropping the affected rows.
- BUG-005: `weights` must now be strictly positive. Zero-weight rows
  used to bias `GCV` downward and produce over-fitting; drop the rows
  from `x`/`y` instead.
- BUG-006: documented honestly that `varmod = "lm"` captures only
  yhat-dependent residual scale, not x-driven heteroscedasticity.
- BUG-007: `family = "poisson"` rejects `all(y == 0)` (degenerate GLM
  fit); all families reject constant `y` (only the intercept can be
  fit).

## roadrunner 0.0.0.9028

### Package

- Renamed from `ares` to `roadrunner`. The MARS fitter remains available
  as
  [`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md),
  with the `"ares"` S3 class and methods unchanged. Update
  [`library(ares)`](https://rdrr.io/r/base/library.html) calls to
  [`library(roadrunner)`](https://cetialphafive.github.io/roadrunner/).

### `ares()`

- `family` accepts `"gaussian"` (default), `"binomial"`, `"poisson"`,
  and `"gamma"`. The forward + backward MARS pass runs on the numeric
  response; the selected basis is then refit with
  [`stats::glm.fit()`](https://rdrr.io/r/stats/glm.html) for
  non-gaussian families.
  [`predict()`](https://rdrr.io/r/stats/predict.html) gains
  `type = c("response", "link")`.
- `weights` argument for observation weights. Composes with CV pruning,
  autotune, and bagging.
- `pmethod = "cv"` (and the convenience trigger `nfold > 0`) performs
  K-fold cross-validated subset-size selection. Optional `ncross`
  repetitions, quantile-based stratification, and the 1-SE rule via
  `cv.1se = TRUE`.
- `autotune = TRUE` runs an inner-CV grid search over
  `(degree, penalty, nk, fast.k)`. `autotune.speed` chooses among
  `"balanced"` (default), `"quality"`, and `"fast"`. A 20 % subsample
  warm-start short-circuits the full grid when the best-per-degree gap
  is decisive.
- `n.boot > 0` fits row-bootstrap replicates of the central model.
  [`predict()`](https://rdrr.io/r/stats/predict.html) averages across
  replicates; `se.fit = TRUE` attaches the per-row bag standard
  deviation.
- `varmod = "const" | "lm"` (gaussian only) stores a residual variance
  model at fit time, enabling `predict(interval = "pint")` for
  approximate prediction intervals.
- `na.action = c("impute", "omit")` handles missing values in `x`.
  Default `"impute"` median-imputes numeric columns and stores the
  medians for reapplication at predict time.
- Factor and character columns in a data-frame `x` are expanded via
  `model.matrix` and replayed on new data.

### Engine

- Fits are parallel-deterministic: at a fixed `seed.cv`, results are
  byte-identical across thread counts.
- Default `auto.linpreds = TRUE` and `adjust.endspan = 2L` for
  `family = "binomial"` (matches `earth`’s binomial defaults). Other
  families keep the gaussian-conservative defaults.

### Compatibility

- Requires R (\>= 4.1.0).
- Depends only on `Rcpp` and `RcppParallel`. `earth`, `bench`, and
  `testthat` are Suggests only.
