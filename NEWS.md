# roadrunner 0.0.0.9032

## Bug fixes (statsclaw 2026-05-13 audit triage, BUG-008..BUG-013)

- **BUG-011 (correctness, medium)**: `ares.formula()` silently dropped
  `subset = ...` (it fell into `...` and went nowhere) and silently
  absorbed `offset(...)` terms as ordinary predictors. Both produced
  wrong fits with zero indication. Fix: (a) add explicit `subset` arg
  to `ares.formula` and pass it through to `model.frame` via the
  lm()-style `match.call()` construction (so NSE inside model.frame
  doesn't trip over the `subset` symbol resolving to the base R
  function); (b) detect `offset()` terms via `attr(terms, "offset")`
  before `model.frame` runs and `stop()` with an actionable message.
  Offset pass-through to the post-hoc GLM refit was punted to a later
  release (touches predict + bag + autotune compose paths). Regression
  test: `tests/testthat/test-bug-011-formula-subset-offset.R`.

- **BUG-008 (correctness, high)**: `predict()` used to return finite WRONG
  values for newdata rows containing `NA` when training used
  `na.action = "omit"`. The pre-existing warning promised "the affected
  rows will return NA predictions" but the code never imposed `NA` --
  `NaN > 0` in the C++ hinge evaluates to `FALSE`, collapsing each
  affected hinge to 0 and yielding a deterministic but wrong prediction.
  Fix: detect NA rows in `xnew` before the C++ basis pass, zero-fill the
  NaN cells, and re-impose `NA` on those rows after the linear-predictor
  compute (including bag mean, bag SE, and the `interval = "pint"` matrix
  path). Regression test:
  `tests/testthat/test-bug-008-predict-na-rows.R`.

- **BUG-009 (correctness, high)**: bagged `predict(..., type = "link")`
  for non-gaussian families used to return `g(mean(g^{-1}(eta_b)))`, i.e.
  the link applied to the response-scale bag mean. By Jensen's
  inequality this is not `mean(eta_b)`, and the simulator audit showed
  divergences of hundreds of log-odds units for binomial bags at
  moderate signal -- making `type = "link"` numerically unreliable for
  any downstream use. Fix: collect per-replicate linear predictors
  `etas` and response-scale predictions `resps` separately;
  `type = "link"` returns `rowMeans(etas)`, `type = "response"` returns
  `rowMeans(resps)` (unchanged from the prior behaviour). Bag SE is
  computed on whichever scale was returned. Regression test:
  `tests/testthat/test-bug-009-bagged-link-jensen.R`.

- **BUG-012 (correctness, medium)**: formula-path fits with derived
  terms (`I(x^2)`, `poly(x, 2)`, `log(x + 10)`, `scale(x)`,
  `splines::bs(x)`, ...) used to fit successfully but `predict()` would
  fail with "newdata is missing columns: I(x^2)" because `predict.ares`
  looked for the *expanded* column name as a literal column of newdata,
  not re-evaluating the original `terms` object on newdata. Fix: when
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
  `model.matrix(~ ., newdata)`'s default `na.action = na.omit` silently
  dropped those rows -- `length(predict(fit, newdata))` was less than
  `nrow(newdata)`. Fix: detect NA in factor / character newdata columns
  up front; error with a clear message naming the column(s), mirroring
  BUG-004's OOV path. Regression test:
  `tests/testthat/test-bug-010-predict-factor-na.R`.

# roadrunner 0.0.0.9031

## Performance

- Autotune now reuses the Householder R / Qty computed during the
  shared forward pass instead of recomputing them inside each
  `mars_backward_only_cpp` call. Saves the O(n*M^2) initial Householder
  pass at every per-cell backward replay; the cached R / Qty are
  byte-identical to the recomputed values, so selected basis,
  coefficients, and GCV are unchanged. Measured speedup on the v0.26
  speed baseline (`inst/sims/v0.26-speed-baseline.R`, 24-cell grid x
  5 reps, nthreads=4): geometric-mean 1.23x across the grid; on the
  autotune cells specifically 1.05-1.15x (gaussian highdim p=20:
  1.09x; gamma n=1000 p=10: 1.11x; binomial n=1500 p=20: 1.15x).
  Determinism invariant (`nthreads=1 == nthreads=N`) preserved.
  Internal C++ entries gain optional `compute_forward_qr` (mars_fit_cpp)
  and `R_in` / `Qty_in` (mars_backward_only_cpp) parameters; the public
  R API is unchanged.

# roadrunner 0.0.0.9030

## New features

- `plot(fit)` now produces a 4-panel diagnostic display in a 2x2 grid,
  modelled on `stats::plot.lm()`: residuals vs fitted, normal Q-Q of
  standardized residuals, scale-location, and residuals vs leverage
  with Cook's-distance contours. Panels 4 (Cook's distance) and 6
  (Cook's vs leverage) are also available via `which`. For binomial /
  poisson / gamma fits, the residual-vs-fitted panel uses deviance
  residuals and the hat matrix uses canonical-link IRLS working
  weights. The training observation weights are now stored on the
  fitted object (`$weights`).

# roadrunner 0.0.0.9029

## Bug fixes

Triage of 2026-05-11 adversarial audit (`audit-2026-05-11/`):

- BUG-001 (high): bagged GLM refits for `family = "binomial" | "poisson" | "gamma"`
  now reuse the same bootstrap indices used to select the basis. The unseeded
  path used to redraw indices from the live RNG, so the post-hoc GLM
  coefficients were fitted on a different bootstrap sample than the basis.
- BUG-002: `predict(bagged_fit)` (with `newdata = NULL`) now returns the bag
  mean, matching `predict(bagged_fit, x_train)`. The training `x` is stored
  on the fit (`out$x`) to support this.
- BUG-003: `varmod = "lm"` prediction intervals now warn and floor at a
  meaningful lower bound (rather than `1e-12`) when extrapolation makes the
  predicted MAD non-positive. In-sample PIs are unchanged.
- BUG-004: `predict()` errors loudly on out-of-vocabulary factor or character
  levels in `newdata`, instead of silently dropping the affected rows.
- BUG-005: `weights` must now be strictly positive. Zero-weight rows used to
  bias `GCV` downward and produce over-fitting; drop the rows from `x`/`y`
  instead.
- BUG-006: documented honestly that `varmod = "lm"` captures only
  yhat-dependent residual scale, not x-driven heteroscedasticity.
- BUG-007: `family = "poisson"` rejects `all(y == 0)` (degenerate GLM fit);
  all families reject constant `y` (only the intercept can be fit).

# roadrunner 0.0.0.9028

## Package

- Renamed from `ares` to `roadrunner`. The MARS fitter remains
  available as `ares()`, with the `"ares"` S3 class and methods
  unchanged. Update `library(ares)` calls to `library(roadrunner)`.

## `ares()`

- `family` accepts `"gaussian"` (default), `"binomial"`, `"poisson"`,
  and `"gamma"`. The forward + backward MARS pass runs on the
  numeric response; the selected basis is then refit with
  `stats::glm.fit()` for non-gaussian families. `predict()` gains
  `type = c("response", "link")`.
- `weights` argument for observation weights. Composes with CV
  pruning, autotune, and bagging.
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
  `predict()` averages across replicates; `se.fit = TRUE` attaches the
  per-row bag standard deviation.
- `varmod = "const" | "lm"` (gaussian only) stores a residual variance
  model at fit time, enabling `predict(interval = "pint")` for
  approximate prediction intervals.
- `na.action = c("impute", "omit")` handles missing values in `x`.
  Default `"impute"` median-imputes numeric columns and stores the
  medians for reapplication at predict time.
- Factor and character columns in a data-frame `x` are expanded via
  `model.matrix` and replayed on new data.

## Engine

- Fits are parallel-deterministic: at a fixed `seed.cv`, results are
  byte-identical across thread counts.
- Default `auto.linpreds = TRUE` and `adjust.endspan = 2L` for
  `family = "binomial"` (matches `earth`'s binomial defaults). Other
  families keep the gaussian-conservative defaults.

## Compatibility

- Requires R (>= 4.1.0).
- Depends only on `Rcpp` and `RcppParallel`. `earth`, `bench`, and
  `testthat` are Suggests only.
