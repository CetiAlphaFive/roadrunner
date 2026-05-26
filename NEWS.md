# roadrunner 0.0.0.9059

## `bgam()` — component-wise P-spline gradient boosting

`bgam()` fits a smooth, strictly additive model via component-wise functional
gradient boosting (Bühlmann & Yu 2003; Schmid & Hothorn 2008). At each of
`mstop` iterations the algorithm selects the single penalised B-spline
base-learner — one per predictor — that most reduces the current loss,
adds a shrinkage-scaled ridge update for that component only, and repeats.
The result is a smooth GAM with built-in variable selection: predictors that
are never selected retain a zero coefficient vector and their selection
frequency (`$selection_frequency`) serves as a natural importance metric.
`bgam()` fills the smooth-additive-model-with-variable-selection gap in the
roadrunner lineup that neither `ares()` (piecewise-linear, interaction-capable)
nor `krls()` (global kernel, no variable selection) covers natively. It scales
as O(n · p · K) per iteration and is particularly effective at n = 1000+, p
moderately large, smooth signal settings.

**Supported families:** `"gaussian"` (squared-error loss) and `"binomial"`
(logistic / binomial deviance, binary response).

**Key arguments:** `mstop` (boosting iterations), `nu` (shrinkage, default
`0.1`), `nknots` (interior knots per predictor, default `20`), `degree`
(B-spline degree, default `3`), `dpen` (penalty order, default `2`),
`df_target` (target effective df per base-learner, default `4`),
`unpenalized` (predictors entered as linear instead of spline),
`weights` (per-observation weights), `n.boot` (bagging replicates),
`autotune` (k-fold CV mstop selection, default `TRUE`).

**`meep()` integration:** `bgam` is an opt-in learner; activate via
`meep(X, y, treatment = D, learners = c("ares", "krls", "bgam"))`. Pass
extra arguments through `bgam_args = list(...)`. `bgam` is family-agnostic:
gaussian path for continuous nuisances, logistic path for propensity scores.
Hyperparameters frozen by `tune = "once"`: `mstop`, `nu`, `nknots`, `degree`,
`dpen`.

**Parity oracle:** `mboost::gamboost()` (Hofner et al. 2014). `mboost` is
listed under `Suggests` and is never a hard dependency; parity tests guard with
`skip_if_not_installed("mboost")`. Loss trajectories agree to `cor > 0.99`
at matched mstop and nu; eta-scale predictions agree to `cor > 0.85` for
gaussian and `cor > 0.80` for binomial (implementations diverge on knot
placement and lambda calibration — see test-spec for rationale).

**Trimmed simulation** (`inst/sims/results/bgam-trimmed-0.0.0.9059.md`):
4-cell × 50-rep Gaussian DGP (`y = sin(2πx1) + 0.5·x2² + N(0, 0.5²)`) shows
bgam tracks ares within Monte Carlo SE on the small-n / high-p cell where the
additive smooth signal is best exploited, and both learners reduce
misspecified OLS RMSE by 30–60%. Full 16-cell × 200-rep run (with
heteroskedastic and binomial DGPs) is deferred to a follow-up.

### Documentation

New man pages added: `bgam.Rd`, `predict.bgam.Rd`, `print.bgam.Rd`,
`summary.bgam.Rd`, `plot.bgam.Rd`.

A dedicated `vignettes/bgam.Rmd` is deferred; the plan is to add it alongside
the full simulation study in the follow-up release.

# roadrunner 0.0.0.9058

## `meep()` krls learner: curve-fitter defaults

`meep()` now invokes its `krls` learner with `vcov = FALSE` and
`derivative = FALSE` by default. The meep ensemble only consumes the
out-of-fold prediction matrix from each nuisance fit; the SE block
(coefficient + fitted-value covariance) and the marginal-effects block
(pointwise + average derivatives + variance of average derivatives)
were computed on every fold of every nuisance and then discarded.

Per-fit cost share at n = 1500 (v9057 single-fit profile):
  vcov block       ~15%  (skipped under vcov = FALSE)
  derivative block ~23%  (skipped under derivative = FALSE)

End-to-end measurement (`meep(learners = "krls", folds = 5)`, 12 threads,
Friedman-1-like DGP with binary treatment):

  n = 800   :  9.5s -> 9.1s   1.04x   (4% saved)
  n = 1500  : 53.6s -> 51.7s  1.04x   (4% saved)

OOF prediction matrices are bit-identical between the two paths
(`max|delta| == 0` on both `y_hat_oof` and `d_hat_oof`); the overall
meep speedup is modest because meep wall is autotune-dominated
(per-cell eig in the inner CV), but the saved work is real, the krls
fits in `fit$folds[[k]]$models$krls` are now `~2 * n^2 + n * p` doubles
lighter, and any downstream uncertainty quantification for DML / CF
nuisance fits should be derived from the residuals / OOF predictions
rather than the individual fit-level SEs.

Opt back in for vcov + derivative via `krls_args`:

    meep(X, y, treatment = D, learners = "krls",
         krls_args = list(vcov = TRUE, derivative = TRUE))

The `ares` / `ols` / `logreg` learner specs are unchanged.

132 meep tests pass; lint clean.

# roadrunner 0.0.0.9057

## `krls()` single-fit: closed-form spectral vcov

The unweighted (`weights = NULL`) `krls(..., vcov = TRUE)` path now
computes coefficient and fitted-value covariances directly from the
eigen-decomposition that was already produced for the lambda search,
using the spectral identity `K = V diag(d) V^T`:

    vcovmatc    = V diag(w)        V^T = tcrossprod(V scaled by sqrt(w))
    vcovmatyhat = V diag(d^2 * w)  V^T = tcrossprod(V scaled by d sqrt(w))

where `w = sigma2 / (d + lambda)^2`. This replaces two R-level
`O(n^3)` dgemm calls (`V %*% (t(V) * w)` and `crossprod(K, vcovmatc %*% K)`)
with two `O(n^3 / 2)` dsyrk calls via `tcrossprod`; no K matrix is
touched on this path. Outputs match the v9056 formulation to within
~1e-15 (machine epsilon) and are symmetric by construction.

Measurements (16-core MKL, Friedman-1, 12 threads,
`inst/sims/v0.26-krls-speed-baseline.R` cells, cumulative vs v9055):

| cell              | v9055   | v9057   | speedup |
| ----------------- | ------- | ------- | ------- |
| n=400  gcv        | 0.033s  | 0.026s  | 1.27x   |
| n=1500 gcv        | 0.403s  | 0.353s  | 1.14x   |
| n=1500 default    | 0.534s  | 0.498s  | 1.07x   |
| n=400  autotune   | 1.456s  | 0.758s  | 1.92x   |
| n=800  autotune   | 8.235s  | 5.526s  | 1.49x   |
| n=1500 autotune   | 50.87s  | 37.92s  | 1.34x   |

Geometric-mean speedup across all 9 baseline cells: **+20.7%** vs v9055.

The weighted (`!is.null(weights)`) path retains the v9056 two-gemm
formulation because the spectral simplification assumes `V` diagonalises
the kernel used in the K-times-K reduction; in the weighted path `V`
diagonalises `K_w = D K D` (not `K`), and the residual `V^T D^-1 V`
terms do not collapse. The non-Gaussian / logistic / Nystrom paths are
unaffected.

`AUDIT-D` thread-determinism is preserved (`nthreads = 1` and
`nthreads = 4` produce bit-identical `K`, `coeffs`, `fitted`, `vcov.c`,
`vcov.fitted`, `avgderivatives`). 639 krls tests pass with
`ARES_FULL_TESTS = 1`; lint clean.

# roadrunner 0.0.0.9056

## `krls()` autotune: progressive sigma-grid pruning

`krls(..., autotune = TRUE)` now drops clearly-losing sigmas from the
inner cross-validation as soon as evidence accumulates. After
`prune_warmup` folds (default 2), any sigma whose mean held-out MSE
exceeds `prune_factor` x the running best (default 1.5) is removed from
the remaining fold evaluations, unless it sits within
`prune_keep_neighbors` grid positions of the running best (default 1,
prevents collateral damage near the minimum on flat MSE landscapes).

Measured on a 16-core MKL box (Friedman-1, `autotune.speed = "fast"`,
12 worker threads), comparing v0.0.0.9056 vs v0.0.0.9055 over the
`inst/sims/v0.26-krls-speed-baseline.R` grid:

| cell              | v9055   | v9056   | speedup |
| ----------------- | ------- | ------- | ------- |
| n=400  autotune   | 1.46s   | 0.80s   | 1.82x   |
| n=800  autotune   | 8.23s   | 5.50s   | 1.50x   |
| n=1500 autotune   | 50.87s  | 38.23s  | 1.33x   |

Selected `sigma`, `lambda`, `coeffs`, and held-out RMSE are bit-identical
to v9055 on every cell tested. The pruning gate is `use_inner_cpp`-only
(Gaussian + ls + unweighted + no user bracket), so the legacy serial
fallback and non-Gaussian / weighted / logistic paths are unchanged.

New options (all with safe defaults):

* `roadrunner.krls.autotune.prune_factor` (default `1.5`). Set to `Inf`
  to disable pruning entirely and recover v9055 behaviour byte-for-byte.
* `roadrunner.krls.autotune.prune_warmup` (default `2L`). Minimum folds
  populated before pruning may activate.
* `roadrunner.krls.autotune.prune_keep_neighbors` (default `1L`).
  Window around the running-best sigma that is protected from pruning.

`fit$autotune` gains `prune_factor`, `prune_warmup`, `prune_keep_neighbors`,
`prune_active`, `pruned_n`, and `prune_log` for inspection.

## New: KRLS speed-tracking sim

`inst/sims/v0.26-krls-speed-baseline.R` patterned on the existing ares
baseline; cells: `(n, p) x setup`, `setup in {default, gcv, autotune_fast}`.
Captures wall, OOS RMSE, and a fit digest as the byte-identity oracle.

# roadrunner 0.0.0.9055

## Tier A audit bug fixes (8 fixes)

This release closes eight bugs identified by the 2026-05-22 audit. Full
plan: `docs/superpowers/audits/2026-05-22-package-audit/01-tier-A-fix-plan.md`.

* **BUG-014** (`ares`): `pmethod = "cv"` crashed with
  `"subscript out of bounds"` when a small fold's forward pass collapsed
  to an intercept-only model (`M_full == 1`). The C++ engine populates
  `$path.subsets` / `$path.coefs` only when `M > 1`, so those lists were
  empty for the M=1 case. The CV loop now guards the index against list
  length and scores the intercept-only prediction directly from
  `$coefficients`, letting the fold contribute a valid size-1 MSE.
* **BUG-015** (`predict.ares`): when the fit carried `$factor_info` (or
  a `$terms`), `model.matrix()` ran through `model.frame()` with the
  global `na.action` (default `na.omit`), silently dropping `newdata`
  rows that had `NA` in numeric columns. The returned vector was shorter
  than `nrow(newdata)`, violating the predict-length contract. Both
  branches now build the model frame with `na.action = na.pass` so NAs
  flow through to the downstream median-impute / NA-row-mask path.
* **BUG-016** (`krls`): `predict(weighted_krls_fit, type = "variance")`
  returned ~0 across all rows. The fit stashed `(V, dvals)` from the
  *weighted* kernel `K_w = D K D` (`D = diag(sqrt(w))`), but the
  predict-side helper built `K_star` from the *unweighted* kernel. The
  helper now column-scales `K_star` by `sqrt(w_tr)` before the
  quad-form, implementing
  `V[f*] = K**(unwt) - K* D (K_w + lam I)^-1 D K*'`.
* **BUG-017** (`krls`): `predict(bagged_krls_fit, newdata = NULL)`
  returned `$fitted` -- the central fit -- instead of the bag mean
  across replicates. The non-NULL branch already does the bagging, so
  the NULL-newdata fast-path was inconsistent with the predict-on-newdata
  behaviour. The fix routes bagged NULL-newdata calls through the
  non-NULL branch (with `newdata = object$X`).
* **BUG-018** (`meep`): `predict.meep(nuisance = "mu0" | "mu1")` returned
  all-NA for any weighted `meep()` fit. The arm-restricted refit inside
  `predict.meep` passed `object$weights` (length n) without slicing to
  the treatment-arm subset, so the base learner errored
  `length(weights) != length(y)`, the error was swallowed, and the
  prediction collapsed to NA. Weights are now sliced to the arm rows.
* **BUG-019** (`meep`): a factor or character outcome with three or
  more levels was silently coerced to integer codes and fit as a
  regression. Multi-class outcomes are out of scope; `meep()` now
  rejects them upfront with a clear message, mirroring the
  treatment-validation pattern.
* **BUG-020** (`logreg`): `null.deviance` disagreed with
  `stats::glm()` for no-intercept formulas (`y ~ 0 + x`). The C++ engine
  always used the weighted mean of y as the null mu, but the canonical
  null model under no intercept is `eta = 0 -> mu = 0.5`. The engine now
  takes a `has_intercept` flag (R-side detected as `"(Intercept)" %in%
  colnames(X)` or `intercept = TRUE`), and uses `mu0 = 0.5` when false.
* **BUG-021** (`plda`): when `nfold > min(table(y))`, the stratified
  fold builder placed each class only into folds `1..length(class)`, so
  some folds received zero test observations and biased the CV error
  toward 0 (selecting the wrong lambda). `plda()` now rejects upfront
  with a clear message recommending the largest valid nfold.

# roadrunner 0.0.0.9054

## meep() — ols() and logreg() are now default learners

`meep()`'s default `learners` is now `c("ares", "krls", "ols", "logreg")`
(previously `c("ares", "krls")`). The two linear learners are applied
*per nuisance by family*: `ols` is fitted only for gaussian (continuous)
nuisances and `logreg` only for binomial (binary) ones, while `ares` and
`krls` remain family-agnostic and apply everywhere.

* Because the Y-model and the D-model can differ in family, applicability
  is resolved per nuisance. A gaussian outcome with a binary treatment
  fits `ols` for `outcome`/`mu0`/`mu1` and `logreg` for `treatment`.
* A learner that does not apply to a nuisance is skipped cleanly: its OOF
  column stays all-`NA` (the ensemble already excludes all-`NA` columns)
  and it is *not* recorded as a fold failure.
* Narrowing `learners` to names that are all family-incompatible with a
  nuisance (for example `learners = "ols"` with a binary outcome) is now
  an error naming the nuisance, its family, and the cause.
* Custom list-learners are never family-filtered.

**This changes default `meep()` output**: the default OOF prediction
matrices now carry `ols`/`logreg` columns and the ensemble weights are
recomputed over the larger learner set.

# roadrunner 0.0.0.9053

## ols() and logreg() — unregularized linear fitters

Two new user-facing functions, `ols()` and `logreg()`: fast, low-dependency
unregularized linear fitters with C++ (Rcpp/Armadillo) engines, classical
and HC robust standard errors, and the package's standard base-R-style API.

* **`ols(x, y, weights, n.boot, varmod, ...)`** — ordinary and weighted
  least squares via a weighted economy-QR solve. Matrix and formula
  interfaces. Returns an S3 object of class `"ols"`. A rank-deficient
  design matrix is rejected with a clear error rather than silently
  dropping columns.
* **`logreg(x, y, weights, n.boot, ...)`** — binary logistic regression by
  iteratively reweighted least squares (Fisher scoring). Matrix and
  formula interfaces; accepts a 0/1 numeric, logical, or two-level-factor
  response. Returns an S3 object of class `"logreg"`. Non-convergence (the
  signature of quasi-complete separation) is reported with a warning.
* **Robust SEs** — classical and HC0–HC3 sandwich covariances; surfaced in
  `summary()` and `predict()` via a `robust` argument.
* **Bagging** — `n.boot > 0` refits on bootstrap-weight resamples; bagged
  `predict()` returns the bag mean. Serial loop, deterministic across
  thread counts.
* **Prediction intervals** — `ols()` supports `varmod = "const" | "lm"`
  for heteroskedastic gaussian prediction intervals (shared `ares()`
  variance-model helper).
* **S3 methods** — `predict`, `print`, `summary`, `plot` for both classes.
* **`meep()` integration** — `"ols"` and `"logreg"` are opt-in `meep()`
  learners (`ols` = gaussian, `logreg` = binomial). The default
  `learners` stays `c("ares", "krls")`. Having no hyperparameters, `tune`
  is a no-op for them: every cross-fitting fold is a plain refit.

# roadrunner 0.0.0.9052

## plda() — Penalized Linear Discriminant Analysis

New user-facing function `plda()`: penalized Fisher's linear discriminant
analysis (Witten & Tibshirani 2011) with L1 and fused-lasso penalties,
multi-class support, built-in CV autotune, and parallel CV fold execution.

* **`plda(x, y, K, lambda, penalty, autotune, ...)`** — matrix and formula
  interfaces. Returns an S3 object of class `"plda"` carrying the discriminant
  matrix, class means, within-class standard deviations, and (when
  `autotune = TRUE`) the full CV result.
* **Penalties** — `penalty = "L1"` (default) applies coordinate-wise soft
  thresholding; `penalty = "fused"` adds a fused-lasso difference penalty
  (Condat 2013 total-variation prox) controlled by `lambda2`.
* **Multi-class** — supports `G >= 2` classes; `K` discriminant vectors up to
  `G - 1`.
* **Autotune** — `autotune = TRUE` (default) cross-validates `lambda` and `K`
  over a log-spaced grid; parallel CV fold loop via TBB with the roadrunner
  determinism invariant (byte-identical results for any `nthreads`).
* **S3 methods** — `predict.plda` (class, posterior, projection),
  `print.plda`, `summary.plda`, `plot.plda` (2-D discriminant projection).

# roadrunner 0.0.0.9051

## meep() — Phase Q-DML: cross-fitted causal ensemble

New user-facing function `meep()`: a pure-R orchestration layer over the
package's `ares()` (MARS) and `krls()` (Kernel Regularized Least Squares)
fitters that produces honest, cross-fitted (out-of-fold) nuisance
estimates for Double Machine Learning and causal-forest workflows.
`meep()` does not estimate a treatment effect itself — it returns the
cross-fitted predictions and residuals, which are meant to be handed
downstream to `grf::causal_forest()` (via `Y.hat` / `W.hat`) or to a DML
estimator.

* **`meep(X, y, treatment = NULL, ...)`** — a grf-style matrix/vector
  interface. Returns an S3 object of class `"meep"` carrying
  `y_hat_oof`, `d_hat_oof`, `mu0_hat_oof`/`mu1_hat_oof`, `y_resid`,
  `d_resid`, the `folds` vector, the per-nuisance OOF prediction
  matrices, the ensemble weights, and per learner-by-nuisance OOF
  performance.
* **Stacking ensemble** — `ensemble = "stack"` fits non-negative least
  squares on the honest OOF prediction matrix (leakage is already purged
  column-wise, so no second nesting layer is needed). `"average"` and
  `"best"` are cheaper fallbacks. The NNLS solver is hand-rolled
  (Lawson–Hanson active-set) — no new package dependencies.
* **`tune = c("once", "per_fold", "none")`** — `"once"` (default)
  autotunes each learner on the full data, freezes the hyperparameters,
  and refits per fold; `"per_fold"` autotunes inside every fold;
  `"none"` is the fast path that calls `ares()` / `krls()` with their
  own defaults and never invokes autotune.
* **Family auto-detection** — a binary outcome routes to
  `ares(family = "binomial")` + `krls(loss = "logistic")`; a binary
  treatment is modelled as a propensity score. Multi-valued treatments
  are rejected with a clear error.
* **Cluster-robust folds** — the `cluster` argument keeps whole clusters
  within a fold.
* **Graceful degradation** — a learner that errors on a fold is dropped
  for those rows, the ensemble weights are renormalized over the
  survivors per row, and the event is logged in `$fold_failures`.
* `predict.meep()`, `print.meep()`, and `summary.meep()` (with an
  overlap diagnostic for binary treatments) S3 methods.
* `grf` added to `Suggests` for the gated integration test.

# roadrunner 0.0.0.9050

## krls() — Phase Q5: true logistic-loss IRLS path (B1)

Closes the biggest functional gap vs `lukesonnet/krls2`: real logistic
ridge regression in the kernel basis, optimised via IRLS on penalised
binomial deviance. Default `loss = "ls"` remains byte-identical to
v0.0.0.9049 (back-compat snapshot in `test-krls-ard-kernel.R` continues
to pass to ULP tolerance, and `loss = "ls"` fits are byte-equal to
v9049 to `tolerance = 0`).

* **`loss = c("ls", "logistic")`** new argument on `krls.default()`.
  Default `"ls"` is the existing kernel ridge LS path. `"logistic"`
  optimises
  `minimise -2 sum(y log p + (1-y) log(1-p)) + lambda c' K c`
  where `p = plogis(K c)`, via a Cholesky-backed Newton step in
  coefficient space (option ii in the spec): each iter solves
  `(K W K + lambda I) c = K W z` with `W = p(1-p)` and working
  response `z = eta + (y - p) / W`. ~5-10 iters typical at default
  tolerance `1e-6`; cap 50. Step-halving (up to 5 per iter) on
  penalised-deviance increase. Numerical guards: W floored at `1e-8`,
  eta clipped to `[-30, 30]`, perfect-separation detection
  (`||c||_inf > 1e6` OR penalised dev `< 1e-10`) warns and returns
  the last finite iterate.

* **Compose rules.** `loss = "logistic"` requires binary y in `{0, 1}`.
  Combinations that don't make sense (or are deferred) error early:
  - `+ approx = "nystrom"` — deferred (large-n use case).
  - `+ lambda.method = "gcv"` — Wahba 2-D extension; deferred.
  - `+ lambda.method = "mll"` — Laplace approx; deferred to Q9.
  - `+ varmod != "none"` — residual variance model is LS-specific.
  - `predict(type = "variance")` — Q6 will derive the IRLS-aware
    posterior variance; current Q5 errors.

  Composes fine with: `whichkernel = "linear"`/`"poly*"`, `ard = "cheap"`,
  `autotune = TRUE` (deviance scoring; force R fallback because the
  Gaussian-specialised C++ inner is LS-specific), `n.boot > 0`
  (bagged averaging on the LINK scale, sigmoid at end for Jensen),
  `lambda.method = "loo"` (closed-form CT-2008), and
  `lambda.method = "cv"` (fold-deviance scoring), `derivative = TRUE`
  (link-scale marginal effects by default), `weights` (case weights
  baked into IRLS `W = w * p(1-p)`).

* **CT-2008 closed-form LOO under logistic.** Cawley & Talbot 2008
  leave-one-out approximation, computed in coefficient space from a
  single `H = K (K W K + lambda I)^{-1} K` solve. The diagonal
  `H_diag` is stashed on the fit (one extra n x n solve per converged
  IRLS). `cv_loo_dev` replaces `Looe` under logistic.

* **Predict gains `type = "class"`** (hard 0/1 prediction). For
  logistic fits: `type = "link"` returns eta; `type = "response"` and
  `type = "prob"` return calibrated probabilities; `type = "class"`
  returns `as.integer(p > 0.5)`. For binary LS fits, `type = "class"`
  uses the existing plogis-on-LS shortcut.

* **Print / summary.** Logistic fits print a
  `"loss: logistic (IRLS, k iter, converged)"` line, use McFadden
  pseudo-R-squared in place of LS R-squared, and append a
  `"(link scale)"` suffix on the average marginal effects header.

## Internal

* New C++ exports: `krls_irls_logistic_cpp` (IRLS solver returning
  coeffs, eta, p, W, deviance, iter, converged, H_diag),
  `krls_logistic_loo_loss_cpp` (CT-2008 closed-form LOO deviance).
* New R helpers: `.krls_logistic_lambdasearch`,
  `.krls_logistic_fit_block`, `.krls_bag_logistic`.
* `.krls_autotune_scalar()` accepts `is_logistic = FALSE`; when TRUE,
  per-cell IRLS fits + deviance scoring replace the LS solve.
* `.krls_run_cheap_ard()` accepts `loss` and forwards through the
  pass-1 isotropic fit so ARD under logistic works end-to-end.
* `slim_krls()` strips `eta_fitted` and `H_diag` along with the
  existing heavy intermediates. `slim_krls(keep_predict = FALSE)`
  retains `loss`, `deviance`, `converged`, `iter`, `cv_loo_dev` for
  the inspection-grade summary.
* `var.avgderivatives` under logistic is `NA` (full IRLS-aware
  derivation deferred to Q6) with a one-shot warning.

## Tests

* `test-krls-logistic.R` (~250 LoC, 13 tests): binary-y validation;
  IRLS convergence (< 20 iter); glmnet deviance match
  (`skip_if_not_installed("glmnet")`); Brier improvement vs LS+plogis;
  CT-2008 LOO matches explicit refit (n = 20 toy, tol 5e-3); CV under
  logistic; autotune + logistic; bagging + logistic; predict type
  consistency; link-scale marginal effects; compose rejections
  (Nystrom, GCV, MLL, varmod); ARD cheap + logistic on sparse signal
  lifts AUC; back-compat — `loss = "ls"` is byte-identical to v9049
  to `tolerance = 0`.

# roadrunner 0.0.0.9049

## krls() — Phase Q2: polynomial kernels + GP variance + MLL

Three bundled feature items implemented behind a single polymorphic
C++ kernel dispatch. Gaussian default fits are byte-identical to
v0.0.0.9048; non-Gaussian kernels open up `K_ij = x_i' x_j` (linear)
and `K_ij = (x_i' x_j + poly_c)^d` for `d` in 1..4 (polynomial).

* **A1 — `whichkernel = c("gaussian", "linear", "poly1..4")` + `poly_c`.**
  Default `"gaussian"` is unchanged. `"linear"` is the inner-product
  kernel; `"poly1..4"` are inhomogeneous polynomial kernels with
  offset `poly_c` (default `1.0`). Non-Gaussian kernels are
  incompatible with `approx = "nystrom"`, `ard != "none"`, and
  vector `sigma` — all error at fit time with clear messages.
  Marginal effects are computed via a closed-form dgemm route
  (`D = d * H * X` with `H = (X X' + c)^(d-1)` for poly, constant
  per column for linear). `var.avgderivatives` is currently `NA` for
  non-Gaussian kernels (full derivation deferred to Phase Q6); a
  one-shot warning surfaces this at fit time. Autotune dispatches:
  Gaussian sweeps `sigma` (9-pt anchor, existing); poly sweeps
  `poly_c` (default `c(0, 0.25, 0.5, 1, 2, 4, 8)`); linear collapses
  to lambda-only (errors if user supplies `autotune.grid`). Non-
  Gaussian autotune runs through the R fallback at the scalar autotune
  helper rather than the Gaussian-specialised `krls_autotune_inner_cpp`
  — expect a 5-10x wall-clock regression vs Gaussian autotune at
  matching grid sizes; this is acceptable for Q2 scope. The fit list
  gains `$whichkernel`, `$kernel_type`, `$poly_c`.

* **A6 — `predict(..., type = "variance")`.** New `type` value returns
  the GP posterior variance `K** - K* (K + lambda I)^{-1} K*'` per
  row of `newdata` via the new `krls_posterior_var_cpp` helper.
  Closed-form computation in the eigen-basis of the fit's kernel:
  reuses cached `dvals` + `V` (now stashed on the fit at v9049). New
  `unscale` arg (default `FALSE`): when `TRUE`, multiplies the
  returned variance by `var(y_train)` so it sits on the raw y scale.
  Incompatible with `approx = "nystrom"` (would require the
  m-length Phi spectrum — deferred). Bagged fits average per-
  replicate posterior variances. Matches a brute-force
  `Kss - diag(K* (K + lI)^{-1} K*')` reference to `tol = 1e-8` on
  n=50. FP-floored at 0.

* **A11 — `lambda.method = "mll"`.** Closed-form Type-II marginal
  log-likelihood as the lambda-selection objective. Pure-R loss
  `.krls_mll_loss(dvals, Vty, lambda, n)` reuses the same eigen-
  basis cache as LOO/GCV. Incompatible with `approx = "nystrom"`.
  Composes with autotune over sigma (Gaussian) and with non-Gaussian
  `whichkernel`. Sigma-by-MLL selection is deferred to Q8.

## Internal

* New C++ exports: `krls_posterior_var_cpp` (closed-form GP
  variance in eigen-basis). Existing kernel exports
  (`krls_kernel_cpp`, `krls_kernel_pred_cpp`, `krls_deriv_cpp`,
  `krls_avg_deriv_var_cpp`) gain optional `kernel_type` + `kernel_c`
  arguments with default `(0, 1.0)` (Gaussian, back-compat).
* New R helpers: `.krls_mll_loss()`, `.krls_predict_variance()`.
* `.krls_autotune_scalar()` accepts `kernel_type` + `poly_c`; when
  non-Gaussian, forces the R-side per-cell fallback (the
  Gaussian-specialised C++ inner is preserved untouched).
* `slim_krls()` strips the new `dvals` field along with the
  existing heavy intermediates. `slim_krls(keep_predict=FALSE)`
  retains `$whichkernel`, `$kernel_type`, `$poly_c` for the
  inspection-grade summary.
* `predict.krls_rr()` accepts new `type = "variance"` and `unscale`
  arguments; legacy fits (no `kernel_type`/`poly_c` field) default
  to Gaussian via fallback at predict time.

## Tests

* `test-krls-kernels-poly.R` (~150 LoC, 17 tests + 2 slow): linear
  kernel byte-equal vs hand-built `XX'`; poly2 byte-equal vs
  `(XX'+c)^2`; linear deriv constant per column; poly2 deriv matches
  finite difference to `tol = 1e-4`; validation gates; autotune
  poly_c sweep; bagging + linear; Gaussian back-compat seal.
* `test-krls-posterior-var.R` (~120 LoC, 8 tests + 1 slow): closed-
  form matches brute force to `tol = 1e-8`; training-point var =
  `1 - diag(H)`; monotonic increase with distance from training (1D
  toy); `unscale` arg; Nystrom errors; bagged averages; poly2
  posterior var diagonal; FP floor at 0.
* `test-krls-mll.R` (~95 LoC, 6 tests + 1 slow): closed-form matches
  hand-computed grid argmin; MLL ≠ LOO chosen lambda on real data;
  Nystrom + MLL errors; composes with autotune; composes with
  poly2.



## krls() — Phase Q1 parity batch

Four small parity items, all R-side (no C++ changes), batched
together to avoid four trivial release cycles.

* **A3 — `Neffective` (effective df).** Fits now carry
  `$Neffective`, the trace of the smoother hat matrix
  `H = K (K + lambda I)^{-1}`. Computed in the eigen-basis as
  `sum(d / (d + lambda))`, where `d` are the kernel eigenvalues
  used for the solve. The Nystrom path uses the m-length
  `Phi' Phi` spectrum so the formula generalises to
  `sum(Sigma2 / (Sigma2 + lambda))`. `summary(fit)` prints an
  `"Effective df: <value>"` line near the existing `R^2` line.
  `Neffective -> 0` as `lambda -> Inf` and `-> n` (or the rank
  of `K`) as `lambda -> 0`.
* **A7 — `subset` on `krls.default()`.** The matrix /
  data-frame default method now accepts `subset = NULL` (a
  no-op, back-compatible) or a logical / integer vector
  selecting rows of `X`. `y` and `weights` are sliced in
  lockstep before any downstream processing. Mirrors the
  formula method's `subset` semantics.
* **A8 — `predict(..., type = "prob")` for binary fits.**
  `predict.krls_rr()` gains a `type = c("response", "link",
  "prob")` argument. `"response"` (default) preserves the prior
  behaviour byte-for-byte. `"prob"` is valid only when the fit
  is binary (`y` has exactly two unique values, tracked on the
  fit via the new `$binary_y` field) and returns
  `plogis(yhat)`. The docstring is explicit: this is a
  calibration shortcut, **not** a true posterior probability;
  the underlying fit remains least-squares loss. `"link"` is a
  synonym for `"response"` reserved for the future
  logistic-loss extension.
* **A9 — `slim_krls()` / `unslim_krls()`.** New exported
  helpers for cheap serialisation. `slim_krls(fit,
  keep_predict = TRUE)` (default) strips the heavy `K`, `V`,
  `Vsq`, `Vty`, `vcov.c`, `vcov.fitted` intermediates while
  preserving everything `predict()` needs; saved size typically
  drops by >50% at `n >= 200`. `slim_krls(fit, keep_predict =
  FALSE)` strips further to an inspection-grade summary
  (`coeffs`, `R2`, `Looe`, marginal effects, sigma, lambda,
  `Neffective`); `predict()` then errors with a clear message.
  `unslim_krls()` is a no-op placeholder — refit `krls()` with
  the stored `$call` to rebuild the dropped fields.

# roadrunner 0.0.0.9047

## krls() — autotune unification + auto-ARD dispatch

* `autotune = TRUE` becomes a hands-free dispatcher. New argument
  `autotune.speed = c("balanced", "quality", "fast")` (default
  `"balanced"`) picks between the v0.0.0.9046 scalar sigma sweep and an
  ARD-dispatched path that routes through the cheap-tier ARD
  orchestrator:
  - `"fast"` keeps the v0.0.0.9046 scalar sigma grid behaviour
    byte-identically.
  - `"balanced"` (default) dispatches through cheap-tier ARD when
    `ncol(X) >= 20` or `ncol(X) >= nrow(X) / 10` (heuristic for high-p /
    sparse-signal regimes); otherwise scalar sigma only.
  - `"quality"` always dispatches through ARD and sweeps a 6-cell grid
    over `ard.alpha in {0.5, 1, 2}` x `ard.imp in {"avgderiv", "vsq"}`,
    picking the winner by inner K-fold CV.
* New argument `autotune.warmstart = TRUE` runs a cheap probe on a 15%
  subsample (capped at 200 rows) to test whether ARD dispatch
  outperforms the isotropic baseline; if the ARD probe does not improve
  held-out MSE by >2%, the ARD branch is dropped. Skipped when
  `n < 200` or `autotune.speed = "fast"`.
* The prior hard rejection of `ard != "none" + autotune = TRUE` is
  lifted. When the user pins `ard = "cheap"` explicitly under
  `autotune = TRUE`, autotune routes through ARD regardless of
  `autotune.speed`.
* `print.krls_rr()` now reports the autotune dispatch decision +
  selected `(alpha, imp)` for ARD-dispatched fits, and the chosen sigma
  for scalar fits.
* `$autotune` output schema gains `$speed`, `$warmstart`,
  `$ard_dispatched`, `$ard_decision_rule`, `$winner_sigma`,
  `$winner_alpha`, `$winner_imp`. `$winner` is preserved on the ARD path
  as the back-compat sentinel `NA_real_`.

### Behaviour change

> **Behaviour change at v0.0.0.9047:** `krls(autotune = TRUE)` now
> dispatches through cheap-tier ARD on high-p / sparse-signal data
> (`ncol(X) >= 20` or `ncol(X) >= nrow(X) / 10`) by default. To restore
> exact v0.0.0.9046 behaviour, pass `autotune.speed = "fast"`. On low-p
> dense data the default `"balanced"` mode is byte-identical to
> v0.0.0.9046.

### Out of scope (P2c / P3 / P5)

* Gradient-based / marginal-likelihood ARD optimisation.
* HSIC pre-screen for high-p feature selection.
* Vector-sigma autotune sweep (requires C++ rewrite of
  `krls_autotune_inner_cpp`).
* Nystrom + autotune-ARD composition.
* Multi-objective tuning (R^2 + sparsity).
* `ssf_grid` beyond `c(1.0)`; multi-threaded ARD cell grid.

# roadrunner 0.0.0.9046

## krls() — automatic two-pass ARD selector (cheap tier)

* New argument `ard = "cheap"` enables a two-pass automatic ARD pipeline:
  pass 1 fits isotropic KRLS at the scale-aware sigma anchor; pass 2
  refits with per-feature lengthscales derived from pass-1 marginal-
  effect importances via `s_k = sigma_iso * (median(imp) / imp_k)^alpha`,
  clipped to `[sigma_iso / cap, sigma_iso * cap]`.
* New args `ard.alpha = 1.0` (mapping exponent), `ard.cap = 100`
  (symmetric multiplicative ceiling), `ard.imp = c("avgderiv", "vsq")`
  (importance source).
* Default `ard = "none"` is byte-identical to v0.0.0.9045.
* Composes with bagging (`n.boot`), CV (`lambda.method = "cv"`), GCV,
  weights, and varmod. Rejects autotune + ARD and Nystrom + ARD at fit
  time (same constraints as manual vector sigma in P2a).
* Compute cost: ~3-5x scalar isotropic. Pass 1 is a throwaway fit;
  its `varmod` / `binary` / `n.boot` are forced off, but `vcov`
  stays on because the engine couples it to `derivative = TRUE`
  (which we need for importance extraction). The pass-1 outputs are
  discarded after the per-feature importance vector is read off.
* Empirical lift: on grf::generate_causal_data aw3 with n_train=750,
  p=100 (sparse signal, mostly noise dimensions), test R^2 rises
  from 0.13 (isotropic) to 0.33 (ard='cheap'), comparable to
  ranger (0.34) on the same data.

### Out of scope (P2c / P3 / P5)

* Gradient-based / marginal-likelihood ARD optimisation.
* HSIC pre-screen for high-p feature selection.
* Autotune over per-feature sigma.
* Nystrom + ARD composition.

# roadrunner 0.0.0.9045

## krls() — per-feature Gaussian bandwidth (manual ARD)

* `krls(..., sigma = c(s1, ..., sp))` now accepts a length-`ncol(X)`
  vector of strictly positive lengthscales. Gaussian kernel becomes
  K_ij = exp(-sum_k (x_ik - x_jk)^2 / sigma_k).
* Scalar sigma path is byte-identical to v0.0.0.9044 fits (FP determinism
  preserved via constant-vector fast detection in kernel worker).
* Marginal effects, vcov of average effects, and predict() all scale
  per-feature by 2/sigma_k and 4/sigma_k^2 respectively.
* `autotune = TRUE` + vector sigma errors at fit time; auto-selection of
  ARD lengthscales deferred to Phase 2b/2c.
* `approx = "nystrom"` + vector sigma errors at fit time; ARD + Nystrom
  is P5 stretch goal.

### Out of scope (P2b/P2c/P3)

* Cheap two-pass ARD selector.
* Gradient-based ARD optimisation.
* HSIC pre-screen.
* Autotune sweep over ARD vectors.

# roadrunner 0.0.0.9044

## krls() — GCV lambda selection

* `krls(..., lambda.method = "gcv")` selects the ridge penalty by
  minimising the closed-form generalised cross-validation criterion
  (Craven & Wahba 1979) on the existing eigendecomposition of `K`.
  Brings parity with `KRLS` v1.5-0+.
* Closed form: no extra kernel evaluation or refit; cost is one
  golden-section search over the same `[L, U]` bracket used by LOO.
* Denominator `(1 - trH/n)^2` is floored at `1e-8` to remain finite
  in the interpolation limit `lambda -> 0`.
* `lambda.method = "loo"` remains the default. `"loo"` and `"cv"`
  paths are byte-identical to v0.0.0.9043 on a fixed-seed problem.

### Out of scope (Phase 2/3)

* ARD-kernel option (per-feature sigma).
* HSIC pre-screen for high-p feature selection.
* Autotune-dispatch routing through GCV (autotune still calls inner
  CV regardless of `lambda.method`).

# roadrunner 0.0.0.9043

## krls() speedup — Phase 2 (Nystrom low-rank approximation)

* New optional argument `approx = "nystrom"` enables the Nystrom
  low-rank approximation. Replaces O(n^3) eigendecomposition on the
  full kernel with O(m^3) + O(n m^2) where m = `nystrom_m`
  (default `ceiling(sqrt(n) * 3)`, e.g. n=2000 -> m=135,
  n=5000 -> m=213).
* Five new optional args: `approx`, `nystrom_m`, `landmarks`,
  `landmark_method`, `landmark_seed`, `nystrom_eps`.
* Default exact path is byte-identical to v0.0.0.9042 — opt-in only.
* New exported helper `get_landmarks(fit)` returns landmark
  coordinates (original or standardized X scale).
* Determinism: at fixed `landmark_seed` (+ `seed.cv` when autotune),
  Nystrom fits are byte-identical across `autotune.nthreads`.

### Empirical wall-clock speedup (EMP-PHASE2, R=5 reps, paired seeds)

| n    | p  | approx   | autotune | wall (s) | speedup | paired RMSE delta (median) |
|------|----|----------|----------|----------|---------|----------------------------|
| 2000 | 10 | exact    | FALSE    | 0.52     | 1.0x    | -                          |
| 2000 | 10 | nystrom  | FALSE    | 0.056    | 9.27x   | -0.02%                     |
| 5000 | 10 | exact    | FALSE    | 10.08    | 1.0x    | -                          |
| 5000 | 10 | nystrom  | FALSE    | 0.329    | 30.64x  | +4.81%                     |

Paired RMSE delta is the per-rep `(RMSE_nystrom - RMSE_exact) /
RMSE_exact` on identical (X, y) draws (paired seed schedule keyed on
`(n, rep)`); n=2000 range [-3.4%, +5.4%], n=5000 range [+4.0%, +8.9%].
Both well within the +10% target on the additive smooth DGP.

### API compatibility

* Zero breaking changes. `krls()` without `approx` defaults to `"exact"`
  and behaves identically to v0.0.0.9042.
* `predict.krls_rr` auto-detects Nystrom fits via `fit$approx` and uses
  the cheap cross-kernel path.

### Out of scope (Phase 3)

* Leverage-score landmark selection.
* Auto-engage Nystrom at large n.
* Sparse kernel truncation.
* Bagging + Nystrom integration.

# roadrunner 0.0.0.9042

## krls() speedup — Phase 1 (shared distance + parallel autotune)

* `krls(..., autotune = TRUE)` is now parallelised over the sigma grid via
  RcppParallel/TBB. The pairwise squared-distance matrix is computed once
  per CV fold and reused across every sigma candidate.
* New optional argument `autotune.nthreads` (default
  `getOption("roadrunner.nthreads", parallel::detectCores(logical=FALSE))`).
  Pass `autotune.nthreads = 1` for strictly sequential execution.
* Determinism contract preserved: fits are byte-identical across
  `autotune.nthreads` values at fixed seed and inputs (each worker writes
  to a unique slot in the output vectors, so no reduction is involved).

### Empirical wall-clock speedup (REQ-001 minimal grid, R=5)

| n    | p  | nthreads | time   | speedup vs v0.0.0.9041 |
|------|----|----------|--------|------------------------|
| 500  | 10 | 1        | (TBD)  | 1.0x (sequential)      |
| 500  | 10 | 4        | (TBD)  | (filled by Task 9)     |
| 1500 | 20 | 1        | (TBD)  | 1.0x (sequential)      |
| 1500 | 20 | 4        | (TBD)  | (filled by Task 9)     |

### API compatibility

* Zero breaking changes. `krls(..., autotune = TRUE)` returns the same
  S3 object with the same fields. `autotune_info` gains
  `nthreads_used` and `sigma_grid_sorted` slots.
* `predict.krls_rr`, `summary.krls_rr`, `print.krls_rr` unchanged.

### Out of scope (Phase 2)

* Nystrom approximation (`approx = "nystrom"`, `nystrom_m`, `landmarks`).
* Predict()-side speedups.
* Multi-response / sparse kernel paths.

# roadrunner 0.0.0.9041

## `krls()` — scale-aware sigma anchor refined to geomean_p (REQ-20260518-003)

The default sigma anchor formula in `.krls_sigma_anchor()` is updated from the
raw median heuristic (`median(d2)`, introduced in v0.0.0.9040) to the
*geomean_p* formula: `sqrt(median(d2) * p)` where `d2` are pairwise squared
Euclidean distances on the standardised training matrix and `p = ncol(X)`.

**Why**: a 15-DGP head-to-head vs `KRLS::krls()` at `n=500, p=10` (iter-0,
REQ-20260518-003) showed the v0.0.0.9040 median anchor over-smoothed locally
nonlinear signals, producing 8 losses and 0 ties in favour of KRLS.  An
anchor sweep over 6 candidate formulae (iter-1) identified geomean_p as the
clear winner (14/15 DGPs).  Verification with the patched anchor (iter-2)
confirmed: **roadrunner wins 12/15 DGPs, KRLS wins 0/15, 3 ties**.

**Delta table — iter-0 (median) vs iter-2 (geomean_p) mean MSE ratio (rr / KRLS)**:

| DGP             | iter-0 ratio | iter-2 ratio | delta   |
|-----------------|-------------|-------------|---------|
| additive        | 1.327        | 0.896        | -0.431  |
| exp-decay       | 1.681        | 0.887        | -0.795  |
| friedman1       | 1.126        | 0.959        | -0.167  |
| friedman3       | 1.084        | 0.957        | -0.127  |
| interaction     | 1.471        | 0.895        | -0.576  |
| mixture         | 1.140        | 0.992 (tie)  | -0.148  |
| poly2           | 1.286        | 0.871        | -0.415  |
| tanh-interaction| 1.485        | 0.952        | -0.533  |
| friedman2       | 0.984 (tie)  | 0.936        | -0.049  |
| sin-sum         | 1.011 (tie)  | 0.982 (tie)  | -0.029  |
| linear          | 0.872        | 0.919        | +0.047  |
| sparse          | 0.865        | 0.917        | +0.052  |
| heterosked      | 0.859        | 0.920        | +0.062  |
| monotone        | 0.960        | 0.925        | -0.035  |
| noise           | 1.000 (tie)  | 1.000 (tie)  |  0.000  |

At `n=500, p=10` on standardised N(0,1) data the geomean_p anchor gives
`sigma ~ 13.5`, splitting the difference between the raw median (~18–19) and
KRLS's fixed `sigma = p = 10`.  The formula is still data-adaptive: on
non-standardised or non-unit-variance inputs it will differ sensibly from both
extremes.

This is a **one-line R change** (`R/krls.R`, `.krls_sigma_anchor()`); the C++
engine and all other logic are unchanged.  At a fixed `(sigma, lambda)` fits
remain byte-identical to all prior versions.  To restore the v0.0.0.9040
median anchor, pass
`sigma = stats::median(as.numeric(stats::dist(scale(X)))^2)`.

# roadrunner 0.0.0.9040

## `krls()` — four overfitting fixes (REQ-20260518-002)

Addresses overfitting confirmed by diagnostic sweeps in REQ-20260518-001.
All four changes are R-level only; the C++ engine is unchanged, so at a
fixed `(sigma, lambda)` fits remain byte-identical to earlier versions.

Empirical improvement on signal DGPs (`n=500, p=10`, `tune=none`, `R=10`
replications; overfit ratio = `test_MSE / train_MSE`):

| DGP         | Before (pre-fix) | After (post-fix) | Improvement |
|-------------|-----------------|-----------------|-------------|
| additive    | 3.35            | 1.58            | 53%         |
| interaction | 4.38            | 2.10            | 52%         |
| sparse      | 2.61            | 1.26            | 52%         |
| linear      | 2.70            | 1.37            | 49%         |
| noise       | 1.06            | 1.04            | —           |

All improvement is attributable to Fix 1 (sigma default) and Fix 2 (lambda
tolerance); Fixes 3 and 4 further stabilise autotune sigma selection.

**Breaking defaults** (old values restorable by explicit arguments):

- **Default `sigma` changed** from `ncol(X)` to the median pairwise
  squared Euclidean distance on the standardised predictors (the
  'median heuristic'). At `n=500, p=10` the oracle sigma is ~20 and the
  median heuristic anchors near that neighbourhood; `sigma = ncol(X)`
  was 10 (2x under-smoothed). Restore old behaviour with
  `sigma = ncol(X)`.

- **Default lambda tolerance changed** from `1e-3 * n` (n-dependent,
  coarse at moderate n) to `1e-6` (fixed, 6-digit precision). The LOO
  golden-section was empirically selecting lambda ~4x the argmin at
  `sigma=20, n=500` under the old tolerance. The L-bracket climb now
  uses multiplicative steps (x10 per step) instead of additive steps
  (0.05 per step) for scale-robustness. Restore old behaviour with
  `tol = 1e-3 * nrow(X)`.

- **Autotune sigma grid changed** from `ncol(X) * c(0.25, 0.5, 1, 2,
  4, 8)` (6 points fixed at `d`) to `sigma_anchor * c(0.125, 0.25,
  0.5, 1, 2, 4, 8, 16, 32)` (9 points centred on the median heuristic
  anchor). Restore old behaviour with `autotune.grid = ncol(X) *
  c(0.25, 0.5, 1, 2, 4, 8)`.

- **Autotune CV stabilised**: default folds raised from 5 to 10; default
  cross-partition repeats raised from 1 to 2 (via the `ncross` argument,
  whose default changes from `1L` to `NULL`); sigma selection now applies
  the 1-SE rule (largest sigma within 1 SE of minimum CV-MSE, biasing
  toward wider kernels when evidence is weak). The `autotune` component
  of the fit gains new fields: `ncross`, `mse_per_fold`, `se_mse`,
  `cv.1se`, `sigma_1se`. Restore old fold count with `nfold = 5`;
  restore old repeat count with `ncross = 1`.

# roadrunner 0.0.0.9033

## New feature: `krls()` -- Kernel Regularized Least Squares

- Adds `krls()`, a from-scratch implementation of the Hainmueller and
  Hazlett (2014) KRLS estimator under the roadrunner roof. The
  algorithm mirrors `KRLS::krls()` exactly (standardisation, Gaussian
  kernel, eigen-basis closed-form solve, golden-section LOO lambda
  search, marginal effects, binary first-difference handling); at
  matched `(sigma, lambda)` fits agree with `KRLS::krls()` to within
  floating-point precision (`< 1e-12` on coefficients, fitted values,
  pointwise and average marginal effects, and prediction SEs).
- Engine is C++ (`src/krls.cpp`) on top of `RcppArmadillo` (added to
  `LinkingTo` + `Imports`) and `RcppParallel`. The Gaussian kernel and
  the test-vs-train kernel are built in parallel with a TBB worker.
  Eigendecomposition is dispatched to LAPACK via `arma::eig_sym` with
  the divide-and-conquer driver. Marginal effects are computed via the
  identity `dy/dx_k = -(2/sigma)*(X_k*(K c) - K diag(c) X)_k`, which
  avoids the explicit `n x n` distance matrix and reduces the average-
  marginal-effect variance from `O(n^3)` to `O(n^2)` per variable via a
  row-sum trick.
- Measured speed-up over `KRLS::krls()` on simulated benchmarks
  (`derivative = TRUE`, `vcov = TRUE`) at `nthreads = default`:
  `n=200 p=3` ~2x; `n=500 p=3` ~5-6x; `n=500 p=10` ~10x;
  `n=1000 p=3` ~6x; `n=1000 p=10` ~10x. Coefficient max-abs error vs
  `KRLS::krls()` is `< 2e-13` on every cell.
- S3 class is `c("krls_rr", "krls")` so `predict()`, `print()`, and
  `summary()` dispatch unambiguously to roadrunner methods even when
  `KRLS` is loaded in the same session. `inherits(fit, "krls")` is
  preserved for downstream-compat checks.
- Tests: `tests/testthat/test-krls.R` (parity vs `KRLS::krls` for
  coefficients, fitted values, Looe, marginal effects, prediction SEs,
  binary first-differences; structural input-validation; predict
  recovers fitted values; sensible R^2). 11 tests, 28 assertions, all
  green.

# roadrunner 0.0.0.9032

## Bug fixes (statsclaw 2026-05-13 audit triage, BUG-008..BUG-013)

- **BUG-013 (usability, low)**: `print.ares()` and `summary.ares()` used
  to be silent about bagging (`n.boot > 0`) and autotune state, so a
  bagged or autotuned fit printed identically to a plain fit. Fix: add
  one-line `Bagging: n.boot = N replicate(s)` and `Autotune: degree=D
  penalty=P nk=K fast.k=F warmstart=T/F` blocks to both printers when
  the corresponding components are present. `summary.ares` also carries
  `$boot` and `$autotune` through to its print method. Regression test:
  `tests/testthat/test-bug-013-print-bag-autotune.R`.

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
