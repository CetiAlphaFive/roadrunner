# Cross-fitted causal ensemble of `ares()` and `krls()`

`meep()` produces honest, cross-fitted (out-of-fold) nuisance estimates
for use in Double Machine Learning (DML) and causal-forest workflows.
For an outcome `y`, an optional `treatment`, and covariates `X`, it
cross-fits an ensemble of the package's
[`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)
(MARS) and
[`krls()`](https://cetialphafive.github.io/roadrunner/reference/krls.md)
(Kernel Regularized Least Squares) learners and returns the out-of-fold
predictions \\\hat E\[Y\mid X\]\\ and (when `treatment` is supplied)
\\\hat E\[D\mid X\]\\, plus the residuals \\Y-\hat Y\\ and \\D-\hat D\\.

## Usage

``` r
meep(
  X,
  y,
  treatment = NULL,
  folds = 5L,
  learners = c("ares", "krls", "ols", "logreg"),
  ensemble = c("stack", "average", "best"),
  arm_models = c("auto", "always", "never"),
  family = NULL,
  treatment_family = NULL,
  tune = c("once", "per_fold", "none"),
  cluster = NULL,
  weights = NULL,
  seed = NULL,
  ares_args = list(),
  krls_args = list(),
  ols_args = list(),
  logreg_args = list(),
  verbose = FALSE,
  ...
)
```

## Arguments

- X:

  An n-by-p numeric matrix or data.frame of covariates. Factor columns
  in a data.frame are expanded by the base fitters.

- y:

  Length-n outcome vector.

- treatment:

  Optional length-n treatment / exposure vector. `NULL` gives an
  outcome-only fit (no `d_hat_oof`, no arm models).

- folds:

  Either a single integer `K` (number of folds), or a length-n integer
  vector of fold ids that is honored verbatim.

- learners:

  Character subset of `c("ares", "krls", "ols", "logreg")`, or a named
  list of adapter closures (extensibility hook). The default is all
  four: `c("ares", "krls", "ols", "logreg")`. `ares` and `krls` are
  family-agnostic and apply to every nuisance. `"ols"` and `"logreg"`
  are unregularized linear learners with no hyperparameters (so `tune`
  is a no-op for them – a plain refit per fold); they are applied *per
  nuisance by family*, `ols` for gaussian (continuous) targets and
  `logreg` for binomial (binary) ones. A learner that does not apply to
  a nuisance is skipped, leaving an all-`NA` OOF column that the
  ensemble excludes; it is not counted as a fold failure. Narrowing
  `learners` to names that are all family-incompatible with a nuisance
  (for example `learners = "ols"` with a binary outcome) is an error.
  Custom list-learners are never family-filtered.

- ensemble:

  How to combine learners: `"stack"` (non-negative least squares on the
  OOF matrix), `"average"` (equal weights), or `"best"` (pick the single
  lowest-OOF-loss learner).

- arm_models:

  `"auto"` fits arm-specific outcome models (\\\mu_0\\, \\\mu_1\\) only
  when the treatment is binary; `"always"` forces them; `"never"`
  suppresses them.

- family:

  Outcome family, `"gaussian"` or `"binomial"`. `NULL` auto-detects.

- treatment_family:

  Treatment family, `"gaussian"` or `"binomial"`. `NULL` auto-detects.

- tune:

  One of `"once"`, `"per_fold"`, `"none"`. See *Tuning*.

- cluster:

  Optional length-n grouping vector. Folds are assigned at the cluster
  level so each cluster lands entirely within one fold.

- weights:

  Optional length-n case weights, forwarded to every base fit and into
  the NNLS stacking objective.

- seed:

  Optional integer seed for fold construction.

- ares_args:

  A named list of extra arguments spliced into every
  [`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)
  call.

- krls_args:

  A named list of extra arguments spliced into every
  [`krls()`](https://cetialphafive.github.io/roadrunner/reference/krls.md)
  call.

- ols_args:

  A named list of extra arguments spliced into every
  [`ols()`](https://cetialphafive.github.io/roadrunner/reference/ols.md)
  call (only relevant when `"ols"` is among `learners`).

- logreg_args:

  A named list of extra arguments spliced into every
  [`logreg()`](https://cetialphafive.github.io/roadrunner/reference/logreg.md)
  call (only relevant when `"logreg"` is among `learners`).

- verbose:

  Logical; if `TRUE`, print progress per nuisance / fold.

- ...:

  Currently unused; reserved.

## Value

An object of S3 class `"meep"`: a list with `y_hat_oof`, `d_hat_oof`,
`mu0_hat_oof`, `mu1_hat_oof`, `y_resid`, `d_resid`, `folds`,
`oof_matrix` (per-nuisance n-by-L matrices), `ensemble_weights`,
`learner_cv_perf`, `learners`, `family`, `treatment_family`, `tune`,
`ensemble`, `n`, `p`, `frozen_hyperparams`, `fold_failures`, `cluster`,
`weights`, `seed`, and `call`.

## Details

`meep()` does **not** estimate a treatment effect. The honest OOF
predictions are meant to be handed to a downstream DML estimator or to
[`grf::causal_forest()`](https://rdrr.io/pkg/grf/man/causal_forest.html)
(see Examples). Because each column of the OOF prediction matrix already
comes from a model that never saw the row it predicts, the stacking
ensemble (`ensemble = "stack"`) can fit non-negative least squares
directly on that matrix without a second layer of nesting – leakage is
already purged column-wise.

## Family auto-detection

When `family` is `NULL`, a binary `y` (exactly two distinct values)
triggers `ares(family = "binomial")` and `krls(loss = "logistic")`;
otherwise the gaussian path is used. The same logic applies to
`treatment_family`: a binary treatment is modelled as a propensity score
via logistic regression, a continuous treatment via gaussian regression.
A treatment with three or more distinct values is rejected (multi-valued
/ continuous-treatment IRM is out of scope).

The family also governs which default learners are fitted for each
nuisance. `ares` and `krls` are family-agnostic and always apply. `ols`
is fitted only for gaussian (continuous) nuisances and `logreg` only for
binomial (binary) nuisances. Because the Y-model and the D-model can
differ in family, applicability is resolved *per nuisance*: a gaussian
outcome with a binary treatment fits `ols` for `outcome`/`mu0`/`mu1` and
`logreg` for `treatment`. A learner that does not apply to a nuisance is
simply skipped – its OOF column stays all-`NA` and it is *not* recorded
as a fold failure.

## Tuning (`tune`)

- `"once"` (default) – autotune each learner on the full data, freeze
  the resulting hyperparameters, and refit those frozen settings within
  every fold. Standard DML practice: fit honesty (cross-fitting) is what
  protects inference; freezing hyperparameters is a deliberate, second-
  order compromise.

- `"per_fold"` – autotune independently inside every fold. Purist
  option; slower.

- `"none"` – the fast path. Call
  [`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)
  and
  [`krls()`](https://cetialphafive.github.io/roadrunner/reference/krls.md)
  with their own default arguments (plus anything passed via `ares_args`
  / `krls_args`). No autotune is invoked at all.

## Graceful degradation

If a learner errors on a fold (for example, a constant nuisance in a
fold subset), its column for that fold becomes `NA`, the ensemble
weights are renormalized over the surviving learners on a per-row basis,
and the event is logged in `$fold_failures`. If *every* learner fails a
fold, `meep()` stops.

## DoubleML recipe

To use these nuisances with the DoubleML package, pass `y_hat_oof` /
`d_hat_oof` as `external_predictions` and reuse `fit$folds` as the
sample split so the cross-fitting partition matches.

## See also

[`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md),
[`krls()`](https://cetialphafive.github.io/roadrunner/reference/krls.md)

## Examples

``` r
# \donttest{
set.seed(1)
n <- 800; p <- 5
X <- matrix(runif(n * p), n, p)
m  <- 0.4 * X[, 1] + 0.3 * sin(3 * X[, 2])
D  <- m + rnorm(n, sd = 0.3)
g  <- X[, 1]^2 + 0.5 * X[, 3]
Y  <- 1.0 * D + g + rnorm(n, sd = 0.5)        # true ATE = 1.0
fit <- meep(X, Y, treatment = D, folds = 5, seed = 42)
fit
#> meep(): cross-fitted causal ensemble
#> Call: meep(X = X, y = Y, treatment = D, folds = 5, seed = 42)
#> 
#>   n = 800,  p = 5
#>   folds: 5-fold cross-fitting
#>   learners: ares, krls, ols, logreg
#>   outcome family: gaussian  |  treatment family: gaussian
#>   ensemble: stack  |  tune: once
#> 
#> Ensemble weights:
#>   outcome  : ares=0.213  krls=0.066  ols=0.720  logreg=0.000
#>   treatment: ares=0.704  krls=0.296  ols=0.000  logreg=0.000
#> 
#> OOF loss (per nuisance, ensemble):
#>   outcome  : mse = 0.378
#>   treatment: mse = 0.09565
# Hand the honest OOF nuisances to a causal forest:
# grf::causal_forest(X, Y, D, Y.hat = fit$y_hat_oof, W.hat = fit$d_hat_oof)
# }
```
