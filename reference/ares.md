# Fast Multivariate Adaptive Regression Splines

Fits a Multivariate Adaptive Regression Splines (MARS) model (Friedman
1991). A forward pass adds pairs of hinge basis functions of the form
`max(0, +/-(x - knot))`; a backward pass prunes terms by minimising the
GCV criterion (or a K-fold CV criterion if requested).

## Usage

``` r
ares(x, ...)

# S3 method for class 'formula'
ares(x, data = NULL, ..., y = NULL)

# Default S3 method
ares(
  x,
  y,
  degree = 1L,
  nk = NULL,
  penalty = NULL,
  thresh = 0.001,
  minspan = 0L,
  endspan = 0L,
  adjust.endspan = 1L,
  auto.linpreds = FALSE,
  fast.k = 10L,
  fast.beta = 1,
  nprune = NULL,
  pmethod = c("backward", "none", "cv"),
  nfold = 0L,
  ncross = 1L,
  stratify = TRUE,
  seed.cv = NULL,
  cv.1se = FALSE,
  autotune = FALSE,
  autotune.speed = c("balanced", "quality", "fast"),
  autotune.warmstart = TRUE,
  n.boot = 0L,
  na.action = c("impute", "omit"),
  family = c("gaussian", "binomial", "poisson", "gamma"),
  weights = NULL,
  varmod = c("none", "const", "lm"),
  trace = 0L,
  nthreads = 0L,
  ...
)
```

## Arguments

- x:

  A numeric matrix or data frame of predictors, or a model formula when
  calling the formula method.

- ...:

  Currently ignored.

- data:

  A data frame. Used only by the formula method.

- y:

  A numeric response vector with length `nrow(x)`. For
  `family = "binomial"`, a 0/1 numeric, logical, or 2-level factor.

- degree:

  Maximum interaction degree. Default `1` (additive). Use `2` or `3` for
  two- or three-way interactions.

- nk:

  Maximum number of basis terms in the forward pass. Default scales with
  `ncol(x)`: `min(200, max(20, 2 * ncol(x))) + 1`.

- penalty:

  GCV penalty per knot. Default `2` for `degree = 1`, `3` otherwise.
  Larger values produce sparser fits.

- thresh:

  Forward-pass early-stop threshold on relative RSS improvement. Default
  `0.001`.

- minspan:

  Minimum gap between knots along a variable. `0` (default) picks an
  automatic value from `n` and `p`.

- endspan:

  Distance from the data ends within which knots are forbidden. `0`
  (default) picks an automatic value.

- adjust.endspan:

  Multiplier applied to `endspan` when the candidate hinge would deepen
  an existing interaction. Default `1` for `gaussian`, `poisson`,
  `gamma`; `2` for `binomial`. Pass an explicit value to override.

- auto.linpreds:

  If `TRUE`, hinge pairs whose best knot sits at a variable's range
  boundary are replaced by a linear term. Default `FALSE` for
  `gaussian`, `poisson`, `gamma`; `TRUE` for `binomial`. Pass an
  explicit value to override.

- fast.k:

  Size of the Fast-MARS candidate cache (Friedman 1993). Larger values
  rescore more candidates each step (slower, slightly more accurate);
  `0` rescores every candidate every step. Default `10`.

- fast.beta:

  Age penalty for stale entries in the Fast-MARS cache. Default `1.0`.
  Only relevant when `fast.k > 0`.

- nprune:

  Maximum number of terms after backward pruning. Default `nk`.

- pmethod:

  Pruning method.

  - `"backward"` (default): minimise GCV over backward-elimination
    subsets.

  - `"none"`: keep all forward-pass terms.

  - `"cv"`: K-fold cross-validated subset selection. Requires
    `nfold > 0` (sets `nfold = 5` if not specified).

- nfold:

  Number of CV folds for `pmethod = "cv"`. Default `0` (no CV; GCV-based
  pruning). When `nfold > 0`, `pmethod` is promoted to `"cv"`.

- ncross:

  Number of CV repetitions (each builds a fresh fold partition; per-size
  mean MSE is averaged). Default `1`.

- stratify:

  If `TRUE` (default), CV folds are quantile-stratified on `y`. Useful
  for small `n` or skewed responses.

- seed.cv:

  Optional integer seed for the CV fold partition. Pass an integer for
  reproducible CV; `NULL` (default) uses the current RNG state.

- cv.1se:

  If `TRUE`, applies the one-standard-error rule under `pmethod = "cv"`:
  among sizes within one SE of the minimum mean CV-MSE, pick the
  smallest. Default `FALSE` (pick the argmin).

- autotune:

  If `TRUE`, runs an inner cross-validated grid search over
  `(degree, penalty, nk, fast.k)` and refits the winner on the full
  data. The chosen settings and full grid scores are returned in
  `$autotune`. Default `FALSE`.

- autotune.speed:

  Speed/quality trade-off for `autotune`. Only used when
  `autotune = TRUE`.

  - `"balanced"` (default): explores a moderate `fast.k` range and picks
    the smallest whose CV-MSE is within 1% of the best.

  - `"quality"`: disables the Fast-MARS cache (most thorough, slowest).

  - `"fast"`: forces an aggressive cache (cheapest, slight accuracy
    trade-off).

- autotune.warmstart:

  If `TRUE` (default), `autotune` first tunes on a small subsample. If
  one cell wins decisively, the full grid is skipped and the winner is
  refit on all rows. Only applies when `autotune = TRUE` and `n >= 200`.

- n.boot:

  Number of bootstrap replicate fits for bagging. Default `0` (no
  bagging). When `> 0`,
  [`predict()`](https://rdrr.io/r/stats/predict.html) averages over the
  replicates plus the central fit and (with `se.fit = TRUE`) returns a
  per-prediction bag standard deviation. Composes with `autotune`: each
  replicate reuses the central fit's chosen hyperparameters.

- na.action:

  Strategy for missing values in `x`. `"impute"` (default) replaces each
  column's `NA`s with that column's training median and stores the
  medians for [`predict()`](https://rdrr.io/r/stats/predict.html) to
  reuse. `"omit"` drops rows with any `NA`. Either action warns. Missing
  values in `y`, or `NaN`/`Inf` anywhere in `x`, are always rejected.

- family:

  Response family for the final coefficient fit.

  - `"gaussian"` (default): identity link, OLS.

  - `"binomial"`: logit link. `y` may be 0/1 numeric, logical, or a
    2-level factor.

  - `"poisson"`: log link. `y` must be non-negative integer-valued.

  - `"gamma"`: log link. `y` must be strictly positive. For non-gaussian
    families, term selection runs on the numeric `y` (fast); the
    selected basis is then refit on the response scale via
    [`stats::glm.fit()`](https://rdrr.io/r/stats/glm.html) with the
    requested family.

- weights:

  Optional non-negative numeric vector of length `nrow(x)` for weighted
  least-squares fitting. Default `NULL` (unweighted). Term selection,
  GCV / CV pruning, and autotune all respect the weights. Negative,
  `NA`, or `NaN` weights are rejected.

- varmod:

  Residual variance model used by
  [`predict()`](https://rdrr.io/r/stats/predict.html) for prediction
  intervals (gaussian only).

  - `"none"` (default): no variance model is stored; `interval = "pint"`
    is unavailable.

  - `"const"`: stores a single residual SD; intervals are
    `yhat +/- qt(level, df) * sigma`.

  - `"lm"`: fits a small linear model of `|resid|` on `yhat` to allow
    simple **yhat-dependent** heteroscedasticity. Captures residual
    scale that changes linearly with the fitted mean. Does NOT capture
    residual scale that depends on a predictor whose contribution to
    `yhat` is small – in that case the \|resid\| ~ yhat slope is close
    to zero and `"lm"` collapses to roughly the same as `"const"`. If
    you suspect x-driven heteroscedasticity (e.g. variance depends on a
    covariate orthogonal to the mean structure), `varmod = "lm"` will
    not help and coverage will degrade in high-variance regions. Ignored
    for non-gaussian families.

- trace:

  Trace level. `0` (default) is silent; `1` reports forward-pass
  progress.

- nthreads:

  Number of threads. `0` (default) uses
  [`RcppParallel::defaultNumThreads()`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html).
  Examples and the vignette cap this at `2` for CRAN compliance.

## Value

An object of class `"ares"`: a list containing the fitted coefficients,
the basis matrix `bx`, the term directions `dirs` and knots `cuts`, the
indices of `selected.terms`, training `rss` and `gcv`, `fitted.values`
and `residuals`, predictor names `namesx`, the call, and echoed control
parameters. Non-gaussian fits additionally include `$family`, `$glm`
(with `deviance`, `null.deviance`, `df.null`, `df.residual`, `aic`,
`converged`, `iter`), and `$linear.predictor`; `$fitted.values` is on
the response scale (probabilities for binomial, positive means for
poisson and gamma). Autotune fits carry `$autotune` and bagged fits
carry `$boot`.

## Details

Two interfaces are provided. The formula method takes a model formula
and a data frame. The default method takes a numeric predictor matrix
`x` and a numeric response vector `y`.

Fits are deterministic across thread counts: at a fixed `seed.cv`,
results are bit-for-bit identical regardless of `nthreads`.

## References

Friedman, J. H. (1991). Multivariate Adaptive Regression Splines.
*Annals of Statistics* 19(1):1-67.

Friedman, J. H. (1993). *Fast MARS*. Stanford University Department of
Statistics Technical Report 110.

## Examples

``` r
fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
print(fit)
#> Call:
#> ares.default(x = as.matrix(mtcars[, -1]), y = mtcars$mpg, nthreads = 2)
#> 
#> MARS fit:  family = gaussian 
#>   Selected terms: 6 of 21 forward-pass terms
#>   RSS: 91.81 
#>   GCV: 6.662 
#>   Degree: 1   Penalty: 2   nthreads: 2 
p <- predict(fit, as.matrix(mtcars[, -1]))
```
