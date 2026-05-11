# Fast Multivariate Adaptive Regression Splines

Fits a MARS model (Friedman 1991) using a fast least-squares forward
pass with a parallelized knot search and a backward subset-selection
that minimizes the GCV criterion. The implementation aims for numerical
parity with [`earth::earth()`](https://rdrr.io/pkg/earth/man/earth.html)
on the gaussian-only core while taking advantage of multi-core CPUs to
reduce wall-clock fitting time.

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

  A numeric matrix or data frame of predictors, OR a model formula when
  calling the formula method.

- ...:

  Additional arguments. Currently ignored.

- data:

  Used by the formula method only. A data frame.

- y:

  A numeric response vector. Length must equal `nrow(x)`.

- degree:

  Maximum interaction degree. Default 1.

- nk:

  Maximum number of basis terms in the forward pass. Default
  `min(200, max(20, 2 * ncol(x))) + 1`.

- penalty:

  GCV penalty. Default `if (degree > 1) 3 else 2`.

- thresh:

  Forward-pass relative-RSS early-stop threshold. Default 0.001.

- minspan:

  Minimum knot span. 0 (default) selects an automatic value.

- endspan:

  Knot offset from the data ends. 0 (default) selects automatic.

- adjust.endspan:

  Multiplier applied to `endspan` when the candidate hinge would deepen
  an existing interaction (parent term has degree at least 1). **Default
  `1` for gaussian / poisson / gamma, `2` for binomial** (auto-flipped
  at v0.0.0.9026 with the same A/B-test justification as
  `auto.linpreds`). Pass `1` or `2` explicitly to override the family
  default.

- auto.linpreds:

  If `TRUE`, and the best forward-pass knot for a candidate hinge sits
  at the boundary of its variable's eligible range, substitutes a linear
  term (`dirs = 2`) for the hinge pair. **Default `FALSE` for gaussian /
  poisson / gamma, `TRUE` for binomial** (auto- flipped at v0.0.0.9026
  based on the inst/sims/v0.25 A/B test: earth- like settings tie or
  improve AUC on all mlbench binary cells, but regress some gaussian
  regression cells). Pass an explicit value to override the family
  default.

- fast.k:

  Fast-MARS priority cache size (Friedman 1993). The forward pass always
  rescores pairs whose parent was added in the previous step; the top
  `fast.k` of the remaining stale pairs (ranked by age-discounted cached
  score) are also rescored, and the rest contribute their cached
  reduction. `0` disables the cache (every pair is rescored every step).
  Default `10`.

- fast.beta:

  Age-penalty for the fast-MARS priority cache: a stale pair's effective
  score is `cached_rss_red / (1 + fast.beta * age)`. Default `1.0`. Only
  matters when `fast.k > 0`.

- nprune:

  Maximum number of terms after backward pruning (default `nk`).

- pmethod:

  Pruning method: `"backward"` (default), `"none"`, or `"cv"` for K-fold
  cross-validated subset selection. `"cv"` requires `nfold > 0` (or
  implicitly sets `nfold = 5` when not specified).

- nfold:

  Number of CV folds (default `0` = no CV; GCV-based backward
  selection). When `nfold > 0`, `pmethod` is promoted to `"cv"`.

- ncross:

  Number of CV repetitions when `nfold > 0`. Each repetition builds a
  fresh fold partition; per-size mean MSE is averaged across
  `nfold * ncross` evaluations. Default `1`.

- stratify:

  If `TRUE` (default), folds are quantile-stratified on `y` so each fold
  has a similar response distribution. Useful on small `n` / skewed
  responses; cheap.

- seed.cv:

  Optional integer seed for fold partitioning. `NULL` (default) uses the
  current RNG state; pass an integer for reproducible CV partitions.

- cv.1se:

  If `TRUE`, applies the one-standard-error rule for `pmethod="cv"`:
  among sizes whose mean CV-MSE is within one standard error of the
  minimum, pick the smallest. Trades a bit of accuracy for parsimony.
  Default `FALSE` (pick `argmin` of mean CV-MSE, i.e. the size that
  minimises holdout error directly).

- autotune:

  If `TRUE`, runs an inner cross-validated grid search over
  `(degree, penalty)` and refits the winner. The grid is
  `degree in {1, 2}` (extended to `{1, 2, 3}` when `nk` is large enough)
  crossed with `penalty in {0.5*d, 1.0*d, 2.0*d, 3.0*d}` at each degree
  `d`. Inner CV uses `nfold` (default 5 when autotune is on and
  `nfold == 0`) and respects `seed.cv`. The chosen `(degree, penalty)`
  and the full grid CV-MSE are returned in `$autotune`. Default `FALSE`
  (current grf-style fast default path). Subsequent versions extend the
  grid (v0.16: `nk`, v0.17: `autotune.speed`, v0.18: `n.boot`).

- autotune.speed:

  One of `"balanced"` (default), `"quality"`, or `"fast"`. Only
  meaningful when `autotune = TRUE`.

  - `"quality"`: forces `fast.k = 0` for every cell (no fast-MARS
    priority cache; every (parent, var) pair rescored every step).
    Matches v0.10 quality at v0.10 wall-clock cost.

  - `"fast"`: forces `fast.k = 5`. Most aggressive cache; cheapest.

  - `"balanced"` (default): sweeps `fast.k in {10, 25, Inf}` (where
    `Inf` means "no cache") inside the autotune grid and picks the
    smallest `fast.k` whose mean CV-MSE is within 1% of the best. Trades
    a marginal accuracy hit for a marginal speed gain.

- autotune.warmstart:

  If `TRUE` (default; only meaningful when `autotune = TRUE`) and
  `n >= 200`, ares first runs autotune on a 15% subsample (capped at 200
  rows). If the subsample has a *decisive* winner – best-per-degree
  CV-MSE more than 5% below the next-best degree's best cell – that
  `(degree, penalty, nk, fast.k)` is used directly to refit on the full
  data, skipping the full-data grid. Cuts autotune wall-clock by ~5x on
  well-separated DGPs.

- n.boot:

  Number of bootstrap replicate fits for bagging. Default `0` (no
  bagging). When `> 0`, the result holds a list of `n.boot` ares fits in
  `$boot$fits`, each fit on a row-bootstrap sample of the data. The
  augmented [`predict()`](https://rdrr.io/r/stats/predict.html) averages
  predictions across replicates plus the central fit and reports
  per-prediction standard deviation in `attr(predictions, "sd")`.
  Bagging composes with `autotune` – each replicate fits with the
  central fit's chosen `(degree, penalty, nk, fast.k)` so the cost is
  `(n.boot + 1) * fit_cost`, not the autotune grid times `n.boot`.

- na.action:

  Strategy for missing values in `x`. Either `"impute"` (default) or
  `"omit"`. `"impute"` replaces each column's `NA`s with that column's
  median (computed from the non-missing rows) and stores those medians
  on the fit so [`predict()`](https://rdrr.io/r/stats/predict.html) can
  reapply them to future newdata. `"omit"` drops every row that has any
  `NA` in `x`. Both actions emit a warning summarising the affected rows
  / columns. Missing values in `y` are always a hard error – there is no
  sensible way to impute the regression target. `NaN` and `+/-Inf` in
  `x` are also rejected outright regardless of `na.action`.

- family:

  Response family. One of `"gaussian"` (default; numeric `y`, identity
  link, OLS coefficients), `"binomial"` (binary `y`, logit link),
  `"poisson"` (non-negative integer-valued `y`, log link), or `"gamma"`
  (strictly positive `y`, log link). For the non-gaussian families the
  forward and backward passes still run on the original numeric `y` as
  if gaussian (matching earth\\s GLM strategy); after backward pruning
  the selected basis is refit on the response via
  [`stats::glm.fit()`](https://rdrr.io/r/stats/glm.html) with the
  appropriate `family` object, and `$coefficients`, `$fitted.values`,
  and `$linear.predictor` come from that GLM. Bagging, CV pruning,
  autotune, and `weights` compose with the GLM families; inner
  model-selection still uses (weighted) Gaussian MSE on the latent scale
  (cheap and adequate for term selection — only the final coefficients
  use IRLS).

- weights:

  Optional non-negative numeric vector of length `nrow(x)` for weighted
  least squares fitting. Default `NULL` (unweighted OLS). When supplied,
  the engine implements WLS via a `sqrt(weights)` transform of the
  design and response, so forward + backward both consume the weights
  (knot selection, GCV pruning, CV pruning, and autotune are all
  weight-aware). Weights are renormalised internally so that
  `mean(weights) = 1`; this keeps the GCV denominator valid (the
  "effective sample size" equals `n`). NA / NaN / negative weights are
  rejected. With `weights = rep(1, n)` the fit is byte-identical to the
  unweighted path.

- varmod:

  Variance-model strategy used by
  [`predict.ares()`](https://cetialphafive.github.io/roadrunner/reference/predict.ares.md)
  to build prediction intervals. One of:

  - `"none"` (default): no variance model is stored.
    [`predict()`](https://rdrr.io/r/stats/predict.html) returns only the
    conditional mean; `interval = "pint"` is unavailable.

  - `"const"`: store a single training residual standard deviation
    `sigma_hat = sqrt(weighted_RSS / df_residual)`. PIs are
    `yhat +/- qt(0.975, df) * sigma_hat`.

  - `"lm"`: fit a small linear model `|resid| ~ yhat` on the training
    fit and use its prediction (multiplied by sqrt(pi/2) so it estimates
    sigma) at the new yhat. Captures simple heteroscedasticity. Only
    meaningful for `family = "gaussian"`. For binomial / poisson / gamma
    the variance is mean-determined; the argument is silently ignored
    (PIs would require family-specific intervals not implemented here —
    see ?stats::predict.glm for that). Defaults to `"none"` (SEs off by
    default).

- trace:

  Trace level. 0 = silent (default), 1 = forward-pass progress.

- nthreads:

  Number of threads. 0 (default) selects
  [`RcppParallel::defaultNumThreads()`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html).
  CRAN-distributed examples and the bundled vignette cap this at 2.

## Value

An object of class `"ares"` – a list with components `coefficients`,
`bx`, `dirs`, `cuts`, `selected.terms`, `rss`, `gcv`, `rss.per.subset`,
`gcv.per.subset`, `fitted.values`, `residuals`, `namesx`, `call`, plus
echoed control parameters. When `family` is non-gaussian the fit
additionally carries `$family`, `$glm` (a small list with `deviance`,
`null.deviance`, `df.null`, `df.residual`, `aic`, `converged`, `iter`,
plus family-specific fields), `$linear.predictor` (length `n`), and
`$fitted.values` on the response scale: probabilities for binomial
(`plogis(lp)`), counts/rates for poisson (`exp(lp)`), positive means for
gamma (`exp(lp)`).

## Details

Two interfaces are provided. `ares.default()` accepts a numeric
predictor matrix `x` and a numeric response vector `y`. `ares.formula()`
accepts a model formula and a data frame.

## References

Friedman, J. H. (1991). Multivariate Adaptive Regression Splines.
*Annals of Statistics* 19(1):1-67.

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
