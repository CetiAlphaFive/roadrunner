# Kernel Regularized Least Squares

Fits a Kernel Regularized Least Squares model with Gaussian kernel,
selecting the ridge penalty by leave-one-out cross-validation via the
closed-form identity of Hainmueller and Hazlett (2014). Marginal effects
and their variances are computed by default and returned on the original
`(X, y)` scale.

## Usage

``` r
krls(X, ...)

# S3 method for class 'formula'
krls(X, data = NULL, subset = NULL, ..., y = NULL)

# Default S3 method
krls(
  X,
  y,
  sigma = NULL,
  lambda = NULL,
  derivative = TRUE,
  binary = TRUE,
  vcov = TRUE,
  weights = NULL,
  L = NULL,
  U = NULL,
  tol = NULL,
  eigtrunc = NULL,
  lambda.method = c("loo", "gcv", "cv"),
  lambda.grid = NULL,
  nfold = 0L,
  ncross = NULL,
  stratify = TRUE,
  seed.cv = NULL,
  cv.1se = FALSE,
  autotune = FALSE,
  autotune.grid = NULL,
  autotune.nthreads = NULL,
  varmod = c("none", "const"),
  n.boot = 0L,
  na.action = c("impute", "omit"),
  approx = c("exact", "nystrom"),
  nystrom_m = NULL,
  landmarks = NULL,
  landmark_method = c("random", "kmeans"),
  landmark_seed = NULL,
  nystrom_eps = 1e-09,
  trace = NULL,
  nthreads = 0L,
  print.level = NULL,
  ...
)

# S3 method for class 'krls_rr'
predict(
  object,
  newdata = NULL,
  se.fit = FALSE,
  interval = c("none", "pint"),
  level = 0.95,
  ...
)
```

## Arguments

- X:

  A numeric matrix or data frame of predictors (`n x p`). Constant
  columns are rejected. Factor / character columns in a data frame are
  expanded via `model.matrix(~ ., x)` (treatment contrasts, intercept
  dropped). Missing values are handled via `na.action` (see below).

- ...:

  Currently unused (caught for forward compatibility).

- data:

  Used only by the formula method. A data frame containing the variables
  referenced by the formula.

- subset:

  Used only by the formula method. An optional integer or logical vector
  restricting rows of `data` used for the fit.

- y:

  A numeric response vector or single-column matrix. Constant `y` is
  rejected.

- sigma:

  Gaussian-kernel bandwidth. Default `NULL`, which sets sigma via the
  geomean_p formula: `sqrt(median(d2) * p)` where `d2` are pairwise
  squared Euclidean distances on the standardised predictors and
  `p = ncol(X)`. Must be a positive scalar if supplied.

- lambda:

  Optional ridge penalty. If `NULL` (default), selected by
  golden-section search on the LOO error.

- derivative:

  Logical. If `TRUE` (default), compute pointwise marginal effects and
  their average per variable. Requires `vcov`.

- binary:

  Logical. If `TRUE` (default), columns of `X` with exactly two unique
  values are treated as binary and their marginal effects are replaced
  by predicted-Y first differences (matches
  [`KRLS::fdskrls`](https://rdrr.io/pkg/KRLS/man/fdskrls.html)).

- vcov:

  Logical. If `TRUE` (default), compute the coefficient covariance and
  the variance of average marginal effects.

- weights:

  Optional vector of observation weights (length `n`, strictly
  positive). Internally normalised to mean 1. Implements weighted KRLS
  via a `D K D` transform where `D = diag(sqrt(w))`.
  `weights = rep(1, n)` is byte-identical to the unweighted path.

- L, U:

  Optional lower / upper bracket for the lambda search. If `NULL`,
  defaults follow
  [`KRLS::krls()`](https://rdrr.io/pkg/KRLS/man/krls.html).

- tol:

  Tolerance for the lambda golden section. Default `1e-6` (fixed,
  independent of `n`; was `1e-3 * n` in earlier versions). A fixed small
  tolerance gives 6-digit lambda precision regardless of sample size.

- eigtrunc:

  Optional eigenvalue truncation cutoff in `(0, 1]`. When set,
  eigenvalues below `eigtrunc * max(d)` are dropped from the solve.
  `NULL` (default) keeps all eigenvalues.

- lambda.method:

  Lambda-selection rule. `"loo"` (default) uses the closed-form
  leave-one-out golden-section search; `"cv"` uses K-fold CV over a grid
  (`nfold > 0` required). `"gcv"` uses the closed-form generalised
  cross-validation criterion (Craven & Wahba 1979) using the same
  eigendecomposition; recommended when `n` is large or LOO behaves
  unstably.

- lambda.grid:

  Optional numeric vector of lambda candidates for
  `lambda.method = "cv"`. `NULL` (default) auto-generates a log-spaced
  grid in `[L, U]`.

- nfold:

  Number of CV folds for `lambda.method = "cv"` or for `autotune`. `0`
  (default) disables CV.

- ncross:

  Number of CV repetitions (each builds a fresh fold partition). Default
  `NULL`, which resolves to `2` for the autotune sigma search
  (repeated-CV stabilisation) and `1` for `lambda.method = 'cv'`.
  Explicit integer values are honoured in both paths.

- stratify:

  If `TRUE` (default), CV folds are quantile- stratified on `y`.

- seed.cv:

  Optional integer seed for the CV fold partition.

- cv.1se:

  If `TRUE`, applies the one-standard-error rule when picking lambda
  under CV (smallest model within 1 SE of the minimum mean CV-MSE).
  Default `FALSE`.

- autotune:

  If `TRUE`, runs an inner repeated-CV grid search over `sigma`
  (default: 9-point multiplicative grid centred on the median-heuristic
  sigma anchor, from `anchor * 0.125` to `anchor * 32`) and refits the
  winner on the full data using the one-standard-error rule (selects the
  largest sigma within 1 SE of the minimum CV-MSE). Default `FALSE`.

- autotune.grid:

  Optional numeric vector of `sigma` candidates for autotune. `NULL`
  (default) uses the 9-point anchor-centred grid
  `sigma_anchor * c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32)` where
  `sigma_anchor` is the geomean_p anchor (`sqrt(median(d2) * p)` on the
  standardised predictors).

- autotune.nthreads:

  Integer. Number of worker threads used to evaluate the autotune
  `sigma` grid in parallel. `NULL` (default) reads from
  `getOption("roadrunner.krls.autotune.nthreads")`, and falls back to
  `1L` if unset. Capped at `length(autotune.grid)`. Only used when
  `autotune = TRUE`.

- varmod:

  Residual variance model used to construct prediction intervals via
  `predict(..., interval = "pint")`. `"none"` (default) disables PIs;
  `"const"` uses a homoscedastic `sigma_hat` estimated from the training
  residuals.

- n.boot:

  Number of bootstrap replicates for bagging. `0` (default) disables
  bagging. When `n.boot > 0`, prediction averages across `n.boot`
  replicate fits.

- na.action:

  How to handle missing values in `X`. `"impute"` (default) replaces NAs
  with the column median (stored on the fit and reapplied at
  [`predict()`](https://rdrr.io/r/stats/predict.html) time). `"omit"`
  drops rows with any NA. Missing `y` is always an error.

- approx:

  Either `"exact"` (default) or `"nystrom"`. The exact path runs the
  original O(n^3) eigendecomposition on the full kernel. `"nystrom"`
  replaces it with a low-rank Nystrom approximation anchored at
  `nystrom_m` landmarks, giving O(m^3) + O(n m^2) cost. At default
  `nystrom_m = ceiling(sqrt(n) * 3)` test RMSE is typically within +10%
  of exact on smooth DGPs.

- nystrom_m:

  Integer. Number of landmarks when `approx = "nystrom"`. Default `NULL`
  resolves to `ceiling(sqrt(n) * 3)` (e.g. n=2000 -\> m=135, n=5000 -\>
  m=213). Tuned empirically to keep test-set RMSE within +10% of exact
  on smooth DGPs. Ignored when `landmarks` is supplied explicitly.

- landmarks:

  NULL (default; auto-draw via `landmark_method`), an integer vector of
  row indices into X, or an m x d numeric matrix of landmark coordinates
  in the original X scale (auto-standardized internally using the
  training centers/scales).

- landmark_method:

  Either `"random"` (default; uniform subsample of training rows) or
  `"kmeans"` (Hartigan-Wong centers on standardized X; not bit-stable
  across R versions). Ignored when `landmarks` is supplied.

- landmark_seed:

  Optional integer seed for the local landmark draw. When supplied, uses
  `.with_seed()` to avoid disturbing the caller's global RNG state.
  Required for byte-identical fits across consecutive
  `krls(approx="nystrom")` calls.

- nystrom_eps:

  Numeric. Relative ridge floor applied to the landmark-kernel
  eigenvalues: `D_reg = max(D, nystrom_eps * max(D))`. Defaults to
  `1e-9`. Stabilizes the m x m eigendecomposition when landmarks are
  near-collinear; the fit object's `nystrom_diagnostics$floored_count`
  reports how many eigenvalues hit this floor (useful for tuning m).

- trace:

  Integer. `0` is silent, `> 0` enables progress diagnostics (currently:
  prints chosen lambda when `> 1`, golden-section progress when `> 2`).
  Replaces `print.level`.

- nthreads:

  Integer. Number of threads to use for the C++ kernel build /
  decomposition. `0` (default) means use
  [`RcppParallel::defaultNumThreads()`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html).

- print.level:

  Deprecated alias for `trace`. Prefer `trace`.

- object:

  A fitted `"krls"` object.

- newdata:

  A numeric matrix or data frame with the same columns as the training
  `X` (or matching the training formula). `NULL` returns
  `object$fitted`.

- se.fit:

  Logical. If `TRUE`, return pointwise standard errors of the
  predictions. Requires the fit was created with `vcov = TRUE`.

- interval:

  Prediction-interval mode. `"none"` (default) returns point predictions
  only; `"pint"` returns lower/upper bounds at confidence level `level`,
  requires the fit was created with `varmod = "const"` (or other
  non-`"none"` `varmod`).

- level:

  Confidence level for `interval = "pint"`. Default `0.95`.

## Value

An object of S3 class `c("krls_rr", "krls")` with components mirroring
[`KRLS::krls()`](https://rdrr.io/pkg/KRLS/man/krls.html): `K`, `coeffs`,
`Looe`, `fitted`, `X`, `y`, `sigma`, `lambda`, `R2`, `derivatives`,
`avgderivatives`, `var.avgderivatives`, `vcov.c`, `vcov.fitted`,
`binaryindicator`. Formula-method fits additionally carry `call`,
`terms`, `xlevels`, `factor_info`, `na.action`, and `na.medians` for
downstream [`predict()`](https://rdrr.io/r/stats/predict.html),
[`update()`](https://rdrr.io/r/stats/update.html), and
[`model.matrix()`](https://rdrr.io/r/stats/model.matrix.html) support.

`Looe` follows the
[`KRLS::krls()`](https://rdrr.io/pkg/KRLS/man/krls.html) scale
convention: it is the sum of squared leave-one-out residuals on the
*standardised* `y` scale multiplied by `sd(y)`, so its units are
`[y^2 / sd_y] = [y]`. This is preserved for downstream compatibility
with code that consumed
[`KRLS::krls()`](https://rdrr.io/pkg/KRLS/man/krls.html) output; it is
**not** the LOO MSE in raw-`y` squared units.

## Details

Three call styles are supported (mirrors
[`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)):

- `krls(X, y, ...)` – matrix / numeric interface (back-compatible).

- `krls(y ~ x1 + x2 + ..., data = df, ...)` – formula interface with
  factor expansion and derived terms (`I(x^2)`, `poly(x, 2)`, etc).

- `krls(df, y, ...)` where `df` is a data frame with numeric + factor /
  character columns – expands categorical columns via
  `model.matrix(~ ., df)` (treatment contrasts, intercept dropped).

The numerical pipeline mirrors
[`KRLS::krls()`](https://rdrr.io/pkg/KRLS/man/krls.html) exactly:

1.  `X` and `y` are standardised (column-centred, unit sd).

2.  The Gaussian kernel `K_ij = exp(-||x_i - x_j||^2 / sigma)` is built
    in parallel C++.

3.  `K` is eigendecomposed. All subsequent solves use the eigen-basis
    closed forms, never inverting `K + lambda I` directly.

4.  `lambda` is selected by golden-section search on the closed-form LOO
    error sum, with the same `(L, U, tol)` bracket as
    [`KRLS::krls`](https://rdrr.io/pkg/KRLS/man/krls.html).

5.  Marginal effects use the closed-form identity
    `dy/dx_k = -(2/sigma) * (X_k * (K c) - K diag(c) X)_k`, computed
    without forming the `n x n` distance matrix.

6.  Output is unstandardised back to the original `(X, y)` scale.

At a fixed `(sigma, lambda)`, fits agree with
[`KRLS::krls()`](https://rdrr.io/pkg/KRLS/man/krls.html) to
floating-point precision (typically `< 1e-12` on coefficients, fitted
values, and marginal effects for `n <= 1000`). Wall-clock time is
roughly `6-10x` faster than
[`KRLS::krls()`](https://rdrr.io/pkg/KRLS/man/krls.html) at `n >= 500`
when marginal effects and variance estimates are requested.

Memory scales as `O(n^2)`: the kernel and its squared eigenvector matrix
are both stored. Expect about `0.4 * n^2 / 1e6` MB of peak working
memory (e.g. ~400MB at `n = 1000`, ~10GB at `n = 5000`).

**Scale-aware sigma default (geomean_p anchor, v0.0.0.9041)**

When `sigma = NULL`, roadrunner sets the Gaussian-kernel bandwidth using
the *geomean_p* formula: let `d2_ij = ||Xs_i - Xs_j||^2` be the pairwise
squared Euclidean distances on the standardised predictor matrix `Xs`;
then `sigma = sqrt(median({d2_ij : i < j}) * ncol(Xs))`. This is the
geometric mean of the raw median heuristic and `p = ncol(X)`, and
empirically lands in the sweet spot between the two extremes. At
`n=500, p=10` on standardised N(0,1) data it yields `sigma ~ 13.5`,
versus `~ 18-19` for the raw median and `10` for
[`KRLS::krls()`](https://rdrr.io/pkg/KRLS/man/krls.html)'s fixed
default.

**Why geomean_p**: a 15-DGP head-to-head vs
[`KRLS::krls()`](https://rdrr.io/pkg/KRLS/man/krls.html)
(REQ-20260518-003, iter-2) with the geomean_p anchor shows roadrunner
**wins 12/15 DGPs, loses 0/15, ties 3/15**. The raw median anchor
(v0.0.0.9040) won only 4/15 DGPs because it over-smoothed locally
nonlinear signals (exp-decay, interaction, tanh-interaction, poly2) by
selecting sigma ~2x too wide. The geomean_p anchor preserves
scale-awareness (adapts to the actual data distribution, unlike the
hard-coded `sigma = p`) while avoiding over-smoothing.

For `n > 500` the pairwise distance matrix is `O(n^2)`; to keep the
default cheap, a 500-row subsample is drawn with a fixed seed
(`set.seed(2718)`) so the result is deterministic within a session.

**Autotune-equals-default equivalence at well-fit settings**

When autotune is enabled at settings where the default sigma is already
near- optimal (e.g. additive and interaction DGPs at `n=500, p=10`), the
autotune grid is centred on `sigma_anchor` and the CV argmin coincides
with the grid centre. The 1-SE rule then selects `sigma_anchor`,
producing a fit identical to the non-autotuned default. This is the
expected behaviour — it confirms that the default sigma is well-chosen
for these DGPs — not a failure of autotune. Autotune yields improvements
when the optimal sigma departs from the anchor (e.g. sparse / linear
DGPs where a much wider kernel is better).

Since v0.0.0.9042 the autotune inner loop is parallelised over sigma
candidates using `RcppParallel` (TBB). The pairwise squared-distance
matrix is computed once per CV fold and reused across every sigma in the
grid (`exp(-D / sigma)` is then a cheap elementwise op). Determinism is
preserved: each worker writes to a unique slot in the output vectors, so
the result is byte-identical regardless of `autotune.nthreads`. Pass
`autotune.nthreads = 1` to force strictly sequential execution.

## Note

At a fixed `(sigma, lambda)` fits remain byte-identical to all earlier
versions of roadrunner KRLS. The four defaults changed in v0.0.0.9040
(`sigma`, `tol`, autotune nfold/ncross, autotune grid) only affect
results when those arguments are left at `NULL` / at their default.
Restore any old default explicitly: `sigma = ncol(X)`,
`tol = 1e-3 * nrow(X)`, `nfold = 5`, `ncross = 1`,
`autotune.grid = ncol(X) * c(0.25, 0.5, 1, 2, 4, 8)`.

The sigma anchor formula was refined in v0.0.0.9041 from the raw median
(`median(d2)`) to the geomean_p formula (`sqrt(median(d2) * p)`). To
restore the v0.0.0.9040 median anchor, pass
`sigma = stats::median(as.numeric(stats::dist(scale(X)))^2)`.

Phase 1 speedup (v0.0.0.9042): parallel autotune assumes single-threaded
BLAS for the small kernel operations inside the parallel region. If you
have set a high process-global BLAS thread count (e.g. via
`RhpcBLASctl::blas_set_num_threads(8)`), consider resetting to 1 around
`krls(..., autotune = TRUE)` calls to avoid oversubscription. The
BLAS-heavy distance computation runs OUTSIDE the parallel region and
benefits from multi-threaded BLAS.

`approx = "nystrom"` is currently incompatible with observation
`weights` and with user-supplied lambda-search bracket `L`/`U`/`tol`.
Both error with a clear message at fit time. `vcov = TRUE` is supported
in the single-fit Nystrom path but the resulting `vcov` is for the
m-length dual coefficients, not the n-length kernel coefficients. Phase
3 may extend these.

Since v0.0.0.9043 `krls()` supports an opt-in Nystrom low-rank
approximation via `approx = "nystrom"`. Replacing the full n x n
eigendecomposition with an m x m one (m = `nystrom_m`, default
`ceiling(sqrt(n) * 3)`) gives ~5x speedup at n=2000 and ~10x at n=5000
with RMSE typically within +10% of the exact fit on smooth DGPs.
Landmarks default to a uniform random subsample of training rows; pass
`landmark_method = "kmeans"` for centroidal landmarks. Fits are
byte-identical at fixed `landmark_seed` and `seed.cv`.

## References

Hainmueller, J. and C. Hazlett (2014). "Kernel Regularized Least
Squares: Reducing Misspecification Bias with a Flexible and
Interpretable Machine Learning Approach." *Political Analysis*
22(2):143–168.

## Examples

``` r
set.seed(1)
n <- 100
X <- matrix(rnorm(n * 3), n, 3)
colnames(X) <- c("age", "income", "score")
y <- sin(X[, 1]) + 0.5 * X[, 2]^2 - 0.3 * X[, 3] +
  rnorm(n, sd = 0.2)

fit <- krls(X, y)
fit
#> Kernel Regularized Least Squares (KRLS)
#>   n = 100   p = 3 
#>   sigma = 3.813   lambda = 0.222   R^2 = 0.9563 
#> 
#> Average Marginal Effects:
#>        age     income      score 
#>  0.6021180 -0.1008827 -0.2742915 
fit$avgderivatives           # average marginal effect per variable
#>           age     income      score
#> [1,] 0.602118 -0.1008827 -0.2742915
summary(fit)
#> Kernel Regularized Least Squares (KRLS)
#>   n = 100  p = 3  sigma = 3.813  lambda = 0.222  R^2 = 0.9563
#> 
#> Average Marginal Effects:
#>         Estimate  Std. Err  t value  Pr(>|t|)    
#> age     0.602118  0.022834  26.3690 < 2.2e-16 ***
#> income -0.100883  0.021381  -4.7182 7.822e-06 ***
#> score  -0.274292  0.020456 -13.4086 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Quartiles of Pointwise Marginal Effects:
#>           age     income      score
#> 25% 0.3618781 -0.6924617 -0.3714887
#> 50% 0.7342455 -0.2581948 -0.3193604
#> 75% 0.8978449  0.5280357 -0.2050795

## Formula interface with mixed numeric + factor predictors.
df <- data.frame(age = X[, 1], income = X[, 2],
                 region = factor(sample(c("N", "S"), n, replace = TRUE)),
                 y = y)
fit_f <- krls(y ~ age + income + region, data = df)

## Predictions on new data with pointwise SEs.
Xnew <- matrix(rnorm(20 * 3), 20, 3)
colnames(Xnew) <- colnames(X)
pr <- predict(fit, Xnew, se.fit = TRUE)
head(pr$fit)
#>             [,1]
#> [1,]  1.62457968
#> [2,]  0.01886781
#> [3,] -0.36472318
#> [4,] -0.94094352
#> [5,]  0.09923114
#> [6,]  1.15424173
head(pr$se.fit)
#>            [,1]
#> [1,] 0.08130310
#> [2,] 0.11459798
#> [3,] 0.07798571
#> [4,] 0.09297983
#> [5,] 0.11479828
#> [6,] 0.09957523
```
