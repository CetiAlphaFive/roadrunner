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
  lambda.method = c("loo", "cv"),
  lambda.grid = NULL,
  nfold = 0L,
  ncross = 1L,
  stratify = TRUE,
  seed.cv = NULL,
  cv.1se = FALSE,
  autotune = FALSE,
  autotune.grid = NULL,
  varmod = c("none", "const"),
  n.boot = 0L,
  na.action = c("impute", "omit"),
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

  Gaussian-kernel bandwidth. Default `ncol(X)`. Must be a positive
  scalar.

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

  Tolerance for the lambda golden section. Default `1e-3 * n`.

- eigtrunc:

  Optional eigenvalue truncation cutoff in `(0, 1]`. When set,
  eigenvalues below `eigtrunc * max(d)` are dropped from the solve.
  `NULL` (default) keeps all eigenvalues.

- lambda.method:

  Lambda-selection rule. `"loo"` (default) uses the closed-form
  leave-one-out golden-section search; `"cv"` uses K-fold CV over a grid
  (`nfold > 0` required).

- lambda.grid:

  Optional numeric vector of lambda candidates for
  `lambda.method = "cv"`. `NULL` (default) auto-generates a log-spaced
  grid in `[L, U]`.

- nfold:

  Number of CV folds for `lambda.method = "cv"` or for `autotune`. `0`
  (default) disables CV.

- ncross:

  Number of CV repetitions (each builds a fresh fold partition). Default
  `1`.

- stratify:

  If `TRUE` (default), CV folds are quantile- stratified on `y`.

- seed.cv:

  Optional integer seed for the CV fold partition.

- cv.1se:

  If `TRUE`, applies the one-standard-error rule when picking lambda
  under CV (smallest model within 1 SE of the minimum mean CV-MSE).
  Default `FALSE`.

- autotune:

  If `TRUE`, runs an inner CV grid search over `sigma` (default grid:
  `ncol(X) * c(0.25, 0.5, 1, 2, 4, 8)`) and refits the winner on the
  full data. Default `FALSE`.

- autotune.grid:

  Optional numeric vector of `sigma` candidates for autotune. `NULL`
  (default) uses the default multiplicative grid.

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
#>   sigma = 3   lambda = 0.07407   R^2 = 0.9735 
#> 
#> Average Marginal Effects:
#>         age      income       score 
#>  0.61762062 -0.09742387 -0.27850374 
fit$avgderivatives           # average marginal effect per variable
#>            age      income      score
#> [1,] 0.6176206 -0.09742387 -0.2785037
summary(fit)
#> Kernel Regularized Least Squares (KRLS)
#>   n = 100  p = 3  sigma = 3  lambda = 0.07407  R^2 = 0.9735
#> 
#> Average Marginal Effects:
#>         Estimate  Std. Err t value  Pr(>|t|)    
#> age     0.617621  0.022220  27.795 < 2.2e-16 ***
#> income -0.097424  0.019662  -4.955 2.988e-06 ***
#> score  -0.278504  0.019242 -14.474 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Quartiles of Pointwise Marginal Effects:
#>           age     income      score
#> 25% 0.4562498 -0.5912371 -0.3714893
#> 50% 0.7490423 -0.1920135 -0.2998952
#> 75% 0.8713291  0.4931363 -0.2067177

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
#> [1,]  1.61265745
#> [2,]  0.11624507
#> [3,] -0.31083262
#> [4,] -0.98481768
#> [5,]  0.05253484
#> [6,]  1.14515989
head(pr$se.fit)
#>            [,1]
#> [1,] 0.08151507
#> [2,] 0.13336955
#> [3,] 0.08682870
#> [4,] 0.09446783
#> [5,] 0.10664122
#> [6,] 0.12060081
```
