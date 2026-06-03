# Penalized linear discriminant analysis

Penalized Fisher's linear discriminant (Witten & Tibshirani 2011) with
L1 or fused-lasso penalties, multi-class support, and built-in
cross-validation.

## Usage

``` r
plda(x, ...)

# Default S3 method
plda(
  x,
  y,
  K = NULL,
  lambda = NULL,
  penalty = c("L1", "fused"),
  lambda2 = NULL,
  autotune = TRUE,
  nfold = 5L,
  lambda_grid = NULL,
  maxit = 100L,
  tol = 1e-06,
  nthreads = 0L,
  ...
)

# S3 method for class 'formula'
plda(formula, data = NULL, ...)
```

## Arguments

- x:

  Numeric predictor matrix, or a formula for the formula interface.

- ...:

  Passed to methods.

- y:

  Factor (or coercible) class label of length `nrow(x)`.

- K:

  Number of discriminant vectors (`<= G-1`). If supplied, `K` is fixed;
  if `NULL` (default), `K` is chosen by cross-validation along with
  `lambda`. When `autotune = FALSE`, `NULL` means `G-1`.

- lambda:

  Magnitude penalty. Required when `autotune = FALSE`.

- penalty:

  `"L1"` (default) or `"fused"`.

- lambda2:

  Fused-lasso difference penalty (used when `penalty = "fused"`).

- autotune:

  If `TRUE` (default), cross-validate `lambda` (and `K`).

- nfold:

  CV folds.

- lambda_grid:

  Optional CV grid.

- maxit:

  MM solver maximum iterations.

- tol:

  MM solver convergence tolerance.

- nthreads:

  Integer. Number of worker threads for the cross-validation autotune
  (`autotune = TRUE`). `0` (the default, taken from
  `getOption("roadrunner.nthreads")` when set) uses
  [`RcppParallel::defaultNumThreads()`](https://rdrr.io/pkg/RcppParallel/man/setThreadOptions.html).
  The CV fold loop is parallelised with fixed per-fold output slots and
  a serial-order reduction, so the fitted discriminants are
  byte-identical regardless of `nthreads`. Ignored when
  `autotune = FALSE` (a single fit is already sequential).

- formula:

  A model formula; response on the left, predictors on the right.

- data:

  A data frame (or environment) containing the formula variables.

## Value

An object of S3 class `"plda"`.
