# Ordinary and weighted least squares

Fits a linear regression model by ordinary least squares, or by weighted
least squares when `weights` are supplied. The solver is a weighted
economy-QR factorisation in C++ (Armadillo/LAPACK); a rank-deficient
design matrix is rejected with a clear error rather than silently
dropping columns.

## Usage

``` r
ols(x, ...)

# S3 method for class 'formula'
ols(
  x,
  data = NULL,
  subset = NULL,
  na.action = stats::na.omit,
  weights = NULL,
  ...,
  y = NULL
)

# Default S3 method
ols(
  x,
  y,
  weights = NULL,
  intercept = TRUE,
  n.boot = 0L,
  seed = NULL,
  varmod = c("none", "const", "lm"),
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

- subset:

  Optional row-subsetting vector for the formula method (integer indices
  or logical), passed through to
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html).
  Ignored by the default method.

- na.action:

  Strategy for missing values, passed to
  [`stats::model.frame()`](https://rdrr.io/r/stats/model.frame.html) for
  the formula method. Default
  [`stats::na.omit()`](https://rdrr.io/r/stats/na.fail.html). For the
  default method, rows with any `NA` in `x` or `y` are dropped with a
  warning.

- weights:

  Optional non-negative numeric vector of length `nrow(x)` for weighted
  least squares. Default `NULL` (unweighted). Negative, `NA`, `NaN`, or
  `Inf` weights are rejected; zero weights are rejected (drop the row
  instead).

- y:

  A numeric response vector with length `nrow(x)`.

- intercept:

  Logical; include an intercept term. Default `TRUE`. Honoured by the
  default method; the formula method follows the formula's own intercept
  specification.

- n.boot:

  Number of bootstrap replicate fits for bagging. Default `0` (no
  bagging). When `> 0`,
  [`predict()`](https://rdrr.io/r/stats/predict.html) averages over the
  replicates plus the central fit and (with `se.fit = TRUE`) returns a
  per-prediction bag standard deviation.

- seed:

  Optional integer seed for the bagging bootstrap draws. Pass an integer
  for reproducible bagging; `NULL` (default) uses the current RNG state.
  The user's RNG stream is restored on exit when a seed is supplied.

- varmod:

  Residual variance model used by
  [`predict()`](https://rdrr.io/r/stats/predict.html) for prediction
  intervals.

  - `"none"` (default): the closed-form prediction interval uses the
    fitted residual standard error directly.

  - `"const"`: stores a single residual SD.

  - `"lm"`: fits a small linear model of `|resid|` on the fitted values
    to allow simple yhat-dependent heteroscedasticity. Shares the
    [`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)
    variance-model helper.

## Value

An object of class `"ols"`: a list with the fitted `coefficients`,
`fitted.values`, `residuals`, `sigma2` and `sigma`, `df.residual`, the
classical variance-covariance matrix `vcov`, the numerical `rank`,
`hatvalues`, `(X'WX)^-1` as `XtXinv`, the design matrix `X` and response
`y`, echoed `weights`, the `call`, and (when requested) `$varmod` and
`$boot` fields.

## Details

Two interfaces are provided. The formula method takes a model formula
and a data frame. The default method takes a numeric predictor matrix
`x` and a numeric response vector `y`.

Optional bagging (`n.boot > 0`) refits the model on bootstrap-weight
resamples; [`predict()`](https://rdrr.io/r/stats/predict.html) then
returns the bag mean. Bagging is a serial loop and is byte-identical
across thread counts.

## References

MacKinnon, J. G. and White, H. (1985). Some
heteroskedasticity-consistent covariance matrix estimators with improved
finite sample properties. *Journal of Econometrics* 29(3):305-325.

## Examples

``` r
fit <- ols(mpg ~ wt + hp, data = mtcars)
print(fit)
#> Call:
#> ols.formula(x = mpg ~ wt + hp, data = mtcars)
#> 
#> linear model (OLS)
#> 
#> Coefficients:
#> (Intercept)           wt           hp  
#>    37.22727     -3.87783     -0.03177  
#> 
#>   Residual SE: 2.593 on 29 degrees of freedom
predict(fit, mtcars[1:3, ])
#>     Mazda RX4 Mazda RX4 Wag    Datsun 710 
#>      23.57233      22.58348      25.27582 
```
