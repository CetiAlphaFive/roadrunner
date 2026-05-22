# Binary logistic regression

Fits a binary logistic regression model by iteratively reweighted least
squares (Fisher scoring). The IRLS step is a weighted economy-QR solve
in C++ (Armadillo/LAPACK), mirroring the
[`ols()`](https://cetialphafive.github.io/roadrunner/reference/ols.md)
engine. The response must be binary; it may be supplied as a 0/1 numeric
vector, a logical vector, or a two-level factor.

## Usage

``` r
logreg(x, ...)

# S3 method for class 'formula'
logreg(
  x,
  data = NULL,
  subset = NULL,
  na.action = stats::na.omit,
  weights = NULL,
  ...,
  y = NULL
)

# Default S3 method
logreg(
  x,
  y,
  weights = NULL,
  intercept = TRUE,
  n.boot = 0L,
  seed = NULL,
  maxit = 25L,
  tol = 0.00000001,
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

  Optional non-negative numeric vector of length `nrow(x)` of prior
  weights for weighted logistic regression. Default `NULL` (unweighted).
  Negative, `NA`, `NaN`, or `Inf` weights are rejected; zero weights are
  rejected (drop the row instead).

- y:

  A binary response of length `nrow(x)`: a 0/1 numeric vector, a logical
  vector, or a factor with exactly two levels (the second level is
  treated as the "success" class).

- intercept:

  Logical; include an intercept term. Default `TRUE`. Honoured by the
  default method; the formula method follows the formula's own intercept
  specification.

- n.boot:

  Number of bootstrap replicate fits for bagging. Default `0` (no
  bagging). When `> 0`,
  [`predict()`](https://rdrr.io/r/stats/predict.html) averages the
  fitted probabilities over the replicates plus the central fit.

- seed:

  Optional integer seed for the bagging bootstrap draws. Pass an integer
  for reproducible bagging; `NULL` (default) uses the current RNG state.
  The user's RNG stream is restored on exit when a seed is supplied.

- maxit:

  Maximum number of IRLS iterations. Default `25`.

- tol:

  Convergence tolerance on the relative change in deviance between IRLS
  iterations. Default `1e-8`.

## Value

An object of class `"logreg"`: a list with the fitted `coefficients`,
`fitted.values` (probabilities), `linear.predictors`, deviance and
working `residuals`, `deviance`, `null.deviance`, `aic`, `df.residual`,
`df.null`, the classical variance-covariance matrix `vcov`, the IRLS
`iter` count and `converged` flag, the numerical `rank`, `hatvalues`,
`(X'WX)^-1` as `XtWXinv`, the design matrix `X` and the (numeric 0/1)
response `y`, echoed `weights`, the `call`, and (when requested) a
`$boot` field.

## Details

Two interfaces are provided. The formula method takes a model formula
and a data frame. The default method takes a numeric predictor matrix
`x` and a binary response `y`.

Optional bagging (`n.boot > 0`) refits the model on bootstrap-weight
resamples; [`predict()`](https://rdrr.io/r/stats/predict.html) then
returns the bag-mean probability. Bagging is a serial loop and is
byte-identical across thread counts.

When the IRLS iteration fails to converge within `maxit` steps – the
usual signature of (quasi-)complete separation – the fit is returned
with `converged = FALSE` and a warning is emitted.

## References

McCullagh, P. and Nelder, J. A. (1989). *Generalized Linear Models*, 2nd
ed. Chapman and Hall.

MacKinnon, J. G. and White, H. (1985). Some
heteroskedasticity-consistent covariance matrix estimators with improved
finite sample properties. *Journal of Econometrics* 29(3):305-325.

## Examples

``` r
df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
fit <- logreg(y ~ wt + hp, data = df)
print(fit)
#> Call:
#> logreg.formula(x = y ~ wt + hp, data = df)
#> 
#> binary logistic regression
#> 
#> Coefficients:
#> (Intercept)           wt           hp  
#>    18.86630     -8.08348      0.03626  
#> 
#>   Null deviance:    43.23 on 31 degrees of freedom
#>   Residual deviance:10.06 on 29 degrees of freedom
#>   AIC: 16.06 
#>   IRLS iterations: 8 (converged) 
predict(fit, df[1:3, ], type = "response")
#>     Mazda RX4 Mazda RX4 Wag    Datsun 710 
#>     0.8423355     0.4047825     0.9702408 
```
