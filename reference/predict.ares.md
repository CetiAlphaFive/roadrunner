# Predictions from an `ares` fit

Predictions from an `ares` fit

## Usage

``` r
# S3 method for class 'ares'
predict(
  object,
  newdata = NULL,
  type = c("response", "link"),
  se.fit = FALSE,
  interval = c("none", "pint"),
  level = 0.95,
  ...
)
```

## Arguments

- object:

  An object of class `"ares"`.

- newdata:

  A numeric matrix or data frame of new predictors. If `NULL` (default),
  returns `object$fitted.values`.

- type:

  Prediction scale. `"response"` (default) returns `plogis(eta)` for
  binomial fits, `exp(eta)` for poisson / gamma, and the fitted mean for
  gaussian fits; `"link"` returns the linear predictor `eta`. Ignored
  for gaussian models (link and response coincide).

- se.fit:

  If `TRUE` and `object` was fit with `n.boot > 0`, the returned vector
  carries an `"sd"` attribute holding the per-prediction bag standard
  deviation. Default `FALSE` (plain numeric vector).

- interval:

  Type of interval to return. `"none"` (default) returns only the fitted
  mean. `"pint"` returns a matrix with columns `c("fit", "lwr", "upr")`
  using the variance model stored at fit time (`varmod = "const"` or
  `"lm"`; gaussian only). If `interval = "pint"` but no variance model
  was stored, the call errors with an informative message.

- level:

  Confidence level for prediction intervals when `interval = "pint"`.
  Default `0.95`.

- ...:

  Additional arguments. Currently unused.

## Value

A numeric vector of length `nrow(newdata)`, or a matrix with
`c("fit", "lwr", "upr")` columns when `interval = "pint"`. When
`se.fit = TRUE` and bag fits are present, the result carries
`attr(., "sd")` with the per-row sample SD across replicates.

## Examples

``` r
fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
p <- predict(fit, as.matrix(mtcars[, -1]))
```
