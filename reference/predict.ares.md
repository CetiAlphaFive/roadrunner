# Predictions from an `ares` fit

Returns predictions for new data from a fitted `ares` MARS model. For
bagged fits (`n.boot > 0`), returns the mean across the central fit and
the bootstrap replicates; with `se.fit = TRUE`, the per-row bag standard
deviation is attached as an `"sd"` attribute. For gaussian fits built
with a residual variance model (`varmod = "const"` or `"lm"`),
`interval = "pint"` returns prediction intervals.

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

  A numeric matrix or data frame of new predictors. `NULL` (default)
  returns `object$fitted.values`.

- type:

  Prediction scale.

  - `"response"` (default): probabilities for binomial, response-scale
    means for poisson and gamma, fitted values for gaussian.

  - `"link"`: the linear predictor. Coincides with `"response"` for
    gaussian fits.

- se.fit:

  If `TRUE` and the fit was built with `n.boot > 0`, the returned vector
  carries an `"sd"` attribute holding the per-row bag standard
  deviation. Default `FALSE`.

- interval:

  Type of interval to return.

  - `"none"` (default): return point predictions only.

  - `"pint"`: return a matrix with columns `c("fit", "lwr", "upr")`
    using the variance model stored at fit time. Requires
    `family = "gaussian"` and `varmod = "const"` or `"lm"`; errors
    otherwise.

- level:

  Confidence level for prediction intervals when `interval = "pint"`.
  Default `0.95`.

- ...:

  Currently unused.

## Value

A numeric vector of predictions, or a matrix with columns
`c("fit", "lwr", "upr")` when `interval = "pint"`. With `se.fit = TRUE`
on a bagged fit, the vector carries an `"sd"` attribute.

## Examples

``` r
fit <- ares(as.matrix(mtcars[, -1]), mtcars$mpg, nthreads = 2)
p <- predict(fit, as.matrix(mtcars[, -1]))
```
