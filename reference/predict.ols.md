# Predictions from an `ols` fit

Returns predictions for new data from a fitted `ols` model, with
optional standard errors and exact closed-form confidence or prediction
intervals.

## Usage

``` r
# S3 method for class 'ols'
predict(
  object,
  newdata = NULL,
  type = "response",
  se.fit = FALSE,
  interval = c("none", "confidence", "prediction"),
  level = 0.95,
  robust = c("none", "HC0", "HC1", "HC2", "HC3"),
  ...
)
```

## Arguments

- object:

  An object of class `"ols"`.

- newdata:

  A numeric matrix or data frame of new predictors. `NULL` (default)
  returns `object$fitted.values`.

- type:

  Prediction scale. Only `"response"` is supported (the linear predictor
  coincides with the response for a gaussian fit).

- se.fit:

  If `TRUE`, the returned vector carries an `"se.fit"` attribute holding
  the per-row standard error of the fitted mean. For bagged fits the
  standard error is the per-row bag standard deviation. Default `FALSE`.

- interval:

  Type of interval to return.

  - `"none"` (default): point predictions only.

  - `"confidence"`: a matrix with columns `c("fit", "lwr", "upr")` for
    the mean response.

  - `"prediction"`: a matrix with columns `c("fit", "lwr", "upr")` for a
    new observation (adds the residual variance, or the variance-model
    variance when the fit was built with `varmod = "const"` or `"lm"`).

- level:

  Confidence/prediction level. Default `0.95`.

- robust:

  Heteroscedasticity-consistent covariance used for the standard errors
  / intervals.

  - `"none"` (default): the classical `sigma2 (X'WX)^-1`.

  - `"HC0"`, `"HC1"`, `"HC2"`, `"HC3"`: a sandwich estimator.

- ...:

  Currently unused.

## Value

A numeric vector of predictions, or a matrix with columns
`c("fit", "lwr", "upr")` when `interval != "none"`. With `se.fit = TRUE`
the result carries an `"se.fit"` attribute.

## Examples

``` r
fit <- ols(mpg ~ wt + hp, data = mtcars)
predict(fit, mtcars[1:5, ], interval = "confidence")
#>                        fit      lwr      upr
#> Mazda RX4         23.57233 22.45623 24.68843
#> Mazda RX4 Wag     22.58348 21.51622 23.65074
#> Datsun 710        25.27582 23.97441 26.57723
#> Hornet 4 Drive    21.26502 20.10932 22.42072
#> Hornet Sportabout 18.32727 17.30889 19.34564
```
