# Predictions from a `logreg` fit

Returns predictions for new data from a fitted `logreg` model on the
link, probability, or class scale, with optional standard errors.

## Usage

``` r
# S3 method for class 'logreg'
predict(
  object,
  newdata = NULL,
  type = c("link", "response", "class"),
  se.fit = FALSE,
  threshold = 0.5,
  robust = c("none", "HC0", "HC1", "HC2", "HC3"),
  ...
)
```

## Arguments

- object:

  An object of class `"logreg"`.

- newdata:

  A numeric matrix or data frame of new predictors. `NULL` (default)
  returns predictions for the training data.

- type:

  Prediction scale.

  - `"link"`: the linear predictor (log-odds).

  - `"response"` (default): the fitted probability.

  - `"class"`: a hard 0/1 (or factor) label obtained by thresholding the
    fitted probability at `threshold`.

- se.fit:

  If `TRUE`, the returned object carries an `"se.fit"` attribute holding
  the per-row standard error. For `type = "link"` the SE is on the
  log-odds scale; for `type = "response"` it is the delta-method SE of
  the probability. `se.fit` is not available for `type = "class"`. For
  bagged fits the SE is the per-row bag standard deviation of the
  predicted probability. Default `FALSE`.

- threshold:

  Probability cut-off for `type = "class"`. Default `0.5`.

- robust:

  Heteroscedasticity-consistent covariance used for the standard errors.

  - `"none"` (default): the classical maximum-likelihood covariance.

  - `"HC0"`, `"HC1"`, `"HC2"`, `"HC3"`: a sandwich estimator.

- ...:

  Currently unused.

## Value

A numeric vector of predictions on the requested scale (a factor when
`type = "class"` and the fit was built from a factor response). With
`se.fit = TRUE` the result carries an `"se.fit"` attribute.

## Examples

``` r
df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
fit <- logreg(y ~ wt + hp, data = df)
predict(fit, df[1:5, ], type = "response")
#>         Mazda RX4     Mazda RX4 Wag        Datsun 710    Hornet 4 Drive 
#>        0.84233554        0.40478253        0.97024082        0.04172803 
#> Hornet Sportabout 
#>        0.06938812 
predict(fit, df[1:5, ], type = "class")
#> [1] 1 0 1 0 0
```
