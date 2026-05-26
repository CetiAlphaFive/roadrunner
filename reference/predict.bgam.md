# Predictions from a `bgam` fit

Returns predictions from a fitted component-wise P-spline boosting
model. New B-spline bases are constructed for `newdata` using the stored
knot sequences; values outside the training range are extrapolated by
the natural boundary behaviour of the de Boor recursion (polynomial
extrapolation). For bagged fits (`n.boot > 0`), the bag mean across all
surviving replicates is returned.

## Usage

``` r
# S3 method for class 'bgam'
predict(
  object,
  newdata = NULL,
  type = c("response", "link", "terms"),
  se.fit = FALSE,
  level = 0.95,
  ...
)
```

## Arguments

- object:

  An object of class `"bgam"`, as returned by
  [`bgam()`](https://cetialphafive.github.io/roadrunner/reference/bgam.md).

- newdata:

  A numeric matrix or data frame of new predictor values with the same
  columns (in any order) as the training data. `NULL` (default) returns
  stored training predictions: `$fitted.values` for non-bagged fits, bag
  mean for bagged fits. Out-of-vocabulary factor levels in `newdata`
  cause an informative error.

- type:

  Prediction scale:

  `"response"`

  :   (default) Fitted probabilities for `family = "binomial"`; fitted
      values (identical to link scale) for `family = "gaussian"`.

  `"link"`

  :   The linear predictor (additive predictor on the link scale).
      Identical to `"response"` for gaussian.

  `"terms"`

  :   An `nrow(newdata)` x `p` numeric matrix of per- predictor
      component contributions. Column `j` holds `B_new_j %*% beta_j`.
      Column names match `object$predictor_names`. Row sums of all
      columns plus `object$intercept_value` equal the link-scale
      prediction.

- se.fit:

  If `TRUE`, returns a list with elements `$fit` (predictions as
  described above) and `$se.fit` (plug-in standard errors summed across
  active base-learners). For predictor j active at observation i, the
  component variance is
  `diag(B_j[i,,drop=FALSE] %*% solve(A_j) %*% t(B_j[i,,drop=FALSE])) * sigma2`
  where `A_j = B_j'B_j + lambda_j * D_j'D_j`. For binomial, delta-method
  SE on the response scale is applied:
  `se_mu = se_eta * mu_hat * (1 - mu_hat)`. Default `FALSE`.

- level:

  Unused; reserved for a future prediction-interval interface. Default
  `0.95`.

- ...:

  Currently unused.

## Value

- `se.fit = FALSE`:

  A numeric vector of length `nrow(newdata)`, or a numeric matrix of
  dimension `nrow(newdata)` x `p` when `type = "terms"`.

- `se.fit = TRUE`:

  A list with elements `$fit` (as above) and `$se.fit` (numeric vector
  of per-observation standard errors).

For bagged fits the bag mean is returned (se.fit is not available for
bagged fits).

## Examples

``` r
set.seed(1)
n <- 200
x <- matrix(rnorm(n * 3), n, 3)
y <- sin(x[, 1]) + 0.5 * x[, 2] + rnorm(n, sd = 0.5)
fit <- bgam(x, y, mstop = 50, autotune = FALSE)
predict(fit, x[1:5, ])
#> [1] -0.32488408  0.77491430 -0.08823158  0.83326512 -0.78510166

# Per-component contributions
terms_mat <- predict(fit, x, type = "terms")
head(terms_mat)
#>              V1         V2           V3
#> [1,] -0.5278110  0.2231641  0.005601494
#> [2,]  0.1220808  0.6562789  0.022393274
#> [3,] -0.6882587  0.6360364 -0.010170570
#> [4,]  1.0240682 -0.1463043 -0.018660156
#> [5,]  0.2259649 -0.9674532 -0.017774634
#> [6,] -0.6785612  0.2965491 -0.019178132

# Plug-in standard errors
p_se <- predict(fit, x[1:5, ], se.fit = TRUE)
p_se$se.fit
#> [1] 0.1363779 0.1451345 0.1600314 0.1360579 0.1144799
```
