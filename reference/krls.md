# Kernel Regularized Least Squares

Fits a Kernel Regularized Least Squares model with Gaussian kernel,
selecting the ridge penalty by leave-one-out cross-validation via the
closed-form identity of Hainmueller and Hazlett (2014). Marginal effects
and their variances are computed by default and returned on the original
`(X, y)` scale.

## Usage

``` r
krls(
  X,
  y,
  sigma = NULL,
  lambda = NULL,
  derivative = TRUE,
  binary = TRUE,
  vcov = TRUE,
  L = NULL,
  U = NULL,
  tol = NULL,
  eigtrunc = NULL,
  print.level = 0
)

# S3 method for class 'krls_rr'
predict(object, newdata, se.fit = FALSE, ...)
```

## Arguments

- X:

  A numeric matrix of predictors (`n x p`). Constant columns are
  rejected. Missing values are not allowed.

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

- print.level:

  Integer. `0` is silent, `1` prints the chosen lambda and
  marginal-effect summaries. Default `0`.

- object:

  A fitted `"krls"` object.

- newdata:

  A numeric matrix with the same columns as the training `X`.

- se.fit:

  Logical. If `TRUE`, return pointwise standard errors of the
  predictions. Requires the fit was created with `vcov = TRUE`.

- ...:

  Currently unused.

## Value

An object of S3 class `"krls"` with components mirroring
[`KRLS::krls()`](https://rdrr.io/pkg/KRLS/man/krls.html): `K`, `coeffs`,
`Looe`, `fitted`, `X`, `y`, `sigma`, `lambda`, `R2`, `derivatives`,
`avgderivatives`, `var.avgderivatives`, `vcov.c`, `vcov.fitted`,
`binaryindicator`.

`Looe` follows the
[`KRLS::krls()`](https://rdrr.io/pkg/KRLS/man/krls.html) scale
convention: it is the sum of squared leave-one-out residuals on the
*standardised* `y` scale multiplied by `sd(y)`, so its units are
`[y^2 / sd_y] = [y]`. This is preserved for downstream compatibility
with code that consumed
[`KRLS::krls()`](https://rdrr.io/pkg/KRLS/man/krls.html) output; it is
**not** the LOO MSE in raw-`y` squared units.

## Details

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

## Predictions on new data with pointwise SEs.
Xnew <- matrix(rnorm(20 * 3), 20, 3)
colnames(Xnew) <- colnames(X)
pr <- predict(fit, Xnew, se.fit = TRUE)
head(pr$fit)
#>            [,1]
#> [1,]  1.8514244
#> [2,]  0.6177231
#> [3,] -0.4742385
#> [4,]  0.9832079
#> [5,] -0.2015325
#> [6,]  0.6944200
head(pr$se.fit)
#>            [,1]
#> [1,] 0.09677354
#> [2,] 0.12094270
#> [3,] 0.05605460
#> [4,] 0.09838911
#> [5,] 0.07140673
#> [6,] 0.14824376
```
