# roadrunner [![roadrunner hex sticker](reference/figures/roadrunner_hex.png)](https://github.com/CetiAlphaFive/roadrunner/blob/main/man/figures/roadrunner_hex.png)

Fast, low-dependency machine learning algorithms in R.

`roadrunner` ships C++ backed implementations of classical ML algorithms
with thin, base-R-style interfaces. Two algorithms today:

- **[`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)**
  – Multivariate Adaptive Regression Splines (MARS).
- **[`krls()`](https://cetialphafive.github.io/roadrunner/reference/krls.md)**
  – Kernel Regularized Least Squares (KRLS).

## Package design

- **Low dependency.** Only `Rcpp`, `RcppArmadillo`, and `RcppParallel`
  at runtime. No tidyverse, no rlang, no S4.
- **Fast.** C++ engines via `Rcpp` with multi-core scoring via
  `RcppParallel`.
- **Deterministic.** Fits are bit-for-bit identical across thread counts
  at a fixed seed.
- **Familiar API.** Base-R style: formula and matrix interfaces and
  standard S3 methods (`predict`, `print`, `summary`, `plot`).

## Install

``` r

# install.packages("pak")
pak::pak("CetiAlphaFive/roadrunner")
```

## Usage

### MARS via `ares()`

``` r

library(roadrunner)

# Regression
fit <- ares(mpg ~ ., data = mtcars, degree = 2)
predict(fit, head(mtcars))

# Classification
y <- as.integer(mtcars$am)
fitc <- ares(as.matrix(mtcars[, -9]), y, family = "binomial")
predict(fitc)

# Hands-free hyperparameter tuning
fitt <- ares(mpg ~ ., data = mtcars, autotune = TRUE,
             autotune.speed = "fast", seed.cv = 1L)
fitt$autotune$degree
```

See the `ares` vignette for autotune, classification, weights, bagging,
and prediction-interval examples.

### KRLS via `krls()`

[`krls()`](https://cetialphafive.github.io/roadrunner/reference/krls.md)
fits a Kernel Regularized Least Squares model (Hainmueller and Hazlett
2014) with closed-form leave-one-out selection of the ridge penalty and
per-observation marginal effects.

``` r

set.seed(1)
n <- 200
X <- matrix(rnorm(n * 3), n, 3)
y <- sin(X[, 1]) + 0.5 * X[, 2]^2 - 0.3 * X[, 3] + rnorm(n, sd = 0.2)

fit <- krls(X, y)             # default sigma = ncol(X); lambda by LOO
fit$avgderivatives            # average marginal effects per variable
predict(fit, X)$fit           # in-sample fitted values

# Predict on new data with pointwise SEs
Xnew <- matrix(rnorm(20 * 3), 20, 3)
pr <- predict(fit, Xnew, se.fit = TRUE)
head(pr$fit); head(pr$se.fit)
```

[`krls()`](https://cetialphafive.github.io/roadrunner/reference/krls.md)
mirrors [`KRLS::krls()`](https://rdrr.io/pkg/KRLS/man/krls.html)
numerically (fits agree to ~1e-13 at matched `sigma` and `lambda`) and
is 6-10x faster on benchmarks at `n >= 500`.

## References

- Friedman, J. H. (1991). Multivariate Adaptive Regression Splines.
  *Annals of Statistics*, 19(1):1-67.
- Hainmueller, J. and Hazlett, C. (2014). Kernel Regularized Least
  Squares: Reducing Misspecification Bias with a Flexible and
  Interpretable Machine Learning Approach. *Political Analysis*,
  22(2):143-168.

## License

MIT (c) Jack T. Rametta
