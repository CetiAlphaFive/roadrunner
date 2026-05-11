# roadrunner

Fast, low-dependency machine learning algorithms in R. First algorithm:
[`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)
— Multivariate Adaptive Regression Splines (Friedman 1991) via `Rcpp` +
`RcppParallel`. Numerically comparable to
[`earth`](https://cran.r-project.org/package=earth) on the gaussian
core, faster on most cells at 4 threads, and parallel-deterministic
(`nthreads = 1` is bit-for-bit identical to `nthreads = N`).

## Install

``` r

# install.packages("devtools")
devtools::install_github("CetiAlphaFive/roadrunner")
```

## Features

- Gaussian / binomial / poisson / gamma response families.
- Observation `weights`.
- Backward (GCV) or `pmethod = "cv"` K-fold cross-validated pruning.
- `autotune = TRUE` — inner-CV grid search over
  `(degree, penalty, nk, fast.k)` with warm-start + successive halving.
- Row-bootstrap `n.boot` bagging with `se.fit`.
- Prediction intervals via `varmod = "const" | "lm"`.
- Built-in NA handling (median impute / drop).
- Factor / character columns expanded via `model.matrix`.

## Usage

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

See the [`ares`
vignette](https://cetialphafive.github.io/roadrunner/vignettes/ares.Rmd)
for autotune, classification, weights, bagging, and prediction-interval
examples.

## License

MIT (c) Jack T. Rametta
