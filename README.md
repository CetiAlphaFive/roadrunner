# roadrunner

Fast, low-dependency machine learning algorithms in R.

`roadrunner` ships C++ backed implementations of classical ML algorithms
with thin, base-R-style interfaces. The first algorithm is `ares()`, a
Multivariate Adaptive Regression Splines fitter; more will follow.

## Package design

- **Low dependency.** Only `Rcpp` and `RcppParallel` at runtime. No
  tidyverse, no rlang, no S4.
- **Fast.** C++ engines via `Rcpp` with multi-core scoring via
  `RcppParallel` (TBB).
- **Deterministic.** Fits are bit-for-bit identical across thread
  counts at a fixed seed.
- **Familiar API.** Base-R style: formula and matrix interfaces,
  standard S3 methods (`predict`, `print`, `summary`, `plot`), and
  argument names that mirror existing R packages where one exists.

## Install

```r
# install.packages("pak")
pak::pak("CetiAlphaFive/roadrunner")
```

## Usage

```r
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

See the `ares` vignette for autotune, classification, weights,
bagging, and prediction-interval examples.

## License

MIT (c) Jack T. Rametta
