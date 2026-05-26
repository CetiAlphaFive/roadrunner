# roadrunner <a href='https://github.com/CetiAlphaFive/roadrunner/blob/main/man/figures/roadrunner_hex.png'><img src='man/figures/roadrunner_hex.png' align="right" height="278" alt="roadrunner hex sticker" /></a>

<!-- badges: start -->
[![R-CMD-check](https://github.com/CetiAlphaFive/roadrunner/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/CetiAlphaFive/roadrunner/actions/workflows/R-CMD-check.yaml)
[![lint](https://github.com/CetiAlphaFive/roadrunner/actions/workflows/lint.yaml/badge.svg)](https://github.com/CetiAlphaFive/roadrunner/actions/workflows/lint.yaml)
[![Codecov test coverage](https://codecov.io/gh/CetiAlphaFive/roadrunner/graph/badge.svg)](https://app.codecov.io/gh/CetiAlphaFive/roadrunner)
<!-- badges: end -->

Fast, low-dependency machine learning algorithms in R. Useful for causal plug-ins (e.g., nuisance fits in DML) or simple predictive applications. 

`roadrunner` ships C++ backed implementations of classical ML algorithms
with thin, base-R-style interfaces. Six core fitters today:

- **`ares()`** -- Multivariate Adaptive Regression Splines (MARS).
- **`krls()`** -- Kernel Regularized Least Squares (KRLS).
- **`plda()`** -- Penalized Linear Discriminant Analysis (L1 / fused-lasso).
- **`ols()`** -- Ordinary and weighted least squares with HC robust SEs.
- **`logreg()`** -- Binary logistic regression by IRLS with HC robust SEs.
- **`bgam()`** -- Component-wise P-spline gradient boosting; smooth additive
  model with built-in variable selection via boosting early stopping.
  Gaussian and binomial families.

Plus **`meep()`** -- a cross-fitted ensemble of `ares()`, `krls()`, `ols()`,
`logreg()`, and opt-in `bgam()`, built for Double Machine Learning and
causal-forest nuisance estimation.

## Package design

- **Low dependency.** Only `Rcpp`, `RcppArmadillo`, and `RcppParallel`
  at runtime. 
- **Fast.** C++ engines via `Rcpp` with multi-core scoring via
  `RcppParallel`.
- **Deterministic.** Fits are bit-for-bit identical across thread
  counts at a fixed seed.
- **Familiar API.** Base-R style: formula and matrix interfaces and
  standard S3 methods (`predict`, `print`, `summary`, `plot`).

## Install

```r
# install.packages("pak")
pak::pak("CetiAlphaFive/roadrunner")
```

## Usage

### MARS via `ares()`

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

### KRLS via `krls()`

`krls()` fits a Kernel Regularized Least Squares model (Hainmueller and
Hazlett 2014) with closed-form leave-one-out selection of the ridge
penalty and per-observation marginal effects.

```r
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

`krls()` mirrors `KRLS::krls()` numerically (fits agree to ~1e-13 at
matched `sigma` and `lambda`) and is 6-10x faster on benchmarks at
`n >= 500`.

### Causal ensembles via `meep()`

`meep()` cross-fits an ensemble of `ares()` and `krls()` and returns
out-of-fold predictions designed to drop into Double Machine Learning `DoubleML` or
causal forest `grf` implementations. 

```r
set.seed(1)
n <- 800
X <- matrix(runif(n * 4, -2, 2), n, 4)
D <- sin(X[, 1]) + 0.5 * X[, 2] + 0.3 * X[, 3]^2 + rnorm(n, sd = 0.3)
Y <- D + cos(X[, 1]) + 0.4 * X[, 2]^2 + 0.5 * X[, 3] + rnorm(n, sd = 0.3)

m  <- meep(X, Y, treatment = D, folds = 5, seed = 42)
m$y_hat_oof   # cross-fitted E[Y | X]
m$d_hat_oof   # cross-fitted E[D | X]

# hand the cross-fitted nuisances to a causal forest
# grf::causal_forest(X, Y, D, Y.hat = m$y_hat_oof, W.hat = m$d_hat_oof)
```

On smooth, structured signal the spline + kernel ensemble fits the
nuisances more tightly than `grf`'s built-in regression forests
(out-of-bag vs out-of-fold R-squared on the toy above):

```r
cf <- grf::causal_forest(X, Y, D, seed = 42)
r2 <- function(p, a) 1 - sum((a - p)^2) / sum((a - mean(a))^2)

data.frame(
  nuisance    = c("E[Y|X]", "E[D|X]"),
  grf_oob_r2  = c(r2(cf$Y.hat, Y),     r2(cf$W.hat, D)),
  meep_oof_r2 = c(r2(m$y_hat_oof, Y),  r2(m$d_hat_oof, D))
)
#>  nuisance grf_oob_r2 meep_oof_r2
#>    E[Y|X]      0.854       0.900
#>    E[D|X]      0.865       0.916
```

### Penalized LDA via `plda()`

`plda()` fits penalized Fisher's linear discriminant analysis (Witten &
Tibshirani 2011) with L1 or fused-lasso penalties, multi-class support,
and built-in cross-validation autotune.

```r
fit <- plda(Species ~ ., data = iris)
predict(fit, iris)        # factor of predicted Species
```

### Linear models via `ols()` and `logreg()`

`ols()` fits ordinary and weighted least squares; `logreg()` fits binary
logistic regression by IRLS. Both have C++ engines, classical and HC0-HC3
robust standard errors, optional bagging, and the standard formula and
matrix interfaces.

```r
fit <- ols(mpg ~ wt + hp, data = mtcars)
summary(fit, robust = "HC3")            # HC3 robust standard errors
predict(fit, mtcars[1:3, ], interval = "confidence")

df  <- data.frame(am = mtcars$am, mtcars[c("wt", "hp")])
lr  <- logreg(am ~ wt + hp, data = df)
predict(lr, df[1:3, ], type = "response")
```

Both also slot into `meep()` as opt-in learners (`learners = c("ols",
"logreg")`).

## References

- Friedman, J. H. (1991). Multivariate Adaptive Regression Splines.
  *Annals of Statistics*, 19(1):1-67.
- Hainmueller, J. and Hazlett, C. (2014). Kernel Regularized Least
  Squares: Reducing Misspecification Bias with a Flexible and
  Interpretable Machine Learning Approach. *Political Analysis*,
  22(2):143-168.
- Witten, D. M. and Tibshirani, R. (2011). Penalized classification using
  Fisher's linear discriminant. *Journal of the Royal Statistical Society,
  Series B*, 73(5):753-772.
- Buehlmann, P. and Yu, B. (2003). Boosting with the L2 loss: Regression
  and classification. *Journal of the American Statistical Association*,
  98(462):324-339.
- Eilers, P. H. C. and Marx, B. D. (1996). Flexible smoothing with
  B-splines and penalties. *Statistical Science*, 11(2):89-121.
- Hofner, B., Mayr, A., Robinzonov, N. and Schmid, M. (2014). Model-based
  boosting in R: A hands-on tutorial using the R package mboost.
  *Computational Statistics*, 29(1-2):3-35.

## License

MIT (c) Jack T. Rametta
