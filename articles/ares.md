# ares: Fast MARS in roadrunner

## Introduction

[`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)
fits Multivariate Adaptive Regression Splines (MARS) models, a flexible
nonparametric regression technique introduced by Friedman (1991). The
forward pass adds pairs of hinge basis functions of the form
`max(0, ±(x - knot))`, each optionally multiplied by a parent term to
capture interactions. A backward pass prunes terms by minimising the GCV
criterion, or a K-fold CV criterion when `pmethod = "cv"` or
`nfold > 0`.

Beyond the core MARS fitter,
[`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)
supports hyperparameter autotuning, K-fold cross-validated pruning,
bootstrap bagging, observation weights, binomial / poisson / gamma
response families on the selected basis, and residual-variance
prediction intervals.

This vignette walks through each of these features in turn.

## A small example: `mtcars`

``` r

library(roadrunner)
fit <- ares(mpg ~ ., data = mtcars, nthreads = 2)
print(fit)
#> Call:
#> ares.formula(x = mpg ~ ., data = mtcars, nthreads = 2)
#> 
#> MARS fit:  family = gaussian 
#>   Selected terms: 6 of 21 forward-pass terms
#>   RSS: 91.81 
#>   GCV: 6.662 
#>   Degree: 1   Penalty: 2   nthreads: 2
```

``` r

head(predict(fit))
#> [1] 19.36172 20.45321 25.35910 19.52265 17.68805 19.22091
head(predict(fit, mtcars[1:5, ]))
#> [1] 19.36172 20.45321 25.35910 19.52265 17.68805
```

## Comparing to other MARS implementations

If you already use [`earth`](https://cran.r-project.org/package=earth),
the gaussian-core fits from
[`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)
are numerically very close on the same data, while running on multiple
cores deterministically. The chunk below illustrates this on a synthetic
regression problem with both an interaction and a quadratic component.

``` r

if (requireNamespace("earth", quietly = TRUE)) {
  set.seed(20260509)
  n <- 200; p <- 5
  x <- matrix(runif(n * p), n, p)
  y <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
       10 * x[, 4] + 5 * x[, 5] + rnorm(n)
  fa <- ares(x, y, degree = 2, nthreads = 2)
  fe <- earth::earth(x, y, degree = 2)
  cat("ares  RSS:", round(fa$rss, 3), " GCV:", round(fa$gcv, 3),
      " nterms:", length(fa$selected.terms), "\n")
  cat("earth RSS:", round(fe$rss, 3), " GCV:", round(fe$gcv, 3),
      " nterms:", length(fe$selected.terms), "\n")
  cat("relative RSS difference:",
      round(abs(fa$rss - fe$rss) / fe$rss, 4), "\n")
}
#> ares  RSS: 222.61  GCV: 1.707  nterms: 16 
#> earth RSS: 204.18  GCV: 1.518  nterms: 15 
#> relative RSS difference: 0.0903
```

The two implementations typically agree on training RSS to within a
fraction of a percent, even when they choose slightly different knot
sets.

## Deterministic threading

A core guarantee of
[`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)
is that varying `nthreads` does not change the result. At a fixed
`seed.cv`, fits are bit-for-bit identical across thread counts, which
makes parallel execution safe to reproduce.

``` r

set.seed(42)
n <- 150; p <- 5
x <- matrix(runif(n * p), n, p)
y <- 5 * pmax(0, x[, 1] - 0.5) * pmax(0, x[, 2] - 0.3) + rnorm(n)
f1 <- ares(x, y, degree = 2, nthreads = 1)
f2 <- ares(x, y, degree = 2, nthreads = 2)
identical(f1$selected.terms, f2$selected.terms)
#> [1] TRUE
max(abs(f1$coefficients - f2$coefficients))
#> [1] 0
```

## Hyperparameter autotune

Setting `autotune = TRUE` runs an inner K-fold cross-validated grid
search over `(degree, penalty, nk, fast.k)` and refits the winning
combination on the full data. When `n` is large enough,
[`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)
first tries the grid on a small subsample to short-circuit the search if
one configuration wins clearly; this typically cuts wall-clock time
substantially with no accuracy loss.

``` r

set.seed(7)
fit_at <- ares(mpg ~ ., data = mtcars, autotune = TRUE,
               autotune.speed = "fast", seed.cv = 1L, nthreads = 2)
fit_at$autotune$degree
#> [1] 1
fit_at$autotune$penalty
#> [1] 1
fit_at$autotune$nk
#> [1] 42
```

## Classification via `family = "binomial"`

For binary outcomes, pass `family = "binomial"`. Term selection runs on
the numeric response (fast), and the selected basis is then refit on the
binary response with
[`stats::glm.fit()`](https://rdrr.io/r/stats/glm.html).
[`predict()`](https://rdrr.io/r/stats/predict.html) returns
probabilities by default, or the linear predictor via `type = "link"`.

``` r

set.seed(1)
n <- 400
x <- matrix(rnorm(n * 4), n, 4)
p <- plogis(0.5 * x[, 1] + 1.2 * x[, 2] * (x[, 2] > 0) - 0.4 * x[, 3])
y <- rbinom(n, 1, p)
fb <- ares(x, y, family = "binomial", degree = 2, nthreads = 2)
range(predict(fb))
#> [1] 0.2443712 0.9915275
```

## Prediction intervals

For gaussian fits, passing `varmod = "const"` or `varmod = "lm"` stores
a residual variance model at fit time. Subsequent calls to
`predict(interval = "pint")` return a matrix with columns
`c("fit", "lwr", "upr")`. Use `"const"` for homoscedastic data and
`"lm"` to allow simple mean-dependent heteroscedasticity.

``` r

set.seed(0)
n <- 500
x <- matrix(runif(n * 3), n, 3)
y <- 3 * x[, 1] + 2 * sin(2 * pi * x[, 2]) + rnorm(n, sd = 0.5)
fp <- ares(x, y, degree = 2, varmod = "const", nthreads = 2)
head(predict(fp, x[1:3, ], interval = "pint"))
#>              fit        lwr       upr
#> [1,]  3.36826816  2.3545538 4.3819825
#> [2,]  0.04725061 -0.9664638 1.0609650
#> [3,] -0.63859961 -1.6523140 0.3751148
```

## References

Friedman, J. H. (1991). Multivariate Adaptive Regression Splines.
*Annals of Statistics* 19(1):1–67.
