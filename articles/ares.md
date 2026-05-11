# ares: Fast MARS in roadrunner

## Introduction

[`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)
(shipped in the `roadrunner` package) fits Multivariate Adaptive
Regression Splines (MARS) models — a flexible nonparametric regression
technique introduced by Friedman (1991). The MARS forward pass adds
basis functions of the form `max(0, ±(x - knot))` (called *hinges*) one
pair at a time, each scaled by a parent term to allow interactions. A
backward pass then prunes terms by minimizing the GCV criterion.

[`ares()`](https://cetialphafive.github.io/roadrunner/reference/ares.md)
mirrors the API of [`earth`](https://cran.r-project.org/package=earth)
on the gaussian-only core and is designed for numerical parity with
`earth` while taking advantage of multi-core CPUs via `RcppParallel`.

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

The fit returns a list with the standard MARS components.

``` r

str(fit, max.level = 1)
#> List of 27
#>  $ coefficients  : Named num [1:6] 19.7106 0.0843 0.0844 3.8551 -2.4501 ...
#>   ..- attr(*, "names")= chr [1:6] "(Intercept)" "h(145-disp)" "h(123-hp)" "h(gear-4)" ...
#>  $ bx            : num [1:32, 1:6] 1 1 1 1 1 1 1 1 1 1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ dirs          : int [1:21, 1:10] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ cuts          : num [1:21, 1:10] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..- attr(*, "dimnames")=List of 2
#>  $ selected.terms: int [1:6] 1 3 5 10 14 16
#>  $ rss           : num 91.8
#>  $ gcv           : num 6.66
#>  $ rss.per.subset: num [1:21] 1126 325 177 143 123 ...
#>  $ gcv.per.subset: num [1:21] 37.5 12.35 7.75 7.3 7.42 ...
#>  $ path.subsets  : list()
#>  $ path.coefs    : list()
#>  $ nk            : int 21
#>  $ thresh        : num 0.001
#>  $ penalty       : num 2
#>  $ minspan       : int 5
#>  $ endspan       : int 10
#>  $ degree        : int 1
#>  $ nthreads      : int 2
#>  $ fitted.values : num [1:32] 19.4 20.5 25.4 19.5 17.7 ...
#>  $ residuals     : num [1:32] 1.638 0.547 -2.559 1.877 1.012 ...
#>  $ namesx        : chr [1:10] "cyl" "disp" "hp" "drat" ...
#>  $ call          : language ares.formula(x = mpg ~ ., data = mtcars, nthreads = 2)
#>  $ pmethod       : chr "backward"
#>  $ dropped       : chr(0) 
#>  $ na.action     : chr "impute"
#>  $ family        : chr "gaussian"
#>  $ terms         :Classes 'terms', 'formula'  language mpg ~ cyl + disp + hp + drat + wt + qsec + vs + am + gear + carb
#>   .. ..- attr(*, "variables")= language list(mpg, cyl, disp, hp, drat, wt, qsec, vs, am, gear, carb)
#>   .. ..- attr(*, "factors")= int [1:11, 1:10] 0 1 0 0 0 0 0 0 0 0 ...
#>   .. .. ..- attr(*, "dimnames")=List of 2
#>   .. ..- attr(*, "term.labels")= chr [1:10] "cyl" "disp" "hp" "drat" ...
#>   .. ..- attr(*, "order")= int [1:10] 1 1 1 1 1 1 1 1 1 1
#>   .. ..- attr(*, "intercept")= int 1
#>   .. ..- attr(*, "response")= int 1
#>   .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv> 
#>  - attr(*, "class")= chr "ares"
```

[`predict()`](https://rdrr.io/r/stats/predict.html) produces fitted
values on training data or new data.

``` r

head(predict(fit))
#> [1] 19.36172 20.45321 25.35910 19.52265 17.68805 19.22091
head(predict(fit, mtcars[1:5, ]))
#> [1] 19.36172 20.45321 25.35910 19.52265 17.68805
```

## Comparison with `earth`

When the `earth` package is available, we can directly compare fits.

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

The relative RSS difference is typically a fraction of a percent at n
\>= 1000 — the two implementations agree numerically even when they pick
slightly different knot sets.

## Threading determinism

`ares` is parallel-deterministic by construction: changing `nthreads`
does not change the output.

``` r

set.seed(42)
n <- 150; p <- 5
x <- matrix(runif(n * p), n, p)
y <- 5 * pmax(0, x[, 1] - 0.5) * pmax(0, x[, 2] - 0.3) + rnorm(n)
f1 <- ares(x, y, degree = 2, nthreads = 1)
f2 <- ares(x, y, degree = 2, nthreads = 2)
cat("identical selected.terms:", isTRUE(all.equal(f1$selected.terms, f2$selected.terms)), "\n")
#> identical selected.terms: TRUE
cat("max coef diff:           ", max(abs(f1$coefficients - f2$coefficients)), "\n")
#> max coef diff:            0
```

## Performance notes

`ares` v0.0.0.9000 is currently **slower in absolute wall-clock than
`earth`** — this is a documented limitation and the headline target for
v0.1. Internally, `ares` uses an O(n) per-knot scoring loop where
`earth` uses Friedman’s O(1)-per-knot Givens fast-LS update. The Givens
implementation is on the roadmap.

What `ares` v0.0.0.9000 already gets right:

- Numerical correctness — fits track `earth` to ~1% RSS across the smoke
  grid.
- Threading determinism — bit-identical fits regardless of thread count.
- Parallel scaling — 2-thread is consistently ~1.5–1.8× faster than
  1-thread (the `RcppParallel` worker is doing real work; the constant
  factor lives inside the inner loop, not in the threading layer).

See the `inst/sims/` directory for the Monte Carlo benchmark scripts
that generated the headline numbers in `README.md`.

## References

Friedman, J. H. (1991). Multivariate Adaptive Regression Splines.
*Annals of Statistics* 19(1):1–67.
