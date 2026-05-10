# ares — Fast Multivariate Adaptive Regression Splines

[![R-CMD-check](https://github.com/jtrametta/ares/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jtrametta/ares/actions/workflows/R-CMD-check.yaml)

**`ares`** fits Multivariate Adaptive Regression Splines (MARS) models in R, using
Friedman's (1991) forward stepwise hinge selection and a backward subset-selection
that minimizes GCV. The implementation is built on
[**Rcpp**](https://www.rcpp.org/) and
[**RcppParallel**](https://rcppcore.github.io/RcppParallel/), and is designed to
produce numerically comparable fits to the well-established
[**`earth`**](https://cran.r-project.org/package=earth) package while taking
advantage of multi-core CPUs.

## Installation

```r
# devtools::install_github("jtrametta/ares")
# or for local development:
# devtools::install("/path/to/ares")
```

## Usage

```r
library(ares)

# Matrix interface
x <- as.matrix(mtcars[, -1])
y <- mtcars$mpg
fit <- ares(x, y, degree = 2, nthreads = 2)
print(fit)
summary(fit)
predict(fit, head(x))

# Formula interface
fit2 <- ares(mpg ~ ., data = mtcars, degree = 2, nthreads = 2)
```

## Status — v0.0.0.9000

This is the initial development release. The core MARS engine is functional and
produces fits that closely match `earth` numerically. Speed parity vs `earth` is a
near-term roadmap item; the current release prioritizes correctness over absolute
wall-clock.

### Numerical parity vs `earth`

Across the smoke-grid Monte Carlo (24 cells × 5 reps = 120 fits per estimator,
spanning Friedman-1, additive, and interaction DGPs at n ∈ {200, 1000} and
p ∈ {5, 10}, degree ∈ {1, 2}):

| Quantity | Median across cells |
| --- | --- |
| `rss_match` (relative RSS difference vs earth) | **1.3%** |
| `gcv_ratio` (ares GCV / earth GCV) | **0.999** |
| `nterms_diff` | **0** |
| Predicted-y RMSE / sd(y) | **7.7%** |

### Wall-clock vs `earth`

| n | p | degree | earth | ares (1t) | ares (2t) | ares-2t / earth |
| --- | --- | --- | --- | --- | --- | --- |
| 500 | 10 | 2 | 0.04 s | 2.96 s | 1.17 s | 33× |
| 1500 | 10 | 2 | 0.09 s | 22.5 s | 6.6 s | 78× |

`ares` v0.0.0.9000 is **slower than `earth` in absolute wall-clock**. Earth uses
hand-tuned C with O(1)-per-knot Givens fast-LS updates; `ares` v0.0.0.9000 uses an
O(n*K*M_q) per-pair scoring loop. The Givens fast-LS inner loop is a **v0.1
milestone**.

What `ares` v0.0.0.9000 already does well:

- **Numerical correctness**: matches earth to ~1% RSS across the smoke grid.
- **Threading determinism**: `nthreads = 1` and `nthreads = 2` produce
  byte-identical fits (selected terms, RSS, coefficients).
- **Parallel scaling**: 2-thread is ~1.75× faster than 1-thread (≈85%
  efficiency) — the `RcppParallel` worker is doing real work.
- **Standard scaffold**: built with `usethis` / `devtools` / `roxygen2` /
  `testthat 3` and passes `R CMD check --as-cran`.

## API

The single user-facing function is `ares()`, with `formula` and `default`
methods. Standard S3 methods: `predict.ares`, `print.ares`, `summary.ares`,
`plot.ares`.

Default arguments mirror `earth`'s gaussian-only core:

| Argument | Default |
| --- | --- |
| `degree` | 1 |
| `nk` | `min(200, max(20, 2*ncol(x))) + 1` |
| `penalty` | 2 (degree=1) or 3 (degree>1) |
| `thresh` | 0.001 |
| `minspan` | 0 (auto, Friedman 1991 formula) |
| `endspan` | 0 (auto) |
| `pmethod` | "backward" |
| `nthreads` | 0 (= `RcppParallel::defaultNumThreads()`) |

The return list mirrors `earth`'s structure: `coefficients`, `bx`, `dirs`,
`cuts`, `selected.terms`, `rss`, `gcv`, `rss.per.subset`, `gcv.per.subset`,
`fitted.values`, `residuals`, plus echoed control parameters.

## References

Friedman, J. H. (1991). Multivariate Adaptive Regression Splines.
*Annals of Statistics* 19(1):1–67.

## License

MIT © Jack Trametta
