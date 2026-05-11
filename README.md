# roadrunner — Fast, Low-Dependency Machine Learning Algorithms

[![R-CMD-check](https://github.com/jtrametta/roadrunner/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jtrametta/roadrunner/actions/workflows/R-CMD-check.yaml)

**`roadrunner`** is a collection of fast, low-dependency implementations of
classical machine learning algorithms. The first algorithm shipped is `ares()`,
a Multivariate Adaptive Regression Splines (MARS) fitter using Friedman's (1991)
forward stepwise hinge selection and a GCV-minimizing backward subset-selection.
The implementation is built on
[**Rcpp**](https://www.rcpp.org/) and
[**RcppParallel**](https://rcppcore.github.io/RcppParallel/), and is designed to
produce numerically comparable fits to the well-established
[**`earth`**](https://cran.r-project.org/package=earth) package while taking
advantage of multi-core CPUs. Additional algorithms will be added as separate
exported functions under the same package roof.

## Installation

```r
# devtools::install_github("jtrametta/roadrunner")
# or for local development:
# devtools::install("/path/to/roadrunner")
```

## Usage

```r
library(roadrunner)

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

## Status — v0.0.0.9028

The MARS engine (`ares()`) is the first algorithm shipped in `roadrunner`. It
matches `earth` numerically on the gaussian-only core, is parallel-deterministic
across thread counts, and supports CV pruning, autotune, and binomial
classification.

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

Median speed ratio (ares-2t / earth) across the inst/sims grid: **~0.93×**
(`ares` is slightly faster than `earth` on most cells). Autotune wall-clock is
sub-second on warm-start hits and ~5–15 s on the high-p p=20 grid.

What `ares()` already does well:

- **Numerical correctness**: matches earth to ~1% RSS across the smoke grid.
- **Threading determinism**: `nthreads = 1` and `nthreads = N` produce
  byte-identical fits (selected terms, RSS, coefficients).
- **Parallel scaling**: 2-thread is ~1.75× faster than 1-thread (≈85%
  efficiency).
- **Standard scaffold**: built with `usethis` / `devtools` / `roxygen2` /
  `testthat 3` and passes `R CMD check --as-cran`.

## API

The MARS fitter is exposed as `ares()`, with `formula` and `default`
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
