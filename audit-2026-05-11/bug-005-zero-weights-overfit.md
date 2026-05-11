# BUG-005 — Zero-weight observations bias GCV downward, causing over-fitting

**Severity:** robustness (medium).
**Files:** `src/ares.cpp:1206-1211` (compute_gcv lambda, both
`mars_fit_cpp` and `mars_backward_only_cpp`), `R/ares.R:1003-1018`
(varmod df).

## Symptom

When some observations have `weights = 0` (and others > 0), the GCV
denominator uses raw row count `n`, not effective sample size
`n_eff = sum(w_i > 0)`. The per-knot penalty is thus implicitly
under-estimated, and the model retains MORE terms than the equivalent
dropped-rows fit.

The same `n`-vs-`n_eff` mismatch propagates to:
- `varmod = "const"` df at `ares.R:1008` (uses raw `n - k`).
- `auto_minspan`/`auto_endspan` defaults at `src/ares.cpp:67-82`
  (use raw `n`).

## Root cause

`compute_gcv` lambda inside `mars_fit_cpp`:

```cpp
auto compute_gcv = [&](double rss_val, int n_terms) {
  double C = double(n_terms) + penalty * (double(n_terms) - 1.0) / 2.0;
  double denom = 1.0 - C / double(n);    // <-- n is raw row count
  if (denom <= 0.0) return std::numeric_limits<double>::infinity();
  return rss_val / (double(n) * denom * denom);
};
```

Under the convention `mean(w) = 1` (enforced by R-side at
`ares.R:541-546`), `sum(w_i) = n`. So with all-positive weights, the
formula is mathematically equivalent to using effective sample size. But
when zero weights are present, `sum(w_i > 0) < n` while `sum(w_i) = n`
still — so `n` is wrong as a count of degrees-of-freedom-bearing
observations.

## Repro

```r
library(roadrunner)
set.seed(120)
n <- 200; p <- 4
x <- matrix(runif(n*p), n, p)
y <- 3*x[, 1] + sin(2*pi*x[, 2]) + rnorm(n, sd = 0.4)
w_z <- runif(n, 0.5, 2)
w_z[1:50] <- 0   # 50 rows have zero weight

fwz <- ares(x, y, weights = w_z, nthreads = 2)
fdr <- ares(x[w_z > 0, ], y[w_z > 0], nthreads = 2)  # equivalent: drop those rows

length(fwz$selected.terms)   # 10
length(fdr$selected.terms)   # 7
```

The 3 extra terms in `fwz` are GCV-pruned out of `fdr` because the
"correct" `n_eff = 150` would penalize them. With `n = 200` the GCV
penalty is too lax.

## Impact

Statistical: over-fitting whenever zero weights are present.
Numerical: `varmod` PI widths use the wrong df, slightly miscalibrating
coverage. Probe 11: sigma_hat 0.4417 vs 0.4300 with df 93 vs 81.

## Suggested fix

Two clean options:

1. **Reject zero weights up front.** Require all weights to be strictly
   positive at `ares.R:337-341`. Document that "to omit a row, drop it
   from x/y; weights = 0 is not a substitute".
2. **Compute and propagate n_eff = sum(w_i > 0)** down to the C++
   engine: pass it as an extra argument to `mars_fit_cpp` and
   `mars_backward_only_cpp`, replace `n` in `compute_gcv`,
   `auto_minspan`, `auto_endspan`. Use it for `varmod` df too.

Option (1) is simpler and lossless (the user can drop rows themselves).
