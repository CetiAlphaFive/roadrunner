# BUG-003 — varmod="lm" silently produces zero-width PIs

**Severity:** correctness (medium).
**File:** `R/predict.R:260-272`.

## Symptom

When `varmod = "lm"` is fitted on a DGP whose residual magnitude
DECREASES with `yhat` (so the lm slope ends up negative), and the model
is then predicted at `yhat` values past the in-sample range, the
predicted `(intercept + slope * yhat)` goes negative. The current code
floors at `1e-12`, so the PI width becomes
`2 * qt(0.975, df) * sqrt(pi/2) * 1e-12 ≈ 4e-12`.

The user sees a (fit, lwr, upr) matrix where `upr - lwr ≈ 4e-12` —
effectively a point estimate disguised as a 95% PI. No warning is emitted.

## Root cause

```r
.ares_make_pi <- function(yhat, object, level = 0.95) {
  vm <- object$varmod
  ...
  sigma_vec <- switch(vm$type,
    const = rep(vm$sigma_hat, length(yhat)),
    lm    = pmax(vm$scale * (vm$intercept + vm$slope * yhat),
                 1e-12)              # floor to keep sigma positive
  )
  cbind(fit = yhat, lwr = yhat - qq * sigma_vec, upr = yhat + qq * sigma_vec)
}
```

The floor at `1e-12` is mathematically necessary to avoid negative SD,
but no warning is raised. When the linear-in-yhat MAD model predicts a
non-positive value, that's a *bug signal*: the model is being asked to
extrapolate beyond its valid range.

## Repro

```r
library(roadrunner)
set.seed(60)
n <- 300; p <- 2
x <- matrix(runif(n * p, 0, 1), n, p)
f <- 10 * x[, 1]
y <- f + rnorm(n, sd = pmax(0.2, 5 - 5 * x[, 1]))   # sigma DECREASES with f
fit <- ares(x, y, varmod = "lm", nthreads = 2)
# fit$varmod: intercept = 4.14, slope = -0.43  -> MAD = 0 at yhat ≈ 9.6

x_big <- matrix(c(seq(0, 3, length.out = 10), rep(0.5, 10)), 10, 2)
pm <- predict(fit, x_big, interval = "pint")
pm[, "upr"] - pm[, "lwr"]   # 7 of 10 rows have width ≈ 3.9e-12
```

Output table (from `/tmp/audit/probe_05c_pi_negative_sigma.R`):

```
        yhat   linear_mad        width
1   8.527920   0.4624170    2.281168e+00
2   3.805652   2.4992808    1.232930e+01
3   6.855161   1.1839308    5.840496e+00
4   9.904670  -0.1314192    3.936407e-12      <-- silent collapse
5  12.954179  -1.4467691    3.936407e-12
...
10 28.201723  -8.0235189    3.936407e-12
```

## Impact

A downstream consumer (e.g. a risk-aware control loop, an uncertainty
quantification dashboard) that uses `predict(..., interval = "pint")` to
make decisions will see a width-zero PI and conclude the model is
EXTREMELY confident at exactly the points where it's actually
extrapolating outside its variance model's training range. This is the
worst possible failure mode for a UQ system.

## Suggested fix

Three options, in order of safety:

1. **Hard error**: stop() when ANY predicted MAD is non-positive, with a
   clear message about extrapolation past the variance model's domain.
2. **Warn + floor at a meaningful lower bound**: warn whenever any
   predicted MAD ≤ 0, and floor at e.g. `min(0.1 * sigma_const, min(
   non-negative linear_mad))`. This at least gives a non-degenerate PI.
3. **Clip predicted yhat to in-sample range before applying the variance
   model**, and warn that the PI at clipped points is the boundary
   estimate. This is essentially what is done in some earth variance
   models.

Option (1) or (2) plus a docstring entry on extrapolation safety would
be ideal.
