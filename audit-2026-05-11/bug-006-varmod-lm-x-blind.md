# BUG-006 — varmod="lm" cannot capture x-driven heteroscedasticity

**Severity:** robustness (low-medium). Documented limitation, but the
docstring overpromises.

**File:** `R/ares.R:1020-1046`, docstring at `R/ares.R:127-134`.

## Symptom

When the true residual scale `sigma(x)` depends on a predictor whose
contribution to `yhat` is small (e.g. through interactions with other
predictors, or through a covariate that affects variance without
affecting mean), `varmod = "lm"` collapses to roughly the same as
`varmod = "const"`. Coverage degrades catastrophically in high-variance
regions.

## Root cause

`varmod = "lm"` regresses `|resid| ~ yhat` and predicts `sigma(x)` as
`sqrt(pi/2) * (intercept + slope * yhat(x))`. This captures
heteroscedasticity ONLY through the linear-in-yhat channel. If
sigma(x) varies with x in a way that's orthogonal to yhat, the slope
ends up near zero.

## Repro

```r
library(roadrunner)
set.seed(101)
n <- 800; p <- 5
x_tr <- matrix(runif(n * p), n, p)
f_tr <- 10 * sin(pi * x_tr[, 1] * x_tr[, 2]) +
        20 * (x_tr[, 3] - 0.5)^2 + 10 * x_tr[, 4] + 5 * x_tr[, 5]
y_tr <- f_tr + rnorm(n, sd = 1 + 4 * x_tr[, 1])   # sigma is x1-driven

fit <- ares(x_tr, y_tr, varmod = "lm", nthreads = 2)
# fit$varmod$slope ~ 0.029  (essentially flat)

# OOS coverage stratified by x1:
set.seed(202)
x_te <- matrix(runif(4000 * p), 4000, p)
y_te <- ... # same DGP
pm <- predict(fit, x_te, interval = "pint", level = 0.95)
hit <- y_te >= pm[, "lwr"] & y_te <= pm[, "upr"]
bins <- cut(x_te[, 1], breaks = seq(0, 1, by = 0.2), include.lowest = TRUE)
tapply(hit, bins, mean)
#   [0,0.2] (0.2,0.4] (0.4,0.6] (0.6,0.8]  (0.8,1]
#     0.996     0.982     0.951     0.886    0.794   <-- catastrophic at x1=high
```

## Impact

Users following the docstring at `R/ares.R:127-134`:

> `"lm"`: fits a small linear model of |resid| on yhat to allow simple
> mean-dependent heteroscedasticity. Captures simple heteroscedasticity;

will reasonably expect that calling `varmod = "lm"` gives them "simple
heteroscedasticity" coverage. In reality, only **yhat-dependent** scale
shifts are captured; **x-dependent** ones are missed entirely.

## Suggested fix

Three options:

1. **Add `varmod = "x"`** (like earth's varmod.method = "x"): regress
   |resid| via a small MARS fit on x. ~50 lines of R, no engine changes.
2. **Tighten the docstring** to say "captures yhat-dependent
   heteroscedasticity ONLY; x-driven heteroscedasticity is NOT captured
   by varmod='lm' — consider transforming the heteroscedastic predictor
   into yhat through e.g. its interaction with a strong mean predictor".
3. **Warn at fit time** if `cor(|resid|, yhat)` is small (e.g. < 0.1)
   while the residual sd varies materially: "varmod='lm' detected near
   zero slope on |resid| ~ yhat; if you expected heteroscedasticity,
   verify that it's correlated with the fitted values, not with x
   directly."

Option (2) is the cheapest and addresses the spec mismatch. Option (1)
is the proper fix.
