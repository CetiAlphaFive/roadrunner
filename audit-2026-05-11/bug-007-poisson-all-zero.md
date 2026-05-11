# BUG-007 — family="poisson" accepts all-zero y silently

**Severity:** robustness (low).
**File:** `R/ares.R:282-296`.

## Symptom

`ares(x, rep(0L, n), family = "poisson")` returns a fit without error.
The fit is degenerate: deviance = 0, fitted means = 0, no signal.

## Root cause

The y-validation block for poisson at ares.R:282-296:

```r
} else if (family == "poisson") {
  ...
  if (!is.numeric(y) && !is.integer(y))
    stop("ares: family = 'poisson' requires numeric y.")
  y <- as.numeric(y)
  yv <- y[is.finite(y)]
  if (any(yv < 0))
    stop("ares: family = 'poisson' requires y >= 0; got negative values.")
  if (any(abs(yv - round(yv)) > 1e-8))
    stop("ares: family = 'poisson' requires integer-valued y; ...")
}
```

`y = c(0, 0, ..., 0)` passes both checks. The forward+backward run as
gaussian, see zero variance everywhere, and emit only the intercept;
then `glm.fit` is called on a 1-column bx with poisson + log-link and
returns a degenerate fit.

## Repro

```r
library(roadrunner)
n <- 100; p <- 3
x <- matrix(runif(n*p), n, p)
y0 <- rep(0L, n)
fit <- ares(x, y0, family = "poisson", nthreads = 2)  # silently succeeds
fit$glm$deviance           # 0
length(fit$selected.terms) # 1 (intercept only)
fit$coefficients           # -Inf (or NaN, replaced by 0 in the engine cleanup)
```

## Impact

Low. The fit is degenerate but doesn't crash. Predict still returns
something. But the user gets no signal that the input is pathological.

## Suggested fix

Add a guard:

```r
if (all(yv == 0))
  stop("ares: family = 'poisson' requires at least one y > 0;",
       " all observed counts are zero.")
```

Same pattern can be applied to `family = "gamma"` if y is constant.
