# BUG-002 — predict(bagged_fit, newdata=NULL) silently returns CENTRAL-only

**Severity:** correctness (medium).
**File:** `R/predict.R:46-67`.

## Symptom

For any fit with `n.boot > 0`:

```r
predict(fit)            # newdata=NULL  -> returns object$fitted.values
                        # (the CENTRAL fit only, NOT the bag mean)
predict(fit, x_training) # newdata=x   -> returns the BAG MEAN
```

These two calls give different answers, diverging by up to ~0.34 on a
5-fold bagged gaussian fit (probe 06B) and ~0.35 in probability on a
binomial bag (probe 06).

## Root cause

```r
predict.ares <- function(object, newdata = NULL, ...) {
  ...
  if (is.null(newdata)) {
    yhat <- if (fam == "binomial" && type == "link") {
      object$linear.predictor %||% object$fitted.values   # <- central only
    } else {
      object$fitted.values                                  # <- central only
    }
    if (isTRUE(se.fit) && !is.null(object$boot)) {
      attr(yhat, "sd") <- .ares_boot_sd(object, NULL)       # <- returns NA
    }
    ...
    return(yhat)
  }
```

The bag-averaging logic only runs when `newdata` is non-NULL (predict.R:
192-251). With `newdata = NULL` we shortcut to `object$fitted.values`
which is the **central** fit (i.e. the post-process at ares.R:643-644:
`out$fitted.values <- drop(out$bx %*% out$coefficients)`, computed once
before bagging fires).

Same issue with `se.fit = TRUE`: `.ares_boot_sd(object, NULL)` returns a
length-n vector of `NA_real_` (see predict.R:276-288 — the helper
explicitly bails out because it can't recover training x from
`object$bx`).

## Repro

```r
library(roadrunner)
set.seed(70)
n <- 200; p <- 4
x <- matrix(runif(n * p), n, p)
y <- 3 * x[, 1] + 2 * x[, 2]^2 + rnorm(n, sd = 0.4)
fb <- ares(x, y, n.boot = 5, seed.cv = 5L, nthreads = 2)

pn   <- predict(fb)     # central only
pxin <- predict(fb, x)  # bag mean
max(abs(pn - pxin))     # 0.342  -- silent disagreement
```

## Impact

The docstring at predict.R:19-21 says:

> For bagged fits (n.boot > 0), returns the mean across the central fit
> and the bootstrap replicates; with se.fit = TRUE, the per-row bag
> standard deviation is attached as an "sd" attribute.

This contract is violated in the `newdata = NULL` path. Users running
`predict(fit)` (the most natural form) get the WRONG kind of prediction
without warning.

## Suggested fix

In the `newdata = NULL` branch, when `!is.null(object$boot)`, fall
through to the bag-averaging logic using training x. The training x is
not stored on the fit, so two options:

1. Store `object$x` at fit time (memory cost: 8 bytes × n × p).
2. Reconstruct an effective basis from `object$bx` for the central path,
   but each bag replicate's basis can't be reconstructed without x. So
   option 1 is the clean fix.
3. Document loudly that `predict(fit)` for bagged returns the central
   only, and `predict(fit, x_train)` for the bag mean. (Hack workaround,
   not a true fix.)

The implementation comment at predict.R:276-289 already notes the
problem: "Reconstruct training x from object$bx is impossible — bag-on
train SD requires the original x. We don't keep it." Solving this
properly requires storing x, which is the right call.
