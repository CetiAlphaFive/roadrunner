# BUG-001 — Unseeded binomial/poisson/gamma bag refits use mismatched indices

**Severity:** correctness (high).
**Files:** `R/ares.R:840-877` (.ares_refit_boot_binomial), `R/ares.R:940-982`
(.ares_refit_boot_glm_loglink).

## Symptom

When `family ∈ {binomial, poisson, gamma}` and `n.boot > 0` and `seed.cv =
NULL` (the default), the post-hoc GLM refit for each bag replicate is
fitted on bootstrap rows that DO NOT match the rows the basis was selected
on in the main bagging loop. The two loops both consume RNG from
`.Random.seed`, but the refit loop runs AFTER the main loop, so its
sample.int() calls produce DIFFERENT bootstrap indices.

## Root cause

```r
.ares_refit_boot_binomial <- function(fits, x_full, y_full, seed_cv,
                                      weights = NULL) {
  n_full <- nrow(x_full)
  has_seed <- !is.null(seed_cv)
  if (has_seed) {                                # <-- ONLY here is set.seed called
    ...
    set.seed(as.integer(seed_cv) + 1009L)
  }
  for (b in seq_along(fits)) {
    idx <- sample.int(n_full, n_full, replace = TRUE)   # <-- consumes RNG either way
    ...
```

In the main bagging loop in `ares.default` (ares.R:654-707), the same
guard exists at lines 665-677: `set.seed()` is only called when seed.cv is
non-NULL. So under unseeded usage:

1. Main bagging loop calls `sample.int(n, n, replace = TRUE)` three times,
   advancing the user's `.Random.seed` past those three draws.
2. The refit loop then calls `sample.int(n, n, replace = TRUE)` three more
   times, continuing from the NEW `.Random.seed` — so it gets indices
   4–6, not indices 1–3.

The basis stored in `fits[[b]]$dirs / cuts / selected.terms` was chosen on
sample 1 (gaussian forward+backward), but the GLM coefficients are then
estimated on sample 4. The two bootstrap samples have *no shared rows on
expectation* beyond what random overlap brings.

## Repro

```r
library(roadrunner)
set.seed(50)
n <- 200; p <- 4
x <- matrix(runif(n*p), n, p)
yb <- as.integer(plogis(2*x[, 1] - 0.5) > 0.5)

# Replicate what ares does in the unseeded path:
set.seed(50)
idx_main <- replicate(3, sample.int(200, 200, replace = TRUE), simplify = FALSE)
# Engine-side mars_fit_cpp consumes NO R RNG. So the refit loop picks up
# wherever the main loop left off.
idx_refit <- replicate(3, sample.int(200, 200, replace = TRUE), simplify = FALSE)
identical(idx_main[[1]], idx_refit[[1]])   # FALSE — bug.
```

In the seeded path (`seed.cv = 1L`), both loops call `set.seed(1L + 1009L)`
right before their sample.int sequence, so they get identical indices.
Only the unseeded path is broken.

## Impact

Statistical interpretation is broken. The bag is no longer a coherent
bootstrap of the full pipeline. Predictions are still produced and SD
attributes still work — but a user reasoning about "what's my bag SD
estimating" gets the wrong answer.

## Suggested fix

Either:
1. **Carry indices forward.** Save `idx_list` in the main bagging loop
   and pass it to the refit function instead of redrawing.
2. **Always seed.** In both loops, set up a sub-RNG via `set.seed` with a
   derived seed (e.g., `set.seed(digest::digest2int(...))` or any
   deterministic offset of a captured global state); the on.exit restore
   pattern is already there.
3. **Reject unseeded bagged GLM fits.** Stop with an instructive error
   message if `n.boot > 0 && family != "gaussian" && is.null(seed.cv)`.

Option (1) is the cleanest — eliminates the entire RNG coupling between
the two loops.
