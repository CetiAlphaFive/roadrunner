# roadrunner v0.0.0.9028 — adversarial audit (2026-05-11)

Scope: `ares()` MARS fitter (R/ares.R, R/predict.R, src/ares.cpp).
Approach: spec-driven probes against Friedman (1991) + earth conventions
+ implementation-internal contract from CLAUDE.md. Read-only audit
(no production code modified).

Test baseline: ARES_FULL_TESTS=1 -> 276/276 PASS. The bugs below are
*not* caught by the existing test suite.

## Confirmed findings, ranked (correctness > robustness > docs)

### CORRECTNESS — 3 findings

#### BUG-001 (CORRECTNESS, high)
**Unseeded binomial/poisson/gamma bag refits use bootstrap rows that
DO NOT match the basis-selection bootstrap.**

- File: `R/ares.R:840-877` (.ares_refit_boot_binomial),
  `R/ares.R:940-982` (.ares_refit_boot_glm_loglink).
- Trigger: `family ∈ {binomial, poisson, gamma}` + `n.boot > 0` +
  `seed.cv = NULL`.
- Effect: The post-hoc GLM refit for each bag replicate is fitted on
  bootstrap rows drawn from the CURRENT R RNG state, not from the
  identical bootstrap rows used in the gaussian-pass main bagging loop
  (which also drew from the live RNG, but has since advanced through
  3 sample.int() calls + intervening C++ work). Result: the basis was
  selected on one bootstrap sample; the GLM coefficients are fitted on
  a different one. Bag predictions are still produced, but the
  statistical interpretation as "bootstrap of the full pipeline" is
  broken.
- Repro: `/tmp/audit/probe_08_bag_binomial_unseeded.R`,
  `/tmp/audit/probe_14b_unseeded_drift.R`.
- Spec citation: CLAUDE.md §"Bagging composes" implies each bag
  replicate is a coherent bootstrap; the bag refit comment at
  ares.R:835-838 says "we rebuild the bag's basis on a fresh
  bootstrap sample (using the same seed offset as the original bagging
  loop)" — but this is only true when seed_cv is non-NULL.
- Fix sketch: always set.seed(_some derived value_) in both loops, or
  carry the indices forward from the main loop.

#### BUG-002 (CORRECTNESS, medium)
**`predict(bagged_fit)` with newdata=NULL returns the CENTRAL fit's
fitted values, not the bag mean.**

- File: `R/predict.R:46-67`.
- Trigger: any fit with `n.boot > 0`, calling
  `predict(fit)` (no newdata) vs `predict(fit, x_training)`.
- Effect: silently inconsistent predictions. Demonstrated divergence
  of up to 0.34 on a 5-fold bagged gaussian fit (probe 06B) and 0.35
  on a binomial classifier (probe 06).
- Same path: `se.fit = TRUE` with newdata=NULL returns
  `attr(yhat, "sd") = NA` (see `.ares_boot_sd` predict.R:276-289),
  but `se.fit = TRUE` with newdata=x_training returns the actual bag
  SD.
- Repro: `/tmp/audit/probe_06_predict_consistency.R`,
  `/tmp/audit/probe_06b_predict_bag_null.R`.
- Spec citation: docstring at predict.R:19-21 says "for bagged fits
  ... returns the mean across the central fit and the bootstrap
  replicates" — this is violated when newdata is omitted.
- Fix sketch: in the newdata=NULL branch (predict.R:47-67), if
  `!is.null(object$boot)`, fall through to the bag-averaging path
  using `object$bx` as the design (it's stored on the central fit) for
  the central, and rebuilding basis on the training x via stored
  `object$bx`-equivalent. The current code shortcuts to
  `object$fitted.values` which is wrong for bagged fits.

#### BUG-003 (CORRECTNESS, medium)
**`varmod = "lm"` PI silently collapses to ~zero width when the
predicted MAD `(intercept + slope * yhat)` goes negative.**

- File: `R/predict.R:265-269`.
- Trigger: `varmod = "lm"` + a DGP where residual magnitude DECREASES
  with yhat + prediction at `yhat` values where the linear MAD model
  hits zero or below. Extrapolation past training yhat range makes
  this routine.
- Effect: the floor `pmax(..., 1e-12)` caps sigma at ~1e-12, producing
  PI widths of ~4e-12. The user gets a "point estimate disguised as a
  PI" with no warning.
- Repro: `/tmp/audit/probe_05c_pi_negative_sigma.R` shows 8 of 10
  extrapolated points with PI width ≈ 4e-12.
- Spec citation: the docstring at ares.R:131-134 says varmod="lm"
  "fits a small linear model of |resid| on yhat to allow simple
  mean-dependent heteroscedasticity" — no warning that this can
  produce degenerate PIs.
- Fix sketch: warn or stop when `intercept + slope * yhat` is non
  positive anywhere in the prediction set; or floor at
  `max(0.01 * fitted_sigma, min(non-neg MAD predictions))`; or
  clip to the in-sample yhat range and warn.

### ROBUSTNESS — 4 findings

#### BUG-004 (ROBUSTNESS, medium)
**OOV factor level in newdata silently DROPS rows from the prediction
output without warning.**

- File: `R/predict.R:99-117`.
- Trigger: data-frame newdata containing a factor/character column
  with a level not present at training time.
- Effect: out-of-vocabulary levels become NA after `factor(..., levels
  = fi$xlevels[[jname]])`. Then `model.matrix(~ ., newdata)` calls
  `na.action = na.omit` (the global default) and silently drops the
  affected rows. Result: `predict()` returns FEWER rows than `nrow(
  newdata)`. Probe 09B: passing 50 rows with 5 OOV returned only 45
  predictions.
- Repro: `/tmp/audit/probe_09b_oov.R`.
- Fix sketch: detect OOV explicitly before model.matrix, either
  (a) error with a clear message, (b) emit NA-prediction rows in the
  output and warn, or (c) impute to the most frequent training level
  with a warning.

#### BUG-005 (ROBUSTNESS, medium)
**Zero-weight observations bias GCV downward, causing extra terms to
be retained vs the equivalent dropped-rows fit.**

- File: `src/ares.cpp:1206-1211` (compute_gcv lambda),
  `R/ares.R:1008` (varmod df), `src/ares.cpp:1208` (denom uses n).
- Trigger: `weights` containing exact zeros (or values numerically
  indistinguishable from zero after sqrt) with other weights > 0.
- Effect: GCV uses raw row count `n` instead of effective sample size
  `n_eff = sum(w_i > 0)`. Probe 12: 50 of 200 rows zeroed -> 10 terms
  retained vs 7 in the dropped-rows fit; sigma_hat 0.4417 vs 0.4300
  (with df 93 vs 81). The penalty per knot is implicitly under
  estimated, the model is over-fit, and PI widths use the wrong df.
- Repro: `/tmp/audit/probe_11_zero_weights.R`,
  `/tmp/audit/probe_12_gcv_weights.R`.
- Spec citation: R-side check at ares.R:339-341 only rejects
  `mean(weights) == 0`, not zero-weight rows individually.
- Fix sketch: either reject zero weights (require all weights > 0),
  or pass `n_eff = sum(w > 0)` to the C++ engine and use it in GCV /
  auto_minspan / auto_endspan / varmod df.

#### BUG-006 (ROBUSTNESS, low-medium)
**`varmod = "lm"` only captures heteroscedasticity whose residual
scale depends on yhat — not on x directly.**

- File: `R/ares.R:1020-1046`.
- Trigger: DGP where `sigma(x)` varies with a predictor that
  contributes little to `yhat`.
- Effect: the |resid|~yhat slope ends up near zero, the lm-varmod
  collapses to const-varmod, and PI coverage degrades catastrophically
  in high-variance regions. Probe 05B: nominal 95% PI achieved 79%
  coverage at x1 ∈ (0.8, 1] when sigma = 1+4*x1 (so sigma_true = 4.6
  in that bin, but lm-varmod predicted sigma ≈ 3.25 everywhere).
- Repro: `/tmp/audit/probe_05b_pi_diagnostic.R`.
- Spec citation: docstring at ares.R:131-133 only says "simple
  mean-dependent heteroscedasticity"; the audit shows that calling
  this even "simple heteroscedasticity" overpromises — it captures
  ONLY yhat-dependent.
- Fix sketch: either (a) add a `varmod = "x"` flavour that regresses
  |resid| ~ x via a small MARS fit (like earth's
  `varmod.method="x"`), or (b) tighten the docstring to say
  "yhat-dependent" and explicitly warn that x-driven
  heteroscedasticity is not captured.

#### BUG-007 (ROBUSTNESS, low)
**`family = "poisson"` accepts all-zero `y`, producing a degenerate
GLM fit silently.**

- File: `R/ares.R:282-296`.
- Trigger: `family = "poisson"` with `y = rep(0L, n)`.
- Effect: the y-validation block only checks `y >= 0` and integer
  valued, so `y = 0` everywhere passes. `glm.fit` returns deviance=0,
  null.deviance=0, coefficients with NaN/zero handling — the fit
  technically "succeeds" but is mathematically degenerate.
- Repro: `/tmp/audit/probe_09_misc_edges.R`.
- Fix sketch: also reject `all(y == 0)` for poisson with a clear
  error.

### DOCS — 2 findings

#### DOC-001
**`ares()` accepts character `y` for binomial but docstring at
ares.R:39-40 only mentions "0/1 numeric, logical, or 2-level factor".**

- File: `R/ares.R:39-40, 260-267`.
- Effect: feature surface is wider than documented. Not a bug; worth
  syncing the docstring.

#### DOC-002
**Bagged + varmod="const" + interval="pint" produces PIs whose width
is CONSTANT across rows** (uses central fit's sigma_hat, doesn't
combine with bag SD).

- File: `R/predict.R:237-250`.
- Spec citation: the comment at predict.R:245-249 documents this
  decision ("PIs use the bag mean for `fit`; bag SD is approximate
  for the variance model component but we don't combine them
  here"). Defensible choice but a user expecting bag-aware PI widths
  may be surprised. Worth promoting that comment to the docstring.

## Verified PASSES (no bug)

These spec-checks succeeded; the engine is correct on these axes:

1. **GCV formula** (probe 01): engine GCV matches Friedman/earth
   formula to 1e-15.
2. **Determinism contract** (probe 03): nthreads=1 == nthreads=N
   byte-identical across {default, degree=2, weights, varmod,
   binomial, n.boot, autotune, cv, poisson, gamma}. 10/10 PASS.
3. **WLS coefficient equivalence** (probes 02, 02C): WLS engine matches
   `lm.wfit` to 5e-15 on full-rank basis; matches OLS on
   row-replicated data when integer weights.
4. **Constant weight scaling invariance** (probe 02A): scaling all
   weights by a constant leaves coefficients, dirs, cuts, selected.terms
   unchanged byte-for-byte.
5. **Post-hoc GLM refit** (probe 04): engine coefs = `stats::glm.fit`
   coefs to 0.000e+00 for binomial, poisson, gamma — including weighted.
6. **predict link/response inversion** (probe 04, 06): max |response
   - plogis(link)| = 0 for binomial; max |response - exp(link)| = 0
   for poisson (single fit).
7. **CV pruning reproducibility** (probe 07): identical fits under
   fixed seed.cv; global RNG NOT advanced when seed.cv is set;
   1-SE size ≤ argmin size as specified.
8. **Backward OLS validity** (probe 10): engine coefs match `lm.fit`
   at the selected basis to 1.776e-15; backward isn't globally optimal
   (expected: it's greedy).
9. **Autotune reproducibility** (probe 13): byte-identical under
   seed.cv; warmstart-fired vs full-grid produce different but
   reproducible winners.
10. **PI empirical coverage on homoscedastic DGP, const varmod**
    (probe 05A): mean=0.945 (nominal 0.95) over 20 seeds. ✓
11. **GCV is correct under POSITIVE weights** with mean(w)=1
    normalisation (probe 12).

## Bottom line

Engine math is solid: GCV, WLS, QR-downdate, GLM refit, determinism.
The bugs are all in **the surface area around the C++ engine** —
the R-side bagging-refit RNG (BUG-001), the predict() consistency
(BUG-002), the varmod design (BUG-003, BUG-006), the OOV / zero-weight
edge cases (BUG-004, BUG-005), and the poisson y-validation
(BUG-007).

Severity ordering for triage:
  - BUG-001 first (statistical correctness of bagged classifiers).
  - BUG-002 next (silent prediction inconsistency).
  - BUG-003 next (degenerate PIs without warning).
  - BUG-004 / BUG-005 next (silent row drops / over-fitting).
  - BUG-006 / BUG-007 last (limitation + degenerate-but-flagged path).
