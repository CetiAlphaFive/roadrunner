# KRLS classification + non-gaussian response — research + plan

**Status**: deferred. LPM (linear-probability-model treatment of 0/1 y)
is acceptable in the meantime — that is the current behaviour, since
[`krls()`](https://cetialphafive.github.io/roadrunner/reference/krls.md)
does not reject 0/1 y.

## The problem

KRLS is fundamentally an **L2 ridge regression in an RKHS**:

``` math
\hat{c} = (K + \lambda I)^{-1} y \quad\Leftrightarrow\quad \min_c \|y - Kc\|^2 + \lambda c^\top K c
```

For binary `y \in \{0, 1\}`, naive KRLS produces a **linear-probability-
style fit** — predictions can fall outside `[0, 1]`, marginal effects
are on the probability scale (fine), but the noise model is gaussian
(wrong) so SEs / CIs are mis-calibrated, and there is no
`predict(..., type = "link")` notion.

For poisson / gamma responses, KRLS gives a mean-fit that is not
constrained to be positive, no MLE-correct dispersion estimate, and no
log-scale predictions.

## Literature

Three established families of fix:

### A. Linear probability KRLS (do nothing different)

- Hainmueller & Hazlett 2014 explicitly note that KRLS-on-binary-y gives
  LPM-style fits, suitable for marginal-effect *interpretation* but not
  probability *prediction*.
- This is what bigKRLS and KRLS::krls effectively do — they just do not
  reject 0/1 y.
- **Pros**: trivial; byte-equal to current implementation if user passes
  0/1 y.
- **Cons**: predictions outside \[0,1\]; no link function; no `pint` PIs
  that respect the bounded scale.

### B. Kernel logistic / GLM regression via IRLS

(Zhu & Hastie 2005, “Kernel Logistic Regression and the Import Vector
Machine”; Wahba et al. 1995.)

- Replace squared loss with deviance. Solve
  ``` math
  \hat{c} = \arg\min_c \ell(y, Kc) + \lambda c^\top K c
  ```
  by iteratively reweighted least squares: each Newton step is a
  weighted ridge in the original kernel space.
- Each IRLS iter: given `mu = link^{-1}(Kc)`,
  `W = diag(d mu/d eta)^2 / Var(y|mu)`,
  `z = Kc + (y - mu)/(d mu/d eta)`, solve
  `(W^{1/2} K W^{1/2} + lambda I) (W^{-1/2} c_new) = W^{1/2} z`.
- 5–15 iterations to converge typically.
- **Pros**: principled MLE in RKHS; correct deviance; correct SEs via
  the Fisher information at convergence; link-respecting predictions.
- **Cons**: weights break the cached eigen structure each iter — each
  iter is O(n^3) (re-eigendecompose `W^{1/2} K W^{1/2}`) OR use a more
  clever update. ~5x cost of gaussian fit.

### C. Post-hoc GLM on eigen features

(Mirrors ares’ family pattern; finite-rank kernel-GLM.)

- Eigendecompose `K = V D V'` (already done for gaussian fit).
- Truncate to `k` largest eigenvalues that explain \>= 99% variance
  (typically `k \approx 50-200` for `n = 1000+`).
- Treat the truncated eigenvectors `V_k` (n x k) as a fixed feature
  matrix.
- Run `stats::glm.fit(V_k, y, family = ...)` with `k << n` parameters
  -\> well-posed, no extra ridge needed inside the GLM (eigtrunc IS the
  ridge).
- Predictions on new data `X_new`: build `K_new` (m x n), project the m
  x k eigen-feature matrix at `X_new`, multiply by `beta_glm` -\> linear
  predictor `eta`. Apply `link^{-1}` -\> response-scale predictions.
- **Pros**: reuses the gaussian-KRLS infrastructure entirely; only one
  extra call to `glm.fit`; fast (k x n vs n x n linear system); aligned
  with ares’ “fit basis as gaussian, post-hoc GLM” pattern.
- **Cons**: not the exact MLE of `min deviance + lambda c'Kc`; it is the
  MLE *given the truncated eigen-basis*. For `k` chosen sensibly this is
  indistinguishable from the full KLR solution at the precision users
  actually care about (~1e-3 on predictions).

## Recommendation

**Go with C (post-hoc GLM on eigen features)**, with B as a stretch
option for users who want strict KLR.

Rationale:

1.  **Mirrors ares pattern exactly.** Syntax-alignment with ares was the
    user’s stated goal for the last alignment chain. ares does post-hoc-
    GLM-on-selected-basis. Same here.
2.  **Cheap.** Adds one
    [`stats::glm.fit`](https://rdrr.io/r/stats/glm.html) call per fit.
    No new C++. No re-eigendecomp per IRLS iter.
3.  **No new theory to validate.** Eigen-feature GLM is a standard
    reduced-rank GP-classification approximation; literature blesses it.
4.  **Marginal effects clean.** Chain rule on existing `krls_deriv_cpp`:
    `d response / d x_k = d link^{-1}/d eta * d eta / d x_k`. Single
    multiply.
5.  **B is a future option.** If a user demands true KLR (e.g. for a
    small-n, high-imbalance binary problem where the eigen-truncation is
    too crude), add `family.method = c("eigen-glm", "irls")`. Until
    then, simpler is better.

------------------------------------------------------------------------

# Phased implementation plan

Mirror ares signature. Default `family = "gaussian"` (current
behaviour). New families opt-in.

## Phase A — binomial (HIGHEST VALUE)

**Scope**:

- Add `family = c("gaussian", "binomial", "poisson", "gamma")` arg to
  `krls.default` (port validation block from `ares.default` lines
  286-360).
- Add `family.eigtrunc.var = 0.99` arg (eigen-feature retention
  threshold; cap at k \<= n-1).
- For binomial: validate `y \in {0, 1}` or 2-level factor / character /
  logical; coerce to 0/1 numeric; stash factor levels for predict
  roundtrip (same as ares).
- Fit pipeline:
  1.  Run existing gaussian KRLS path on **centered y** (`y - mean(y)`)
      to get coefs / eigendecomp / sigma / lambda.
  2.  After the gaussian fit converges, compute `K %*% c` once -\> that
      is the gaussian-scale linear predictor.
  3.  **Refit step**: take `V_k` (eigvectors with cumulative variance
      \>= 0.99), run `stats::glm.fit(V_k, y, family = binomial())`.
      Store `beta_glm` (length k).
  4.  Replace `out$coeffs` with the **kernel-coefficient equivalent**
      `c_glm = V_k beta_glm` (length n). This makes
      `K %*% c_glm = V_k beta_glm = eta_hat`, the linear predictor.
      Stash the original gaussian solve under `out$coeffs_gauss`.
  5.  Stash `out$family = "binomial"`, `out$y_levels`,
      `out$family_fit = list(beta = beta_glm, k_used = k, V_k = V_k)`.
  6.  `fitted` is the **response-scale** prediction: `plogis(eta_hat)`.
      For non-gaussian families, skip the y-standardisation entirely.
- `predict.krls_rr`:
  - Add `type = c("response", "link")` arg (default `"response"`).
  - Build `K_new`, compute `eta_new = K_new %*% c_glm`.
  - `type = "link"`: return `eta_new`.
  - `type = "response"`: return `plogis(eta_new)`.
  - `se.fit`: derived from `vcov(glm.fit)` projected through
    `K_new V_k`. Standard GLM delta method.
  - For 2-level factor responses, optional: when caller passes a factor
    as `y`, predict returns a factor when `type = "response"` is rounded
    — but default behaviour is to return probabilities (matches
    `predict.glm(type="response")`).
- **Marginal effects**: chain rule.
  - Currently `krls_deriv_cpp` returns `d eta / d x_k` (gaussian KRLS
    deriv). For binomial: `d p / d x_k = p (1 - p) . d eta / d x_k`.
  - `binary = TRUE` first-differences: replace per-pair predictions with
    `p(x_high) - p(x_low)` on the response scale.
- **Tests**:
  - Recovery: train on logistic DGP
    `y = 1{plogis(0.5*x1 + 0.5*x2) > U(0,1)}`, n=500. Check AUC \>= 0.85
    on holdout; calibration plot lines up.
  - vs `glm(y ~ poly(x, 2), family = binomial)` on a known nonlinear DGP
    — krls should match or beat AUC.
  - vs ares `family = "binomial"` on same DGP — krls is allowed to
    differ but ROC should be in the same ballpark (within 0.02 AUC).
  - Round-trip 2-level factor y -\> factor predictions match input level
    encoding.
  - `predict(fit, type = "link")` returns gaussian-scale eta;
    `type = "response"` returns probabilities in \[0, 1\].

**Scope estimate**: ~200 LOC R, ~150 LOC tests. No new C++.

## Phase B — poisson (log link)

**Scope**:

- Validation block (port from `ares.default` lines 328-349):
  non-negative, integer-valued (with tolerance).
- Fit pipeline same as Phase A but
  `glm.fit(V_k, y, family = poisson())`.
- `out$fitted = exp(eta_hat)`.
- `predict(type = "response") = exp(K_new %*% c_glm)`,
  `type = "link" = eta`.
- Marginal effects: `d lambda / d x_k = lambda . d eta / d x_k`.
- `binary` first-differences: `exp(eta(x_high)) - exp(eta(x_low))` on
  rate scale.

**Tests**:

- Recovery on counts DGP `y ~ Poisson(exp(0.5*x1 - 0.3*x2))`, n=500.
  RMSE on rate scale better than
  `glm(y ~ poly(x, 2), family = poisson)`.
- vs ares poisson on same DGP.

**Scope estimate**: ~80 LOC R, ~80 LOC tests.

## Phase C — gamma (log link)

Same as Phase B with `family = Gamma(link = "log")`. Validation: y \> 0
strict.

**Scope estimate**: ~30 LOC R, ~50 LOC tests.

## Phase D — family-specific prediction intervals

KRLS’ current `interval = "pint"` errors for non-gaussian. Add:

- **binomial**: Wald or Wilson on `eta`, then squash:
  - `eta_lo, eta_hi = eta +/- z * se(eta)`
  - `lo, hi = plogis(eta_lo), plogis(eta_hi)`
  - Wilson interval is the safer choice for finite-sample binomial
    coverage near 0/1.
- **poisson**: Wald on log scale:
  - `lo, hi = exp(eta +/- z * se(eta))`
- **gamma**: same as poisson on log scale.

`varmod` for non-gaussian: probably skip — the GLM already has
dispersion baked in. Document that `varmod` only applies to gaussian
fits.

**Tests**: empirical coverage at nominal 95% in \[0.92, 0.98\] on a
clean DGP per family.

**Scope estimate**: ~60 LOC R, ~80 LOC tests.

## Phase E — marginal effects on response scale (chain rule)

Currently `derivatives` and `avgderivatives` are computed on the eta
scale. For non-gaussian, users probably want them on the response scale
(interpretable as “probability per unit x” for binomial).

- Add `marginal.scale = c("response", "link")` arg, default `"response"`
  for non-gaussian, `"link"` ignored for gaussian.
- Multiply existing `derivmat_s` by the link-derivative evaluated at
  each row’s fitted `mu`.
- For binomial: factor of `mu (1 - mu)`.
- For poisson / gamma log link: factor of `mu`.
- `var.avgderivatives` gets the same chain rule applied to the analytic
  variance derivation.

**Tests**: at known DGP, `avgderivatives` on response scale should be
close to numerical-derivative-of-fitted; on link scale they are the raw
gaussian-KRLS derivatives.

**Scope estimate**: ~80 LOC R, ~80 LOC tests.

## Phase F — multinomial (PUNT unless demand)

If a user needs multi-class:

- Add `family = "multinomial"` with
  [`nnet::multinom`](https://rdrr.io/pkg/nnet/man/multinom.html) as a
  `Suggests` dep.
- Same eigen-feature pattern: `nnet::multinom(y ~ V_k - 1)` for k-class.
- Predictions: softmax over K classes.
- Marginal effects: K x p matrix per observation.

**Scope**: ~150 LOC R, ~120 LOC tests, adds nnet dep.

**Recommendation**: do not ship in this round. Add only if a user files
an issue.

------------------------------------------------------------------------

## Sequencing

`A -> B -> C -> D -> E -> F (deferred)`

Phase A is the value driver (binomial is what users actually ask for
first). Phases B, C, D, E ride on the eigen-feature infrastructure A
builds and are mostly mechanical.

Total scope if A-E ship: ~450 LOC R, ~440 LOC tests, 0 new C++, 0 new
deps.

## Open questions to resolve before kickoff

1.  **eigen-feature variance threshold default**: 0.99 is conservative
    (keeps many eigvecs). For binary classification of n=1000 with a
    clear signal, k=20 might be plenty. Could expose `family.k` directly
    OR `family.eigtrunc.var = 0.99` (relative threshold) — pick the
    latter since it auto-adapts to spectrum.

2.  **What happens to `out$coeffs` for non-gaussian fits?**

    - **Option 1 (preferred)**: `out$coeffs` is always the kernel-
      coefficient vector that gives `K %*% out$coeffs = eta_hat`. Add
      `out$coeffs_gauss` for the original gaussian solve. This way
      `predict.krls_rr` internally works the same for all families —
      just multiply by `K_new`, then apply link.
    - **Option 2**: `out$coeffs` is undefined for non-gaussian; only
      `out$family_fit$beta` is meaningful. Worse — breaks the
      `K %*% coeffs` interpretation users may rely on.

3.  **Should `lambda.method = "cv"` (Phase 3 from last alignment)
    interact with non-gaussian fits?** The gaussian eigen-fit picks
    lambda by LOO MSE; for binomial we might want lambda chosen by
    CV-deviance. Likely yes — extend Phase 3’s CV mechanism to score by
    deviance when family is non-gaussian. Small change, but worth
    flagging.

4.  **`vcov` for non-gaussian**: GLM produces its own vcov on
    `beta_glm`. Project back via `V_k`:
    `vcov(c_glm) = V_k vcov(beta_glm) V_k^T`. Document that this is the
    asymptotic GLM covariance, not the gaussian-ridge sandwich.
