# roadrunner - handoff brief (post-rename roadmap)

## What this package is

`roadrunner` is an R package collecting fast, low-dependency machine
learning algorithm implementations. The first algorithm shipped is
`ares()`, a Multivariate Adaptive Regression Splines (MARS) fitter built
with `Rcpp` + `RcppParallel` (TBB). Cross-platform. Deterministic across
threads (byte-identical fits at `nthreads=1` vs N). Additional algorithms
will plug in as separate exported functions under the same package roof.

The package was renamed from `ares` to `roadrunner` at v0.0.0.9028. The
MARS fitter remains exposed as `ares()`, the S3 class is still `"ares"`,
and `predict.ares` / `print.ares` / `summary.ares` / `plot.ares` are
unchanged. The working directory stays `/home/jack/Dropbox/ares`.

**Author**: Jack Trametta <jtrametta@gmail.com>. **License**: MIT. **Repo**:
local git only on `main`.

`earth` (`earth::earth`) is a useful diagnostic but not a parity target.

## Current state - v0.0.0.9028

### Test suite
- 276/276 testthat green under ARES_FULL_TESTS=1 (202 prior + 51 weights /
  poisson / gamma + 23 from v0.0.0.9026-v0.0.0.9027). Determinism contract
  holds across every new code path: `nthreads=1` and `nthreads=N` produce
  byte-identical fits at fixed seed, and the unweighted (default) path is
  byte-identical to pre-weights v0.0.0.9024.
- `R CMD check`: 0 errors, 0 warnings, 1 note (system non-portable
  compile flags, pre-existing).

### v0.0.0.9025-v0.0.0.9027 earth-parity roadmap (this session)

**v0.0.0.9025 — observation weights + poisson / gamma families (Phase A)**:
- `weights` arg on ares(). Engine-level WLS via a sqrt(w) transform of
  the design and response: B_internal[:, t] = sqw .* B_orig[:, t],
  Y_internal = sqw .* y. Forward + backward both consume weights
  (KnotScanner prefix sums and the QR/downdate work in transformed
  space). GCV denom uses normalised w (mean=1, so n_eff=n). Output bx
  is divided by sqw at the C++ boundary so predict() sees the
  original-scale basis. `weights = rep(1, n)` is byte-identical to
  the unweighted path. CV pruning and autotune subset weights per fold
  and score test-side weighted MSE. Bagging resamples weights with the
  bootstrap and renormalises per replicate.
- `family = "poisson"` (log link) and `family = "gamma"` (log link).
  Mirror of the binomial pattern: forward + backward run as gaussian;
  the selected basis is refit via stats::glm.fit() under the
  appropriate GLM. predict() honours type = {"response", "link"} for
  log-link families. Weights pass through to the post-hoc GLM.

**v0.0.0.9026 — binomial default flips (Phase B)**:
- `auto.linpreds` defaults to TRUE and `adjust.endspan` defaults to 2L
  **only for family = "binomial"**. Gaussian / poisson / gamma keep
  the conservative defaults (auto.linpreds = FALSE, adjust.endspan = 1L)
  because the v0.20-mlbench grid showed regressions on regression DGPs.
  Decision was driven by inst/sims/v0.25_ab_binomial.R: A/B test on
  twonorm / threenorm / ringnorm / waveform at n = {500, 1500} x 6 reps.
  Flipped defaults tied or improved mean AUC on all 8 (dgp, n) cells.
- Result CSV: inst/sims/results/v0.25-ab-binomial.csv.

**v0.0.0.9027 — prediction intervals (Phase D)**:
- `varmod = c("none", "const", "lm")` arg on ares(); default "none"
  so SEs are off by default (explicit user requirement).
  - "const": store sigma_hat from (weighted) training RSS over the
    residual df. PI = yhat +/- qt(1 - alpha/2, df) * sigma_hat.
  - "lm": fit |resid| ~ yhat closed-form (weighted closed form when
    weights are supplied) and apply sqrt(pi/2) to convert
    mean-absolute-deviation back to a sigma estimate. Captures simple
    heteroscedasticity.
- `predict.ares(interval = "pint", level = 0.95)` returns a matrix
  with columns c("fit", "lwr", "upr"). Errors loudly if the fit was
  built without a varmod, or if family != "gaussian". Heavy test in
  tests/testthat/test-varmod.R verifies empirical 95%-PI coverage on
  a clean gaussian DGP at n = 1000 lies in [0.90, 0.99].

### Punted in this session

**Item 4 (multi-response gaussian) — punted to roadmap**.
Rationale: the KnotScanner's fast-LS prefix sums (T_pr, S_pr, T_pxr, S_pxr,
T_Qp, T_Qpx etc.) and AVX2 inner loops are vectorised over Mq (basis
dimension), and the residual is a single n-vector. A K-response extension
needs either (a) outer loop over responses (loses AVX over Mq, hits the
parallel-over-pairs path K times) or (b) restructuring inner loops to
vectorise over K. The scope estimate exceeded the "~300 lines new code"
ceiling in the roadmap once the predict / CV / autotune / bagging
compose-paths are accounted for. The user's prompt allows punting if too
big; the WLS engine surgery and four new feature surfaces (weights,
poisson, gamma, varmod) already shipped successfully in one session
without breaking the determinism invariant. Defer until there's user
demand.

### Inherited (still in place)

#### Pruning options
- `pmethod = "backward"` (default) - GCV-based subset selection.
- `pmethod = "none"` - keep all forward-pass terms.
- `pmethod = "cv"` - K-fold cross-validated subset selection. Promoted
  automatically when `nfold > 0`. Now weight-aware (CV folds use
  weighted MSE scoring for consistency with the WLS objective).

#### Autotune (`autotune = TRUE`)
- Inner-CV grid search over `(degree, penalty, nk, fast.k)`. Cells
  share a fold partition; successive halving + warm-start subsample
  short-circuit + shared-forward fast path. Weight-aware: train weights
  are subset per fold and renormalised to mean = 1 inside the C++
  engine; test-side scoring is weighted MSE.

#### Family
- "gaussian" (default), "binomial", "poisson", "gamma". All non-gaussian
  families use the earth strategy: forward + backward as gaussian on the
  selected basis, then post-hoc stats::glm.fit() with the appropriate
  family object.

### v0.20 mlbench-style benchmark (vs earth, ranger)
inst/sims/results/v0.20-mlbench.csv (baseline) and
inst/sims/results/v0.23-mlbench.csv (current). Highlights:
- ares default beats earth on 6/10 cells, ties most others.
- autotune wins on highdim cells by picking degree=2/3 correctly.
- ranger is ~3-4x worse on every regression cell - ares is the right
  MARS pick for continuous structured DGPs.
- Autotune wall-clock now: low-p sub-second (warm-start short-circuit
  + shared-forward), highdim p=20 ~5-15s (was 60-120s pre-v0.21).

### Determinism gotcha (still in place)
- The `-ffp-contract=off` flag was dropped at v0.0.0.9012 (FMA
  enabled). Run-to-run fits are still bit-identical.

## Speed status

- ares default vs earth: median ratio ~0.93x at 4t (ares slightly
  faster). 11/18 cells in inst/sims faster than earth at 4t.
- Autotune wall-clock: sub-second when warm-start fires; ~1-3s on
  small n full grid; ~5-15s on highdim p=20 full grid (post v0.23
  gates). Further large-p speed work belongs to deg=3 / `auto_linpreds`
  defaults — out of scope for current roadmap.
- Weighted fits have negligible overhead vs unweighted (one extra
  multiply per column added to B at C++ load).

## Scope discipline (locked-in)

- Response families: gaussian (default), binomial, poisson, gamma. All
  non-gaussian use post-hoc GLM refit on the selected basis.
- **Multi-response (y as matrix)**: NOT supported. See punt rationale
  above. The implementation cost exceeded the scope ceiling.
- `pmethod` in `{"backward", "none", "cv"}`. No exhaustive/forward/seqrep.
- No `fast.k` heuristic past the default (10) is intended; user can
  override via `autotune.speed`.
- Deps: hard = `Rcpp`, `RcppParallel`. Suggests = `earth`, `testthat`,
  `bench`. **No tidyverse, no rlang, no S4, no caret/tidymodels.**
- Variance model (`varmod`) is R-side only; no engine changes.
- No external tuning deps. License stays MIT.

## CRAN gotcha

Name collision with archived CRAN package `ARES`. Local install
unaffected. Before any CRAN submission, rename to `aresMARS`.

## Where things live

- `R/ares.R` - main function, formula + matrix interfaces, CV
  orchestrator (.ares_cv_fit), autotune orchestrator (.ares_autotune),
  binomial / poisson / gamma post-hoc GLM helpers
  (.ares_refit_binomial, .ares_refit_glm_loglink, bag-replicate
  variants), varmod helper (.ares_fit_varmod).
- `R/predict.R` - bagging-aware predict() with se.fit, type = response/
  link for binomial/poisson/gamma, and interval = "pint" PIs for
  gaussian + varmod fits.
- `R/print.R` - S3 methods.
- `src/ares.cpp` - C++ engine. Three exports, all weight-aware:
  `mars_fit_cpp` (forward + backward; new arg `weights_in`),
  `mars_basis_cpp` (basis builder for predict; weight-agnostic,
  always rebuilds untransformed basis),
  `mars_backward_only_cpp` (backward replay; new arg `weights_in`).
- `tests/testthat/` - 11 test files, 276 tests. New since v0.0.0.9024:
  test-weights.R, test-family-glm.R, test-varmod.R.
- `inst/sims/` - Monte Carlo + bench scripts. Latest results in
  `results/v0.23-mlbench.csv` (gaussian high-p benchmark),
  `results/v0.24-binomial.csv` (mlbench classification bench), and
  `results/v0.25-ab-binomial.csv` (binomial default-flip A/B test).
- `brain-contributions-v0.20.md` (root, Rbuildignored) — 10 distilled
  knowledge entries proposed for the statsclaw brain seedbank.
  Awaiting `/contribute` upload from a brain-connected session.

## Commit history (post-tuning roadmap)

```
93c6daa v0.0.0.9027: prediction intervals via residual variance model
e6beddc v0.0.0.9026: flip binomial defaults to earth-like
72557d1 v0.0.0.9025: weights arg + poisson/gamma families (Phase A)
5c4b22a v0.0.0.9024: family = "binomial" classification (Phase A)
e07139a v0.0.0.9023: nfold defaults to 3 on high-p autotune
a4623de v0.0.0.9022: drop fast.k=0 from balanced fk_grid on high-p
0045135 v0.0.0.9021: cap autotune nk-grid at 2x for high-p
e302aa5 docs: update CLAUDE.md handoff for v0.0.0.9020 + mlbench results
fa30989 v0.0.0.9020 (followup): tighten warmstart decisiveness rule
284635c v0.0.0.9020: Phase 3 - shared forward pass across autotune grid
a0e0f72 v0.0.0.9019: Phase 3 - autotune.warmstart (subsample pre-fit)
bb3f858 v0.0.0.9018: Phase 2 - n.boot bagging (earth has no bag)
4f48699 v0.0.0.9017: Phase 2 - autotune.speed knob
5983f30 v0.0.0.9016: Phase 2 - autotune nk grid + successive halving
1878772 v0.0.0.9015: Phase 2 - autotune (degree x penalty grid)
6aebff8 v0.0.0.9014: Phase 1 polish - 1-SE rule + per-fold CV diagnostics
29c27bc v0.0.0.9013: Phase 1 - pmethod="cv" with nfold + ncross
```

## Working conventions

- Standard devtools/usethis/roxygen2/testthat 3 chain.
- Always run `devtools::document()` after roxygen edits.
- Always run `devtools::test()` before commits to C++ or R/.
- Default `devtools::test()` is quick (~3s, 149 tests). It skips the
  heavy autotune + bagging + coverage suites via the `skip_if_quick()`
  helper. Run `Sys.setenv(ARES_FULL_TESTS=1); devtools::test()` (or
  `ARES_FULL_TESTS=1 Rscript -e devtools::test()`) to exercise the
  full 276-test suite (~5 min) before commits that touch autotune /
  bagging / CV / weights.
- Determinism is a hard invariant - any change that breaks
  `nthreads=1 == nthreads=N` byte-for-byte is a bug.

## Likely next sessions

1. **Multinomial family** (stretch). Either K-1 binomial via softmax
   decomposition OR `nnet::multinom` on the selected basis (would add
   a Suggests dep). Mostly mirror of the poisson / gamma plumbing.
2. **Multi-response gaussian** (item 4 from v0.0.0.9027 roadmap).
   Punted in v0.0.0.9027 due to KnotScanner / AVX2 surgery scope.
   Revisit if user demand materialises.
3. **Family-specific prediction intervals** for binomial / poisson /
   gamma. Currently `interval = "pint"` errors for non-gaussian
   families. The natural extensions are:
   - binomial: Wilson / Agresti-Coull on `plogis(eta +/- z*se)`.
   - poisson: Wald on log scale, exp() back to response.
   - gamma: log-normal-ish; needs shape estimate.
4. **`/contribute`** the 10 entries in `brain-contributions-v0.20.md`
   from a brain-connected session.
5. **Further high-p autotune speed** — remaining bottleneck is
   `auto_linpreds = FALSE` default at p=20 deg=2/3 inflating forward
   emit count; the v0.21-v0.23 gates attacked symptoms, not the
   default. Revisit if a user reports highdim wall-clock issues.
6. **GitHub publish** - user has not yet authorized remote push.
7. **Rename for CRAN** - `aresMARS` is the proposed name (avoids
   collision with archived `ARES` package).
