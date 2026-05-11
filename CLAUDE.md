# ares - handoff brief (post-tuning roadmap)

## What this package is

`ares` is an R package providing a fast Multivariate Adaptive Regression
Splines (MARS) implementation. Single user-facing function `ares()`. Built
with `Rcpp` + `RcppParallel` (TBB). Cross-platform. Deterministic across
threads (byte-identical fits at `nthreads=1` vs N).

**Author**: Jack Trametta <jtrametta@gmail.com>. **License**: MIT. **Repo**:
local git only on `main`.

`earth` (`earth::earth`) is a useful diagnostic but not a parity target.

## Current state - v0.0.0.9023

### Test suite
- 164/164 testthat green (was 38 at v0.0.0.9012). Determinism contract:
  `nthreads=1` and `nthreads=N` produce byte-identical fits at fixed seed.
- `R CMD check`: 0 errors, 3 warnings (Makevars GNU extensions, qpdf
  missing, system non-portable compile flags - all pre-existing,
  documented).

### Pruning options
- `pmethod = "backward"` (default) - GCV-based subset selection.
- `pmethod = "none"` - keep all forward-pass terms.
- `pmethod = "cv"` (v0.13+) - K-fold cross-validated subset selection.
  Promoted automatically when `nfold > 0`. Args: `nfold` (default 0;
  5 when pmethod="cv"), `ncross` (CV repetitions), `stratify`,
  `seed.cv`, `cv.1se` (1-SE rule).

### Autotune (`autotune = TRUE`)
- **Phase 2 (v0.15-v0.18)**: inner-CV grid search over
  `(degree, penalty, nk, fast.k)`. Cells scored via shared fold
  partition; per-cell mean holdout MSE. Successive halving drops
  cells whose fold-1 MSE is >1.5x the running best. Tie-break
  prefers parsimony.
- **Phase 3 (v0.19-v0.20)**: warm-start probes a 20% subsample
  (capped 200 rows) and short-circuits to the subsample winner if
  decisive (>=10% gap between best-per-degree). Shared-forward
  groups cells by `(degree, nk, fast.k)` and runs forward once per
  group per fold, replaying GCV-backward via
  `mars_backward_only_cpp` per cell.
- `autotune.speed`: `"balanced"` (default), `"quality"` (forces
  `fast.k = 0`), `"fast"` (forces `fast.k = 5`).
- `autotune.warmstart`: TRUE (default).
- `n.boot > 0`: row-bootstrap bagging. Composes with autotune (each
  replicate uses the central fit's tuned hyperparams).

### v0.21-v0.23 high-p autotune speedups
- v0.21: cap `nk_grid` at `2*nk_eff` when `nk_eff >= 31` (p >= 15).
- v0.22: drop `fast.k = 0` from balanced `fk_grid` on high-p (no-cache
  cell never sits within 1% of best).
- v0.23: `nfold` defaults to `3` (not 5) on high-p autotune; inner-CV
  variance unchanged, per-fold cost halved.
- Cumulative effect on highdim p=20 n=500: 138s -> 9.3s (~15x).
- Holdout MSE on highdim cells preserved (parsimony pick actually
  *improves* friedman1 p=20 n=500 by ~15%).

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

## Scope discipline (locked-in)

- Gaussian response only. No GLM, no multi-response, no obs weights.
- `pmethod` in `{"backward", "none", "cv"}`. No exhaustive/forward/seqrep.
- No `fast.k` heuristic past the default (10) is intended; user can
  override via `autotune.speed`.
- Deps: hard = `Rcpp`, `RcppParallel`. Suggests = `earth`, `testthat`,
  `bench`. **No tidyverse, no rlang, no S4, no caret/tidymodels.**
- No external tuning deps. License stays MIT.

## CRAN gotcha

Name collision with archived CRAN package `ARES`. Local install
unaffected. Before any CRAN submission, rename to `aresMARS`.

## Where things live

- `R/ares.R` - main function, formula + matrix interfaces, CV
  orchestrator (.ares_cv_fit), autotune orchestrator (.ares_autotune).
- `R/predict.R` - bagging-aware predict() with se.fit support.
- `R/print.R` - S3 methods.
- `src/ares.cpp` - C++ engine. Three exports:
  `mars_fit_cpp` (forward + backward),
  `mars_basis_cpp` (basis builder for predict),
  `mars_backward_only_cpp` (backward replay - v0.20 shared forward).
- `tests/testthat/` - 8 test files, 164 tests.
- `inst/sims/` - Monte Carlo + bench scripts. Latest results in
  `results/v0.23-mlbench.csv` (`v0.20-mlbench.csv` is the
  pre-Phase-4 reference).
- `brain-contributions-v0.20.md` (root, Rbuildignored) — 10 distilled
  knowledge entries proposed for the statsclaw brain seedbank.
  Awaiting `/contribute` upload from a brain-connected session.

## Commit history (post-tuning roadmap)

```
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
- Default `devtools::test()` is quick (~2s, 65 tests). It skips the
  heavy autotune + bagging suites via the `skip_if_quick()` helper.
  Run `Sys.setenv(ARES_FULL_TESTS=1); devtools::test()` (or
  `ARES_FULL_TESTS=1 Rscript -e devtools::test()`) to exercise the
  full 164-test suite (~5 min) before commits that touch autotune /
  bagging / CV.
- Determinism is a hard invariant - any change that breaks
  `nthreads=1 == nthreads=N` byte-for-byte is a bug.

## Likely next sessions

1. **`/contribute`** the 10 entries in `brain-contributions-v0.20.md`
   from a brain-connected session.
2. **Further high-p autotune speed** — remaining bottleneck is
   `auto_linpreds = FALSE` default at p=20 deg=2/3 inflating forward
   emit count; the v0.21-v0.23 gates attacked symptoms, not the
   default. Revisit if a user reports highdim wall-clock issues.
3. **GitHub publish** - user has not yet authorized remote push.
4. **Rename for CRAN** - `aresMARS` is the proposed name (avoids
   collision with archived `ARES` package).
