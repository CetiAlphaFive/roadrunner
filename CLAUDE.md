# ares - handoff brief (post-tuning roadmap)

## What this package is

`ares` is an R package providing a fast Multivariate Adaptive Regression
Splines (MARS) implementation. Single user-facing function `ares()`. Built
with `Rcpp` + `RcppParallel` (TBB). Cross-platform. Deterministic across
threads (byte-identical fits at `nthreads=1` vs N).

**Author**: Jack Trametta <jtrametta@gmail.com>. **License**: MIT. **Repo**:
local git only on `main`.

`earth` (`earth::earth`) is a useful diagnostic but not a parity target.

## Current state - v0.0.0.9020

### Test suite
- 154/154 testthat green (was 38 at v0.0.0.9012). Determinism contract:
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

### v0.20 mlbench-style benchmark (vs earth, ranger)
inst/sims/results/v0.20-mlbench.csv. Highlights:
- ares default beats earth on 6/10 cells, ties most others.
- autotune wins on highdim n=500 (1.059 vs earth 1.242, ranger 1.711)
  and n=1500 (0.946 vs 1.103, 1.556) by selecting deg=2 / deg=3.
- ranger is ~3-4x worse on every cell - ares is the right MARS pick
  for continuous structured DGPs.
- Autotune wall-clock: ~0.1-1.4s for friedman1/2/additive/interaction
  at typical n. Highdim p=20 hits ~60-120s - degree=3 grid is heavy.

### Determinism gotcha (still in place)
- The `-ffp-contract=off` flag was dropped at v0.0.0.9012 (FMA
  enabled). Run-to-run fits are still bit-identical.

## Known issue - speed target only met conditionally

- ares default vs earth: median ratio ~0.93x at 4t (ares slightly
  faster). 11/18 cells in inst/sims faster than earth at 4t.
- Autotune wall-clock: sub-second when warm-start fires; ~1-3s on
  small n full grid; up to 100s on highdim full grid.
- Sub-second autotune target met for warm-start path; full-grid
  on large p needs further work (out of scope of v0.20 roadmap).

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
- `tests/testthat/` - 8 test files, 154 tests.
- `inst/sims/` - Monte Carlo + bench scripts.
  Latest results in `results/v0.20-mlbench.csv`.

## Commit history (post-tuning roadmap)

```
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
- Determinism is a hard invariant - any change that breaks
  `nthreads=1 == nthreads=N` byte-for-byte is a bug.

## Likely next sessions

1. **Speed up large-p autotune** - highdim p=20 deg=3 grid is slow.
   Options: drop deg=3 on noisy data, parallelize cells across
   threads.
2. **Brain contribution upload** - `brain-contributions.md` from the
   v0.0.0.9000 session has 5 entries pending.
3. **GitHub publish** - user has not yet authorized remote push.
4. **Rename for CRAN** - only when ready to submit.
