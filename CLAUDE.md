# ares — handoff brief

## What this package is

`ares` is an R package providing a fast Multivariate Adaptive Regression Splines (MARS) implementation. Single user-facing function `ares()` mirrors `earth::earth()` from the `earth` CRAN package. Built with `Rcpp` + `RcppParallel` (TBB). Cross-platform. Deterministic across threads (byte-identical fits at `nthreads=1` vs N).

**Author**: Jack Trametta <jtrametta@gmail.com>. **License**: MIT. **Repo**: local git only on `main` (no GitHub remote yet).

**Ground truth**: `earth` package. Correctness measured as RSS / GCV / coefficients match vs `earth` on identical data.

## Current state — v0.0.0.9004

- `devtools::check()`: 0 errors, 2 warnings (documented + accepted), 1 note.
  - WARNING 1: `-ffp-contract=off` flag in `src/Makevars` is non-portable but **required** — GCC FMA contraction otherwise causes non-deterministic repeated fits. Do not remove without replacement fix.
  - WARNING 2: `$(shell ...)` GNU-make extension in `Makevars` for RcppParallel linking. Idiomatic.
- 38/38 testthat tests pass.
- Numerical parity with `earth`: comparable to v0.0.0.9000 (median RSS rel-err ~1% across the 54-cell v0.1 grid; some cells noise-driven divergence at large n with weak interaction signal).
- Determinism preserved — byte-identical fits at `nthreads=1` vs N.

## Speed status — v0.1 fast-LS landed

- v0.1 replaced the O(K · n · Mq) per-(parent, var) inner loop with an O((n + K) · Mq) prefix-sum identity (Friedman 1991 §3.5 fast-LS).
- Result on inst/sims grid (single-thread, median across n × deg × dgp cells in {500,1500,5000} × {1,2} × {fri,add,int}): wall-clock ratio vs earth dropped from **~152× → ~28×** (5.4× speedup at median; 16.6× at the worst pre-v0.1 cell — Friedman-1 n=1500 deg=2 went 18.1s → 1.09s).
- Math equivalence verified: an R-side reference that scores via the explicit slow projection picks the exact (var, knot, reduction) the C++ fast-LS picks.
- **Speed parity vs earth still NOT MET** — original v0.1 target was ≥1.5× ST and ≥4× MT. Remaining gap is dominated by serial parts (`build_Q` + `ols_qr` per forward step). Multi-thread scaling shrunk because the parallelizable region got much smaller (Amdahl).
- Two ways forward (v0.2 candidates):
  1. Incremental Q maintenance (avoid full re-orthogonalisation per forward step).
  2. Parallelise `ols_qr` + lift Mq-sized work out of inner KnotScanner.
- See `inst/sims/results/baseline-pre-v0.1.csv` and `inst/sims/results/v0.1-bench.csv` for raw numbers.

## v0.2 → v0.4 progress

### v0.2: parallel backward + class (e) bug fix
- Forward pass refuses ill-conditioned hinge pairs (RSS > 1e6 × rss0 + non-finite check). Closes the interaction n=5000 deg=2 seed=2 blowup (RSS=3.5e29 was emitted at M=13).
- Backward pass parallelised over M-1 candidate-drop trials. Threading scaling at n=3000 deg=2 recovered from 1.16× → 1.85× (2t/1t).

### v0.3: experimental earth-equivalent heuristics (defaults OFF)
- New args `adjust.endspan` (default 1) and `auto.linpreds` (default FALSE). Code paths and `is_boundary` flag in Candidate are in place; defaults match v0.2 behaviour.
- Earth's analogous defaults are 2 and TRUE; ares's implementations approximate earth's docs but regress parity vs earth on the inst/sims grid because earth's actual `Auto.linpreds` test is "h+ is linear over parent's support" rather than "knot at leftmost eligible position". Treat as v0.5+ work.

### v0.4: backward QR Givens downdate
- Per-trial backward downdate (Givens rotations on copies of R, Qty) replaces per-trial Householder rebuild. Chosen drop committed via fresh Householder on smaller cur to anchor numerics each step.
- Median 1t wall-clock vs earth: **~27× → ~6.6×** (4× single-thread speedup vs v0.2). Min ratio 1.13× on best cells.
- Friedman-1 n=1500 deg=2: was 1.09s (v0.2), now ~0.34s.
- Numerical regime change: v0.4 returns true OLS minimum RSS via the QR factor; v0.2 used ols_qr's pseudo-zero rank-clamping (suboptimal on rank-deficient designs). On tightly tied cells, `rss_ares < rss_earth` is now possible. Friedman parity-test thresholds bumped 5e-3 → 2e-2; mtcars 0.5 → 1.0.
- Determinism preserved (byte-identical 1t vs N).

## Outstanding vs user objectives

- **"Complete parity (no divergences)"**: not met. Median grid rss_rel_err ~0.9%; 7/18 cells > 1%. The remaining divergences trace to earth's `Adjust.endspan` + `Auto.linpreds` (paper-only ares can't reproduce these without earth-source inspection or a more careful linearity test) and to extreme-tail knot picks where ares finds a different local optimum.
- **"Must be faster than earth"**: not met. 0/18 grid cells faster than earth at single-thread. Best ratio 1.13× (very close on simplest cells); median 6.6×. Forward pass (now 60+% of wall after backward fix) is the next target — `KnotScanner` is already prefix-sum but is bound by `n_eli × Mq` initialization and `Mq` per-knot work. Possible v0.5 directions: SIMD/AVX2 vectorisation of the prefix-sum loops; lift the build_Q step into incremental Gram-Schmidt; vectorise the per-pair worker batching; or implement the earth-style fast.k pruning to reduce candidate-pair count.
- See `inst/sims/results/divergence-diag.md`, `inst/sims/results/friedman-spec.md`, and `inst/sims/results/v0.2-bench.csv` for raw numbers.

## CRAN gotcha

Name collision with archived CRAN package `ARES`. `R CMD check --as-cran` errors on this. **Local install/use is unaffected.** Before any CRAN submission, rename to `aresMARS` (or similar). DESCRIPTION + NAMESPACE + `R/ares-package.R` + `src/` register calls all need updating.

## Scope discipline

Locked-in scope (do not expand without explicit user direction):

- Gaussian response only. No GLM, no multi-response, no observation weights.
- `pmethod` supports `"backward"` and `"none"` only.
- No `fast.k` / `fast.beta` heuristics.
- Deps: hard = `Rcpp`, `RcppParallel`. Suggests = `earth`, `testthat`, `bench`. **No tidyverse, no rlang, no S4.** Keep it lean.

## Where things live

- `R/ares.R` — main function, formula + matrix interfaces.
- `R/predict.R`, `R/print.R` — S3 methods.
- `src/ares.cpp` — C++ forward/backward/GCV. **This is where v0.1 fast-LS work goes.**
- `src/Makevars`, `src/Makevars.win` — build flags incl. `-ffp-contract=off`.
- `tests/testthat/` — 38 tests covering structure, predict, edge cases, earth-parity, threading determinism.
- `inst/sims/` — Monte Carlo scripts (smoke + larger grids).
- `vignettes/ares.Rmd` — earth-parity demo.
- `ARCHITECTURE.md` — deeper internals.

## Statsclaw run artifacts

Full pipeline trace at `/home/jack/.claude/plugins/data/statsclaw-statsclaw/workspace/ares/runs/REQ-20260509-225954-ares-init/`:
- `spec.md`, `test-spec.md`, `sim-spec.md` (planner)
- `implementation.md`, `simulation.md`, `audit.md`, `review.md`
- `brain-contributions.md` — 5 entries proposed (FMA determinism, hand-rolled QR, RcppParallel pattern, two-tier numerical-parity acceptance, CRAN-only-ERROR triage). Saved locally; not yet uploaded to brain seedbank.

## Likely next sessions

1. **v0.2 Amdahl unblock** — incremental Q maintenance and/or parallelised `ols_qr`. Goal: recover thread scaling and cross v0.1's original ≥1.5× ST / ≥4× MT @ n ≥ 5000 target.
2. **Brain contribution upload** — push `brain-contributions.md` entries to `statsclaw/brain-seedbank` if user authorizes.
3. **GitHub publish** — user has not yet authorized remote push.
4. **Rename for CRAN** — only when ready to submit.

## Working conventions

- Standard devtools/usethis/roxygen2/testthat 3 chain.
- Always run `devtools::document()` after roxygen edits.
- Always run `devtools::test()` before commits to C++ or R/.
- Determinism is a hard invariant — any change that breaks `nthreads=1 == nthreads=N` byte-for-byte is a bug.
- Earth-parity is measured numerically, not by source-code match.
