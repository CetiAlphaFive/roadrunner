# ares — handoff brief

## What this package is

`ares` is an R package providing a fast Multivariate Adaptive Regression Splines (MARS) implementation. Single user-facing function `ares()` mirrors `earth::earth()` from the `earth` CRAN package. Built with `Rcpp` + `RcppParallel` (TBB). Cross-platform. Deterministic across threads (byte-identical fits at `nthreads=1` vs N).

**Author**: Jack Trametta <jtrametta@gmail.com>. **License**: MIT. **Repo**: local git only on `main` (no GitHub remote yet).

**Ground truth**: `earth` package. Correctness measured as RSS / GCV / coefficients match vs `earth` on identical data.

## Current state — v0.0.0.9002

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

## v0.2 progress (current)

- Class (e) blowup fixed: forward pass refuses ill-conditioned hinge pairs (RSS guard at 1e6 × rss0 + non-finite check). Multi-seed sweep on interaction n=5000 deg=2 confirms blowup is gone.
- Backward pass parallelised (`RcppParallel::parallelFor` over M-1 candidate-drop trials). Threading scaling at n=3000 deg=2 recovered from 1.16× → 1.85× (2t/1t).
- Median 1t wall-clock vs earth: ~27× (basically unchanged from v0.1; backward parallelism is a multi-thread win, not a single-thread win).
- Median 4t wall-clock: 0.24s vs earth's 0.020s ≈ 12× ratio.
- **Neither user objective met yet:** complete-parity (no divergences) and "must-be-faster-than-earth" both open.
- Closed-form Cholesky downdate ΔRSS_j = β_j² / [(R'R)^{-1}]_jj was tried and reverted: numerically unstable on the rank-deficient designs in the bench grid. Helpers `householder_qr_R / solve_R_back / chol_diag_inverse / qr_downdate_col` are still in `src/ares.cpp` (lines ~399-510) and ready for v0.3 when paired with a pivoted/regularised factor.
- See `inst/sims/results/divergence-diag.md` and `inst/sims/results/friedman-spec.md` (both produced by v0.2 planning subagents) for outstanding parity gaps.

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
