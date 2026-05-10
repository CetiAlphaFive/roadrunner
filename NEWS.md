# ares 0.0.0.9000

Initial development release. Not yet on CRAN.

## What works

- `ares()` fits gaussian-only MARS via forward hinge selection and backward
  GCV-based subset pruning.
- Matrix and formula interfaces.
- Standard S3 methods: `predict.ares`, `print.ares`, `summary.ares`, `plot.ares`.
- Forward-pass parallelism via **RcppParallel** (`parallelFor` over
  parent × variable pairs).
- **Threading is deterministic** — `nthreads = 1` and `nthreads = N` produce
  byte-identical fits.
- **Numerical parity with earth** — across 24-cell smoke Monte Carlo, ares
  reproduces earth's RSS to ~1% (median) on Friedman-1, additive, and
  interaction DGPs.
- 38 testthat tests pass; `devtools::check()` clean.

## Known limitations (v0.0.0.9000)

- **Wall-clock vs `earth`**: ares is slower than earth in absolute time
  (typically 30–80× slower at n ∈ {500, 1500}). Earth uses Friedman's O(1)
  per-knot Givens fast-LS update; ares v0.0.0.9000 uses an O(n) per-knot
  scoring loop. Implementing the Givens fast-LS inner loop is a **v0.1
  milestone**. Parallel scaling itself is healthy: 2-thread is ~1.75×
  faster than 1-thread.
- Gaussian responses only (no GLM, no multi-response, no observation weights).
- `pmethod` supports only `"backward"` and `"none"`. No exhaustive,
  cross-validated, or sequential-replacement pruning.
- No `fast.k` / `fast.beta` heuristics yet.

## Roadmap

- v0.1: Friedman fast-LS Givens-rotation inner loop (target: ≥ 1.5× single-thread
  vs earth, ≥ 4× multi-thread on n ≥ 5000).
- v0.2: Cross-validated pruning (`pmethod = "cv"`).
- v0.3: GLM responses (binomial, poisson, gaussian-GLM).
- v0.4: Observation weights and stratified CV.
