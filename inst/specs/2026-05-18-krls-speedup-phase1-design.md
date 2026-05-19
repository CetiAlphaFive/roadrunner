# Design — roadrunner::krls() speedup, Phase 1: shared-distance + parallel autotune

**Date**: 2026-05-18
**Status**: Design approved (brainstorming session)
**Target version**: v0.0.0.9042
**Repo HEAD at design time**: `2c83d4a` (v0.0.0.9041)

---

## 1. Motivation

`roadrunner::krls()` is currently O(n³) eigendecomposition × 9 sigma cells × 10 folds × 2 ncross = 180 sigma-cells per autotune fit. At n=500 p=10 a single autotune fit takes ~6 s. At n=1500 p=20 it becomes painful.

The head-to-head sweep in REQ-20260518-003 confirmed correctness (12-0-3 vs `KRLS::krls`). Now we go after wall-clock without changing fits.

Three observations drive Phase 1:
- The pairwise squared-distance matrix `D` is independent of sigma. Current code recomputes it for every sigma cell.
- The 9 sigma cells per (ncross, fold) are mutually independent and CPU-bound.
- The TBB backend used by `ares()` is already linked; reusing it requires no new dependency.

Phase 2 (Nystrom approximation, opt-in default at n > threshold) is a separate design and ship.

---

## 2. Goals and non-goals

### Goals
- Median wall-clock speedup ≥ **3.0×** vs v0.0.0.9041 sequential on `autotune=TRUE` fits at (n=500, p=10) and (n=1500, p=20).
- Preserve byte-identical fits at fixed seed across `autotune.nthreads ∈ {1, ..., N}`.
- Preserve all current API contracts. No new required arguments.
- Cross-platform (Linux, macOS, Windows) via RcppParallel.

### Non-goals (Phase 1)
- Nystrom approximation.
- Random Fourier features.
- Predict()-side speedups.
- Multi-response / sparse kernel paths.
- Changes to the non-autotune fit path (`autotune=FALSE`).
- Changes to `lambda.method='cv'` outside of autotune.

---

## 3. Architecture

### Current autotune flow
For each (sigma in autotune.grid) × (cross in 1..ncross) × (fold in 1..nfold):
```
build K(sigma) [O(n²p)]  →  eigen(K) [O(n³)]  →  lambda LOO search  →  MSE on fold
```
180 cells. K rebuilt from raw X each cell.

### New flow
```
for cross in 1..ncross:
  for k in 1..nfold:
    D_tr = pairwise_sqdist(X_train_k, X_train_k)         # ONCE, [O(n_tr² p)]
    D_te = pairwise_sqdist(X_test_k,  X_train_k)         # ONCE, [O(n_tr n_te p)]

    parallel for sigma in sigma_grid (TBB):
      K_tr = exp(-D_tr / sigma)                          # [O(n_tr²)] no p factor
      K_te = exp(-D_te / sigma)                          # [O(n_tr n_te)]
      eig(K_tr) -> (V, d)                                # [O(n_tr³)]  hot loop
      lambda  = lambdasearch(V, d, y_tr)                 # tol=1e-6, log-step
      alpha   = solve(V, d, lambda, y_tr)
      yhat_te = K_te @ alpha
      mse_arr[cross, k, sigma_idx] = mean((y_te - yhat_te)^2)
      lam_arr[cross, k, sigma_idx] = lambda
```

After the per-fold sweep, aggregate `mse_arr` and apply the existing 1-SE rule (selecting the largest sigma within 1 SE of the minimum). Refit at the chosen sigma using the existing non-autotune path (zero change there).

### Key invariants
- `D` matrix computed once per fold partition, reused across all sigmas.
- Eigendecomposition still happens per (fold, sigma) cell. No algebraic shortcut across sigma exists.
- Parallel scope is the sigma dimension only. Folds and ncross stay sequential (cache locality + clean reduction).
- Public API unchanged; `fit` object structure identical to v0.0.0.9041.

---

## 4. API surface

### Public
Zero breaking changes.

### One new optional knob (additive)
- `autotune.nthreads` (integer, default `getOption("roadrunner.nthreads", parallel::detectCores(logical=FALSE))`). Caps TBB worker count for the sigma-grid parallel pass.
  - `1` → strictly sequential. Matches pre-patch behavior byte-for-byte at fixed inputs.
  - `>= 2` → TBB parallel over sigma.
  - Mirrors `ares()`'s existing `nthreads` arg.

### Internal C++ exports — new (in `src/krls.cpp`)
- `krls_pairwise_sqdist_cpp(X_a, X_b)` → `(nrow(X_a), nrow(X_b))` numeric matrix.
- `krls_autotune_inner_cpp(D_tr, D_te, y_tr, y_te, sigma_grid, lambda_args, nthreads)` → list with `mse_per_sigma`, `lambda_per_sigma`.

### Internal C++ exports — unchanged
`krls_kernel_cpp`, `krls_kernel_pred_cpp`, `krls_eig_cpp`, `krls_solve_cpp`, `krls_loo_loss_cpp`, `krls_deriv_cpp`, `krls_avg_deriv_var_cpp`, `krls_vsq_cpp` — all preserved.

### R-level
- `.ares_autotune` body refactored. Loops over (ncross, fold) sequentially; each fold calls `krls_autotune_inner_cpp` once.
- New helper `.krls_call_autotune_inner(...)` wraps the C++ call with validation.
- `autotune_info` return gains `nthreads_used` and `sigma_grid_sorted`.

### Roxygen
- New `@param autotune.nthreads`.
- `@details` block mentions parallel sigma sweep + determinism guarantee.
- `@note` documents `getOption("roadrunner.nthreads", ...)` resolution.

### No new dependencies
No new `Imports` or `Suggests`. RcppParallel/TBB already linked.

---

## 5. C++ internals

### `krls_pairwise_sqdist_cpp(X_a, X_b)` (~30 LOC)
- Returns squared Euclidean distance matrix `D[i, j] = ||X_a[i, ] - X_b[j, ]||²`.
- Identity: `D = ||X_a||² + ||X_b||²ᵀ - 2 X_a X_bᵀ`.
- Implementation: single `dgemm` via Armadillo `X_a * X_b.t()` plus two `sum-of-squares` row vectors.
- Symmetric shortcut: if `&X_a == &X_b`, skip second sum.
- Numerical floor: clamp tiny negatives to 0 (FP rounding produces `-1e-15` on diagonal).

### `krls_autotune_inner_cpp(...)` (~180 LOC)
- Outer: `RcppParallel::parallelFor` over `[0, nsigma)` with grain size 1.
- Per worker, per sigma s:
  1. `K_tr = exp(-D_tr / s)` (n_tr × n_tr)
  2. `K_te = exp(-D_te / s)` (n_te × n_tr)
  3. `arma::eig_sym(K_tr) -> (V, d)` — deterministic.
  4. Cache `Vsq = V % V` and `Vty = V.t() * y_tr`.
  5. Golden-section LOO lambda search inlined (same formula as `krls_loo_loss_cpp`, tol=1e-6, log-step L).
  6. Solve for alpha; `yhat_te = K_te * alpha`; store `mse_te(s)` and `lambda(s)`.
- Worker output: per-sigma slot in pre-allocated `arma::vec mse_te(nsigma)` and `arma::vec lambda(nsigma)`. Each thread writes a unique slot. No reduction needed → race-free → byte-identical result regardless of `nthreads`.

### Determinism strategy
- No random draws inside the parallel region.
- BLAS thread contention is mitigated by structure: the only large `dgemm` call (`krls_pairwise_sqdist_cpp`) runs OUTSIDE the parallel region. Inside, the heavy op is `arma::eig_sym` (LAPACK `dsyevr`), which is single-threaded by default. Tiny elementwise `exp(-D/s)` and matrix-vector products do not spawn BLAS threads.
- If a user has set process-global BLAS threads high (e.g. `RhpcBLASctl::blas_set_num_threads(8)`), recommend setting to 1 around `krls()` autotune calls. Documented in `@note`.
- `lambda_args` carries the bracket logic (tol=1e-6, log-step L). Within-sigma lambda is identical to v0.0.0.9041 sequential output.
- Eigendecomposition is deterministic on fixed input (LAPACK `dsyevr` is deterministic).
- Test: `expect_identical(autotune_nthreads_1, autotune_nthreads_4)` on full grid output.

### Per-sigma cost (n_tr=450, n_te=50 at nfold=10)
- Kernel build: O(n_tr²) ≈ 200K flops (down from O(n_tr² p) ≈ 4M).
- Eigen: O(n_tr³) ≈ 90M flops — still hot.
- LOO lambda search: ~25 iters × O(n_tr) per iter — cheap.
- → Eigen dominates. Distance-reuse alone ~1.2× speedup; parallel sigma stacks on top.

### Real speedup estimate
| n | p | nthreads | est. speedup |
|---|---|---|---|
| 500 | 10 | 1 | 1.2× (distance reuse) |
| 500 | 10 | 4 | 3.5–4× |
| 500 | 10 | 8 | 5–6× (9 sigmas, diminishing returns) |
| 1500 | 20 | 4 | 3–4× |

---

## 6. R harness changes

### File and scope
`R/krls.R`, `.ares_autotune` body (~lines 463-700 in v0.0.0.9041). ~80 LOC delta.

### New loop
```r
sigma_grid_sorted <- sort(autotune.grid)
mse_arr <- array(NA_real_, c(ncross, nfold, length(sigma_grid_sorted)))
lam_arr <- array(NA_real_, c(ncross, nfold, length(sigma_grid_sorted)))

for (cross in seq_len(ncross)) {
  folds <- make_folds(n, nfold, seed = cross_seed[cross])
  for (k in seq_len(nfold)) {
    idx_te <- folds[[k]]
    X_tr   <- X[-idx_te, , drop = FALSE]
    X_te   <- X[ idx_te, , drop = FALSE]
    y_tr   <- y[-idx_te]
    y_te   <- y[ idx_te]

    D_tr <- krls_pairwise_sqdist_cpp(X_tr, X_tr)
    D_te <- krls_pairwise_sqdist_cpp(X_te, X_tr)

    inner <- krls_autotune_inner_cpp(
      D_tr, D_te, y_tr, y_te,
      sigma_grid_sorted, lambda_args,
      nthreads = autotune.nthreads
    )

    mse_arr[cross, k, ] <- inner$mse_per_sigma
    lam_arr[cross, k, ] <- inner$lambda_per_sigma
  }
}

mse_per_sigma <- apply(mse_arr, 3, mean, na.rm = TRUE)
se_per_sigma  <- apply(mse_arr, 3, sd,   na.rm = TRUE) / sqrt(ncross * nfold)
min_idx       <- which.min(mse_per_sigma)
threshold     <- mse_per_sigma[min_idx] + se_per_sigma[min_idx]
sigma_1se_idx <- max(which(mse_per_sigma <= threshold))
sigma_chosen  <- sigma_grid_sorted[sigma_1se_idx]
```

### Validation block
- `autotune.nthreads` integer ≥ 1; default `getOption("roadrunner.nthreads", parallel::detectCores(logical=FALSE))`.
- Clamp upper at `length(sigma_grid)`.
- If `RcppParallel::defaultNumThreads()` returns 1, force `autotune.nthreads = 1` and emit a one-time message.

### Refit at chosen sigma
Unchanged path: `krls_kernel_cpp` → `krls_eig_cpp` → existing `.krls_lambdasearch` → `krls_solve_cpp`.

### Unchanged
`krls.formula`, `predict.krls_rr`, `summary.krls_rr`, `print.krls_rr`, derivative path, S3 class, `lambda.method='cv'` non-autotune path.

---

## 7. Testing

### New file
`tests/testthat/test-krls-autotune-parallel.R`. ~12 tests.

### Correctness (parity with v0.0.0.9041 sequential)
1. **EQUIV-1** — byte-identical at `nthreads=1` vs snapshot baseline in `inst/testdata/autotune-baseline-9041.rds`.
2. **EQUIV-2** — `nthreads=1 ≡ nthreads=4` byte-identical (determinism contract).
3. **EQUIV-3** — shuffled `autotune.grid` input → identical fit.
4. **EQUIV-4** — (ncross × nfold) ∈ {(2,2), (2,5), (3,10)} all match baselines.

### Distance helper
5. **DIST-1** — matches `as.matrix(dist(X))^2` to tol 1e-10.
6. **DIST-2** — non-square cross-matrix matches naive double loop.
7. **DIST-3** — numerical floor (no negative leakage on near-zero rows).

### Parallel C++ inner
8. **INNER-1** — vs R-level reference reimplementation, `mse_per_sigma` and `lambda_per_sigma` match to tol 1e-9.
9. **INNER-2** — `sigma_grid` length 1 still works.
10. **INNER-3** — `nthreads > nsigma` clamps cleanly.

### Regression
11. **REG-1** — full krls suite `devtools::test(filter="krls")` 0 fails.
12. **REG-2** — `KRLS::krls` parity at matched (sigma, lambda), tol 1e-9 (existing test).

### Empirical re-validation — EMP-PHASE1
Sim script `inst/sims/krls-speedup-phase1.R`. **Mandatory sim rules (this session)**:
- Smoke test: R=2, 1 DGP, both `nthreads ∈ {1, 4}`; assert finite + identical output.
- Intermediate CSV per cell.
- Wall-clock cap 5 min.

Cells:
| n | p | nthreads | baseline (v0.0.0.9041) | new |
|---|---|---|---|---|
| 500 | 10 | 1 | ✓ | ✓ |
| 500 | 10 | 4 | — | ✓ |
| 1500 | 20 | 1 | ✓ | ✓ |
| 1500 | 20 | 4 | — | ✓ |

R=5 reps per cell. Metrics: wall-clock, test_MSE (identical expected vs v0.0.0.9041).

### Acceptance
- Median speedup at (500, 10, 4) ≥ **3.0×** vs v0.0.0.9041 sequential.
- Median speedup at (1500, 20, 4) ≥ **3.0×**.
- `test_MSE_new == test_MSE_baseline` byte-identical when `nthreads = 1`.
- `|test_MSE_new(nthreads=4) - test_MSE_new(nthreads=1)| < 1e-9`.

### R CMD check
0E / 0W / no new notes (2 pre-existing pass through).

---

## 8. Rollout

### Version and commit shape
v0.0.0.9042. Single commit covering:
- `src/krls.cpp` (~210 LOC new)
- `R/krls.R` (~80 LOC delta in `.ares_autotune` + roxygen)
- `tests/testthat/test-krls-autotune-parallel.R` (new, ~12 tests)
- `inst/sims/krls-speedup-phase1.R` (new)
- `NEWS.md` entry
- `DESCRIPTION` bump
- `man/krls.Rd` regen
- `inst/testdata/autotune-baseline-9041.rds` (snapshot for EQUIV-1)

### Pipeline
Full statsclaw workflow — planner → builder ∥ tester → scriber → reviewer → shipper.

### Definition of done
- Commit shipped, pushed to `origin/main`.
- `NEWS.md` v0.0.0.9042 entry with speedup table.
- `inst/sims/results/krls-speedup-phase1.csv` checked in.
- Open issue or note in NEWS.md for Phase 2 (Nystrom) follow-up.

---

## 9. Risks

| Risk | Likelihood | Mitigation |
|---|---|---|
| TBB nesting with BLAS inside eigen | medium | Set BLAS threads = 1 inside parallel block via `RcppParallel::TBBControl` or `omp_set_num_threads(1)` in worker entry |
| FP non-determinism from BLAS thread races | medium | Force BLAS = 1 thread. Verify via EQUIV-2 test |
| Eigen result differs across Armadillo versions | low | Pin via `Imports: RcppArmadillo (>= 0.12)`; snapshot test against installed version |
| Windows TBB linkage | low (existing `ares` works) | RcppParallel handles cross-platform |
| `autotune.nthreads` default surprises users on shared servers | medium | Default via `getOption("roadrunner.nthreads", ...)`; documented in `@note` |
| Memory blow-up: `D_tr` is `n² × 8` bytes (n=5000 → 200MB) | low for current scope | Sized for n ≤ 2000 (32MB). Larger n flagged with one-time message suggesting Nystrom Phase 2 |
| Stale autotune snapshot baselines | high | Refresh `autotune-baseline-9041.rds` BEFORE building, freeze at HEAD `2c83d4a` |

---

## 10. Out of scope (Phase 2, separate ship)

- Nystrom approximation (`approx = "nystrom"`, `nystrom_m`, `landmarks` args)
- Random Fourier features
- Predict-side speedups
- Multi-response / sparse kernel paths
- C++ rewrite of the non-autotune fit path
- Replacing the eigen path with a direct ridge solve (algorithmic change, deferred)

---

## Approval

Approved during brainstorming session 2026-05-18. Next step: invoke `superpowers:writing-plans` to produce the implementation plan.
