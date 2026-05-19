# Design — roadrunner::krls() Phase 2: Nystrom low-rank approximation

**Date**: 2026-05-19
**Status**: Design approved (brainstorming session)
**Target version**: v0.0.0.9043
**Repo HEAD at design time**: `e474b83` (v0.0.0.9042)
**Reference impl**: `/tmp/krls-reference/R/nystrom.R` (j-hai/KRLS, 380 LOC pure-R)

---

## 1. Motivation

Phase 1 (v0.0.0.9042) added a C++ autotune inner with byte-identical determinism across `nthreads`, but empirical wall-clock speedup capped at ~1× on multi-threaded BLAS systems. At n ≥ 800, `arma::eig_sym` on the full n×n kernel matrix saturates all cores via MKL/OpenBLAS, and adding TBB-level parallelism over sigma cells just redistributes the same CPU budget.

The real bottleneck is the O(n³) eigendecomposition itself, not parallelism. Nystrom approximation attacks it directly: replace the full kernel `K` (n×n) by a low-rank approximation anchored at m << n landmarks. The eigendecomposition then lives in m-space:

```
Phi      = C @ U @ diag(D_reg^{-1/2})           (n x m feature map)
beta_hat = (Phi' Phi + lambda I)^{-1} Phi' y    (m-dim ridge)
alpha    = U @ (D_reg^{-1/2} * beta)            (m-dim "dual" coefficient)
f(x)     = K(x, Z) @ alpha                      (predict via cross-kernel)
```

where C = K(X, Z), W = K(Z, Z), and (U, D) is the eigendecomposition of the m×m W matrix.

Cost: O(m³) eigen + O(nm²) feature map + O(nm) per LOO evaluation. At m = ⌈√n × 2⌉, e.g. n=2000 → m=90 → m³ ≈ 730K flops vs n³ = 8B flops. ~10,000× algorithmic reduction at n=2000.

This is not bottlenecked by BLAS; even fully MKL-saturated, the algorithm is faster because there's simply less work.

---

## 2. Goals and non-goals

### Goals
- Single-fit speedup ≥ **5×** at n=2000 p=10 vs `approx = "exact"` (v0.0.0.9042).
- Single-fit speedup ≥ **10×** at n=5000 p=10 vs exact.
- Autotune speedup ≥ **3×** at n=2000 with the default 9-sigma grid.
- Predict-time speedup ≥ **3×** at n_test=1000 once fit (cheaper kernel product).
- Test-set RMSE within **+10%** of exact at default m = ⌈√n × 2⌉ on the 15-DGP head-to-head suite.
- Determinism contract: at fixed `landmark_seed`, fits are byte-identical across `autotune.nthreads`.
- Zero breaking API changes. `krls()` without new args defaults to exact and behaves identically to v0.0.0.9042.

### Non-goals (Phase 2)
- Auto-engaging Nystrom at large n (opt-in only).
- Replacing exact mode (it stays the default).
- Multi-response / sparse kernel paths.
- Online / streaming landmark updates.
- Leverage-score-based landmark selection (random + kmeans only).
- Bagging integration (autotune + Nystrom, plus bagging, is a Phase 3 question).

---

## 3. API surface

### Public — additive, no breaking changes

`krls(X, y, ...)` gains five new optional arguments:

| Argument | Type | Default | Purpose |
|---|---|---|---|
| `approx` | character | `"exact"` | One of `"exact"` (Phase 1 behavior) or `"nystrom"` |
| `nystrom_m` | integer or NULL | NULL → `ceiling(sqrt(n) * 2)` | Number of landmarks |
| `landmarks` | NULL, integer vec, or numeric matrix | NULL | User-supplied landmarks. Integer = row indices into X; matrix = explicit landmark coordinates in original X scale |
| `landmark_method` | character | `"random"` | One of `"random"` (uniform subsample) or `"kmeans"` (Hartigan-Wong on standardized X). Ignored when `landmarks` is supplied |
| `landmark_seed` | integer or NULL | NULL | Seed for the local landmark draw. NULL = use global RNG state. When supplied, uses `.with_seed()` to avoid disturbing caller's RNG. |
| `nystrom_eps` | numeric | `1e-9` | Relative ridge floor on landmark-kernel eigenvalues. Stabilizes the m×m eigendecomp |

### S3 class

Fit object stays `class = "krls_rr"`. New `fit$approx` field set to `"exact"` or `"nystrom"`. Nystrom fits additionally carry:
- `fit$landmarks` — m×d standardized landmark coordinates
- `fit$landmark_indices` — integer indices into training X if landmarks came from rows, else NULL
- `fit$landmark_method` — `"random"`, `"kmeans"`, `"user_indices"`, or `"user_matrix"`
- `fit$nystrom_m` — integer m actually used
- `fit$nystrom_eps` — the relative ridge floor value used
- `fit$W_eigen` — list with `values` and `vectors` for the m×m anchor eigendecomp
- `fit$Dinvsqrt` — length-m vector for the feature-map scaling
- `fit$Sigma2` — length-m squared singular values of Phi (for vcov / diagnostics)
- `fit$nystrom_diagnostics$floored_count` — number of eigenvalues hit by the ridge floor
- `fit$nystrom_diagnostics$D_min_raw`, `D_max_raw` — pre-floor spectrum extremes

### `predict.krls_rr` — extended

Detects `fit$approx == "nystrom"` and builds the cross-kernel via the stored landmarks. Predictions are `K(X_new, Z) @ alpha` where `alpha` is the m-length dual coefficient. Same return shape: list with `$fit` (numeric vector) plus optional `$se.fit` etc.

### New helper export

`get_landmarks(fit, scale = c("original", "standardized"))` — accessor that returns the m×d landmark coordinates. Default unstandardizes via training centers/scales so users can pass them back into a new `krls(..., landmarks = ...)` call.

### Internal C++ exports — new

| Export | Purpose |
|---|---|
| `krls_nystrom_fit_cpp(X_tr_std, Z_std, y_tr, sigma, lambda_args, eps, compute_vcov)` | Single-fit Nystrom path. Returns alpha, lambda, fitted, optional vcov_alpha, plus W_eigen + Dinvsqrt + Sigma2 for the fit object |
| `krls_nystrom_autotune_inner_cpp(X_tr_std, Z_std, X_te_std, y_tr, y_te, sigma_grid, lambda_args, eps, nthreads)` | Autotune inner: per sigma builds C/W, eigen, lambda search, MSE on held-out. Returns `mse_per_sigma`, `lambda_per_sigma`, `nthreads_used` |
| `krls_nystrom_predict_cpp(X_new_std, Z_std, alpha, sigma)` | Predict via cross-kernel times alpha. Used by `predict.krls_rr` |

Existing C++ exports unchanged.

### No new top-level dependencies

`Suggests` gains nothing. `RhpcBLASctl` already discussed for Phase 1 may land here too if BLAS-thread coordination proves needed at large n.

---

## 4. Architecture

### Fit flow (single fit)

```
krls.default(X, y, approx="nystrom", ...):
  1. Standardize X -> Xs (existing path)
  2. Resolve landmarks via .resolve_landmarks() -> Z_std (m x d), method, indices
  3. Dispatch to:
       a. NON-AUTOTUNE: krls_nystrom_fit_cpp(Xs, Z_std, y, sigma, ...)
       b. AUTOTUNE: per (cross, fold) call krls_nystrom_autotune_inner_cpp,
          aggregate MSE matrix, 1-SE rule, refit at chosen sigma via
          krls_nystrom_fit_cpp on full Xs
  4. Build fit object with $approx="nystrom" + Nystrom bookkeeping
```

### Predict flow

```
predict.krls_rr(fit, newdata):
  if fit$approx == "exact":
    existing path  (unchanged)
  else:
    Xn_std = standardize(newdata, fit$X_means, fit$X_sds)
    yhat = krls_nystrom_predict_cpp(Xn_std, fit$landmarks, fit$coeffs, fit$sigma)
    return list(fit = yhat, ...)
```

### Per-fold landmarks in autotune

Per the Q5 decision: landmarks drawn once per (cross, fold) cell, shared across all sigmas within that fold (landmarks are X-coordinates, no sigma dependence). The C++ inner builds `C(sigma) = K(X_tr, Z, sigma)` and `W(sigma) = K(Z, Z, sigma)` per sigma elementwise from cached pairwise sqdist matrices `D_train_Z = ||X_tr - Z||²` and `D_Z_Z = ||Z - Z||²`, computed once per fold.

---

## 5. C++ internals

### `krls_nystrom_fit_cpp` (~200 LOC)

Single Nystrom fit at fixed sigma.

```cpp
Rcpp::List krls_nystrom_fit_cpp(const arma::mat& X_tr, const arma::mat& Z,
                                const arma::vec& y_tr, double sigma,
                                Rcpp::List lambda_args, double eps,
                                bool compute_vcov);
```

Steps:
1. Build `C = exp(-||X_tr - Z||² / sigma)` (n×m) via `krls_pairwise_sqdist_cpp` (Phase 1).
2. Build `W = exp(-||Z - Z||² / sigma)` (m×m).
3. `eig_sym(D, U, W)` → m-vector D, m×m U.
4. `D = max(D, 0)`; `Dmax = max(D)`. If `Dmax <= 0`, `Rcpp::stop("anchor W numerically zero")`.
5. `D_reg = max(D, eps * Dmax)`; `Dinvsqrt = 1 / sqrt(D_reg)`.
6. `Phi = C @ U @ diag(Dinvsqrt)` (n×m).
7. SVD of Phi: `(U_phi, sigma_phi, V_phi)`. Cache `Sigma2 = sigma_phi^2`.
8. Golden-section LOO lambda search on `f(lambda) = sum(((y - yfit) / (1 - diagS))^2)` where `w = Sigma2 / (Sigma2 + lambda)`, `yfit = U_phi @ (w * U_phi' y)`, `diagS = (U_phi^2) @ w`.
9. At chosen lambda: solve `beta`, recover `alpha = U @ (Dinvsqrt * beta)`.
10. Optional vcov: closed form via `Sigma2`-weighted outer products of `V_phi` and `alpha_map`.

Return: `Rcpp::List` with `coeffs` (m-vector), `fitted_std` (n-vector), `lambda`, `Looe`, `landmarks` (just echoes Z), `W_eigen$values`, `W_eigen$vectors`, `Dinvsqrt`, `Sigma2`, optional `vcov_alpha_std`, plus diagnostics `floored_count`, `D_min_raw`, `D_max_raw`.

### `krls_nystrom_autotune_inner_cpp` (~180 LOC)

```cpp
Rcpp::List krls_nystrom_autotune_inner_cpp(
  const arma::mat& X_tr, const arma::mat& Z, const arma::mat& X_te,
  const arma::vec& y_tr, const arma::vec& y_te,
  const arma::vec& sigma_grid, Rcpp::List lambda_args, double eps,
  int nthreads);
```

Mirrors Phase 1's `krls_autotune_inner_cpp` shape but with Nystrom internals.

Per-fold setup (caller's responsibility OR inside the function):
- `D_tr_Z = krls_pairwise_sqdist_cpp(X_tr, Z)` — n×m, computed once.
- `D_Z_Z = krls_pairwise_sqdist_cpp(Z, Z)` — m×m, computed once.
- `D_te_Z = krls_pairwise_sqdist_cpp(X_te, Z)` — n_te × m, computed once.

To save one allocation, the function can take the three pre-computed distance matrices instead of X_tr / X_te / Z. Final signature TBD by implementer — either pattern works; recommend passing the X matrices and computing distances inside for caller simplicity.

Per-sigma worker (TBB parallelFor with `numThreads=nthreads` per the Phase 1 lesson):
1. `C = exp(-D_tr_Z / s)` (n×m).
2. `W = exp(-D_Z_Z / s)` (m×m).
3. `K_te = exp(-D_te_Z / s)` (n_te × m).
4. `eig_sym(D, U, W)`; `D_reg = pmax(D, eps * max(D))`; `Dinvsqrt`.
5. `Phi = C @ U @ diag(Dinvsqrt)`.
6. `svd(Phi)` → `U_phi`, `Sigma2`.
7. Golden-section lambda search (same as fit path).
8. `alpha = U @ (Dinvsqrt * beta)`.
9. `yhat_te = K_te @ alpha`.
10. `mse_te(i) = mean((y_te - yhat_te)^2)`.

Each worker writes to a unique `mse(i)` and `lam(i)` slot → race-free → byte-identical across `nthreads`.

CRITICAL: use `RcppParallel::parallelFor(0, nsigma, worker, grainSize=1, numThreads=n_workers)` (Phase 1 lesson — `tbb::task_arena` wrap is bypassed by RcppParallel's internal arena).

### `krls_nystrom_predict_cpp` (~30 LOC)

```cpp
arma::vec krls_nystrom_predict_cpp(const arma::mat& X_new, const arma::mat& Z,
                                   const arma::vec& alpha, double sigma);
```

`K_new = exp(-||X_new - Z||² / sigma)` (n_new × m) then `yhat = K_new @ alpha`. Single export, trivial.

---

## 6. R harness changes

### `.resolve_landmarks` helper

R-side port of the reference helper (~80 LOC). Handles four cases: NULL (auto-draw via random or kmeans), integer indices, numeric matrix (standardized via training centers/scales), and validation. Uses `.with_seed()` wrapper to localize random draws.

### `krls.default` signature

Append five new optional arguments to the existing signature:
```r
krls.default <- function(X, y, ...,
                         approx = c("exact", "nystrom"),
                         nystrom_m = NULL,
                         landmarks = NULL,
                         landmark_method = c("random", "kmeans"),
                         landmark_seed = NULL,
                         nystrom_eps = 1e-9,
                         ...)
```

Validate via `match.arg`. Top-of-function gate:
```r
approx <- match.arg(approx)
landmark_method <- match.arg(landmark_method)
nystrom_eps <- .validate_nystrom_eps(nystrom_eps)
```

### Routing logic

After standardization, before kernel build:
```r
if (approx == "nystrom") {
  lm <- .resolve_landmarks(landmarks, landmark_method, nystrom_m,
                           Xs, X_centers, X_scales, landmark_seed)
  Z_std <- lm$matrix
  if (autotune) {
    # Autotune Nystrom path: per-fold landmarks, C++ inner
    # (Per Q5: re-resolve landmarks per (cross, fold) cell)
    # ...
  } else {
    nys <- krls_nystrom_fit_cpp(Xs, Z_std, y, sigma, lambda_args,
                                nystrom_eps, compute_vcov = vcov)
    # Assemble fit object with Nystrom bookkeeping
  }
} else {
  # Existing exact path - unchanged
}
```

### Autotune Nystrom — per-fold landmark resolution

Inside the autotune (cross, fold) loop, when approx == "nystrom":
```r
for (cross_i in seq_len(ncross)) {
  fold_seed <- as.integer(cross_seeds[cross_i])
  folds <- .ares_make_folds(n, nfold, fold_seed)

  for (k in seq_len(nfold)) {
    idx_te <- folds[[k]]; X_tr <- Xs[-idx_te, ]; X_te <- Xs[idx_te, ]
    y_tr <- y[-idx_te]; y_te <- y[idx_te]
    n_tr <- nrow(X_tr)

    # Per-fold landmark resolution
    nystrom_m_eff <- if (is.null(nystrom_m)) ceiling(sqrt(n_tr) * 2) else nystrom_m
    fold_landmark_seed <- if (is.null(landmark_seed)) NULL else landmark_seed + cross_i * 1000L + k
    lm <- .resolve_landmarks(landmarks, landmark_method, nystrom_m_eff,
                             X_tr, X_centers, X_scales, fold_landmark_seed)
    Z_std <- lm$matrix

    inner <- krls_nystrom_autotune_inner_cpp(
      X_tr, Z_std, X_te, y_tr, y_te,
      sigma_grid_sorted, lambda_args, nystrom_eps,
      nthreads = autotune.nthreads
    )

    mse_arr[cross_i, k, ] <- inner$mse_per_sigma
    lam_arr[cross_i, k, ] <- inner$lambda_per_sigma
  }
}
```

Notes:
- `fold_landmark_seed` derives from `landmark_seed + cross_i * 1000L + k` so each fold has a distinct deterministic seed (when `landmark_seed` is set) but the overall scheme is reproducible.
- When `landmarks` is user-supplied (indices or matrix), it's used as-is and the per-fold draw is skipped. Caveat: user-indices into the full training set may not be valid in the fold's training subset — the helper must validate and either remap or error. Implementation detail for the planner.

### `predict.krls_rr` extension

Add a branch at the top:
```r
if (object$approx == "nystrom") {
  Xn_std <- scale(newdata, center = object$X_means, scale = object$X_sds)
  yhat_std <- krls_nystrom_predict_cpp(Xn_std, object$landmarks,
                                       object$coeffs, object$sigma)
  yhat <- yhat_std * object$y_sd + object$y_mean  # de-standardize
  return(list(fit = yhat))
}
# Existing exact-path code unchanged
```

### `get_landmarks` exported helper

New ~20 LOC R function. Errors if `fit$approx != "nystrom"`. Default returns landmarks in original X scale by reversing the standardization.

---

## 7. Testing

### New file `tests/testthat/test-krls-nystrom.R` (~14 tests)

**Kernel + landmarks**

1. **NYS-XK-1** — `krls_nystrom_predict_cpp` cross-kernel matches `as.matrix(dist(rbind(Xnew, Z)))^2`-derived reference at sigma=5.
2. **LAND-RAND-1** — `.resolve_landmarks(NULL, "random", 50, Xs, ...)` returns m=50, no duplicate indices, all in 1:n.
3. **LAND-RAND-SEED-1** — Same `landmark_seed` produces identical indices across two calls; different seeds produce different indices.
4. **LAND-IDX-1** — Integer vector input returns matching X_std rows; out-of-range indices error.
5. **LAND-MAT-1** — Numeric matrix input is standardized correctly; ncol mismatch errors.

**Fit parity**

6. **NYS-FIT-1** — `krls_nystrom_fit_cpp` matches reference R impl (`/tmp/krls-reference/R/nystrom.R::.fit_krls_nystrom`) at fixed (sigma, lambda, landmarks) to `tol = 1e-8` on coefs and fitted.
7. **NYS-FIT-2** — Repeated calls with same `landmark_seed` produce byte-identical `fit$sigma`, `fit$lambda`, `fit$coeffs`.

**Autotune**

8. **NYS-AT-1** — `krls(..., autotune=TRUE, approx="nystrom", landmark_seed=42)` produces a finite `fit$sigma`, `fit$lambda` from the sigma grid.
9. **NYS-AT-EQUIV** — `autotune.nthreads=1` vs `autotune.nthreads=4` byte-identical with `landmark_seed=42` and `seed.cv=42`.
10. **NYS-AT-FALL** — `approx="nystrom"` + `weights` (or `L`/`U`/`tol` user-supplied) — must either fall back to R loop cleanly OR error with a clear message.

**Predict**

11. **NYS-PRED-1** — `predict(fit, newdata=X[1:10,])` returns length-10 numeric on a fitted Nystrom model; matches R-reference predict.

**Determinism**

12. **NYS-DET-1** — At fixed `landmark_seed`, two consecutive `krls(approx="nystrom")` calls produce byte-identical fits.
13. **NYS-DET-2** — Setting `landmark_seed` does NOT disturb the caller's `.Random.seed`.

**Accessor**

14. **NYS-GETLM-1** — `get_landmarks(fit, scale="original")` returns m×d matrix; `get_landmarks(fit, scale="standardized")` matches `fit$landmarks`. Errors on exact fit.

### Empirical sim — EMP-PHASE2

New file `inst/sims/krls-nystrom-phase2.R`. Mandatory session rules apply:
- Smoke test first (1 cell, R=2).
- Intermediate CSV per cell.
- Wall-clock cap 5 min.

**Cells**: 8 = {n=2000 p=10, n=5000 p=10} × {approx=exact, approx=nystrom} × {autotune=FALSE, autotune=TRUE}. R=5 reps per cell.

**Metrics per replicate**:
- wall_s, test_RMSE, sigma_hat, lambda_hat, nystrom_m_used, landmark_method_used

**Acceptance**:
- Single-fit speedup ≥ 5× at (n=2000, autotune=FALSE)
- Single-fit speedup ≥ 10× at (n=5000, autotune=FALSE)
- Autotune speedup ≥ 3× at n=2000
- Predict-time speedup ≥ 3× at n_test=1000
- Nystrom RMSE within +10% of exact at default m

### Regression

- Full krls suite passes (`devtools::test(filter="krls")` 0 fail)
- Phase 1 tests (`krls-autotune-parallel`) still pass — exact path unaffected
- R CMD check 0E/0W/no new notes

---

## 8. Rollout

### Version and commit shape

v0.0.0.9043. Multi-commit ship — too large for one coherent commit. Suggested decomposition (lines up with implementation plan tasks):

1. `krls_pairwise_sqdist_cpp` reuse (no new commit; Phase 1 already shipped it).
2. `.resolve_landmarks` R helper + LAND-* tests.
3. `krls_nystrom_predict_cpp` + NYS-XK-1 + helper export `get_landmarks` + NYS-GETLM-1.
4. `krls_nystrom_fit_cpp` + NYS-FIT-1, NYS-FIT-2.
5. R harness routing for non-autotune Nystrom + NYS-DET-1, NYS-DET-2, NYS-PRED-1.
6. `krls_nystrom_autotune_inner_cpp` + autotune harness + NYS-AT-1, NYS-AT-EQUIV, NYS-AT-FALL.
7. Roxygen / NEWS / DESCRIPTION / man regen / EMP-PHASE2 sim + speedup numbers / final commit.

### Pipeline

Full statsclaw workflow OR subagent-driven-development with per-task implementer + reviewer pair.

### Definition of done

- All 14 new tests + EMP-PHASE2 pass.
- Speedup targets met empirically (5× at n=2000 single fit; 10× at n=5000; 3× autotune).
- RMSE within +10% of exact at default m.
- Determinism: NYS-AT-EQUIV passes `expect_identical` across `autotune.nthreads`.
- v0.0.0.9043 entry in NEWS with measured speedup + RMSE deltas.
- Pushed to origin/main.

---

## 9. Risks

| Risk | Likelihood | Mitigation |
|---|---|---|
| C++ port has FP-order drift vs reference R impl | medium | NYS-FIT-1 parity test at tol=1e-8; if fails, document tolerance + investigate. The reference uses base `eigen()` (LAPACK `dsyevr`); Armadillo `eig_sym` uses the same LAPACK routine — should agree. |
| `nystrom_eps` floor too aggressive → silent loss of accuracy at small m | low | Diagnostics `floored_count` exposed in fit object; default 1e-9 matches reference. |
| Autotune Nystrom + per-fold landmark draw inflates CV variance | medium | NYS-AT-EQUIV pins determinism via `landmark_seed`; if CV variance is high in EMP-PHASE2, raise default ncross or document trade-off. |
| RMSE degrades > 10% on hard DGPs (interaction, exp-decay) at default m | medium | Bump default m → `ceiling(sqrt(n)*3)` if EMP-PHASE2 reveals this. Easy knob. |
| `landmarks` as user-supplied indices breaks in autotune fold subsets | high | Per-fold path must either: (a) treat user indices as global-X-row indices and skip the fold's training-subset check (with documented caveat), OR (b) error if user supplies indices with autotune. Implementer chooses safer route. |
| Predict path with stored landmarks fails on factor-expanded newdata | low | Existing `predict.krls_rr` already handles factor expansion via `model.matrix`; Nystrom branch reuses this code before the `krls_nystrom_predict_cpp` call. |
| `vcov` support adds significant LOC; deferring saves ~80 LOC | medium | First Phase-2 ship: `vcov = NULL` for Nystrom paths; document as "exact only for v0.0.0.9043". Add in v0.0.0.9044 if requested. |
| Memory: `C` is n×m × 8B (e.g. n=5000, m=200 → 8 MB); `Phi` similar | low | Within bounds for target n. |

---

## 10. Out of scope (Phase 3 candidates)

- Leverage-score / RLS landmark selection
- Auto-engage Nystrom at large n via heuristic
- Sparse kernel paths (truncate small `C` entries)
- Multi-response Nystrom
- Online landmark updates
- Bagging + Nystrom integration
- vcov for Nystrom fits

---

## Approval

Approved during brainstorming session 2026-05-19. Next step: invoke `superpowers:writing-plans` to produce the implementation plan.
