# pLDA — Penalized Linear Discriminant Analysis · Design Spec

**Date:** 2026-05-20
**Package:** roadrunner
**Status:** Approved (brainstorming) — pending implementation plan

## 1. Summary

Add `plda()`, a fourth fitter to roadrunner: penalized linear discriminant analysis implementing Witten & Tibshirani (2011), "Penalized classification using Fisher's linear discriminant," *JRSS-B* 73(5), 753–772. Supports L1 (lasso) and fused-lasso penalties, native multi-class, built-in k-fold CV autotune. Same architecture, restrictions, and objectives as `ares()`, `krls()`, `meep()`.

## 2. Algorithm

Reference: Witten, D. M. & Tibshirani, R. (2011). Penalized classification using Fisher's linear discriminant. *Journal of the Royal Statistical Society B*, 73(5), 753–772.

Given n×p feature matrix X, class labels y ∈ {1,…,G}, n_g per class:

1. **Within-class standardization.** Center each feature by its global mean; scale by within-class standard deviation σ̂_w,j. The method assumes a **diagonal** within-class covariance Σ̂_w = diag(σ̂²_w,j) — this keeps it well-posed for p ≫ n.

2. **Penalized Fisher criterion.** For discriminant vector β_k (k = 1,…,K):
   maximize_β  βᵀ Σ̂_b β − P(β)   subject to  βᵀ Σ̂_w β ≤ 1
   where Σ̂_b is the between-class covariance.
   - **L1:** P(β) = λ Σ_j |β_j|
   - **Fused:** P(β) = λ Σ_j |β_j| + λ₂ Σ_j |β_j − β_{j−1}|

3. **Minorize-maximize (MM).** The between-class term βᵀΣ̂_b β is convex; minorize it by its linearization at the current iterate. Each MM step reduces to a penalized least-squares problem in a working response:
   - L1 → coordinate soft-thresholding (closed form per coordinate).
   - Fused → fused-lasso signal approximator (FLSA), solved by the Johnson (2013) / Condat (2013) 1-D total-variation dynamic program.
   Iterate to `tol` or `maxit`.

4. **Multi-class via deflation.** After β_k is found, deflate: project Σ̂_b onto the orthogonal complement (in the Σ̂_w metric) of the discriminant vectors found so far, repeat. Up to K ≤ G−1 vectors.

5. **Prediction.** Project X onto the K discriminant vectors → n×K scores. Project class centroids likewise. Assign each observation to the nearest centroid in projected space (Euclidean = Mahalanobis in the diagonal-Σ_w metric). Posterior = softmax over negative half-squared projected distances.

## 3. Architecture

Two-file engine layout, identical to `krls` and `ares`.

### R/plda.R
- `plda()` — S3 generic
- `plda.formula(formula, data, ...)` — formula interface; builds model matrix, extracts factor response
- `plda.default(x, y, ...)` — matrix interface; core entry
- `.plda_cv()` — k-fold CV autotune over (λ, K) grid
- `.plda_standardize()` — within-class moment helper (R-side glue; heavy compute in C++)
- input validation, class S3 construction

### src/plda.cpp
C++17 + RcppArmadillo + RcppParallel. Rcpp exports:
- `plda_fit_cpp` — full fit for a fixed (λ, λ₂, K, penalty): within-class moments, MM discriminant solver, deflation loop. Returns discriminant matrix (p×K), projected centroids, per-class means, σ̂_w, iteration counts.
- `plda_project_cpp` — project a new X onto stored discriminant vectors; used by `predict`.
- Internal (not exported): soft-threshold sweep, FLSA dynamic program, MM inner loop.

### R/predict.R
- `predict.plda(object, newdata, type = c("class","posterior","projection"))`. OOV factor-level error message consistent with `predict.ares`.

### R/print.R
- `print.plda` — call, classes, penalty, λ/K, #nonzero features per discriminant
- `summary.plda` — per-class counts, discriminant sparsity, CV-selected λ/K, CV error
- `plot.plda` — scatter of first two discriminant projections coloured by class; K=1 → per-class density strip

### Supporting
- `tests/testthat/test-plda-*.R` — see §6
- `_pkgdown.yml` — new "Penalized LDA: plda()" reference section (avoid the topic-index CI break)
- `DESCRIPTION` — add `penalizedLDA` to **Suggests**; bump version
- `NEWS.md`, `README.md`, `CLAUDE.md` — document the new fitter
- `man/` — roxygen-generated

## 4. API

```r
plda(x, y,
     K           = NULL,                 # # discriminant vectors; default min(G-1, autotune)
     lambda      = NULL,                 # L1 / fused magnitude penalty
     penalty     = c("L1", "fused"),
     lambda2     = NULL,                 # fused difference penalty (penalty="fused")
     autotune    = TRUE,
     nfold       = 5L,
     lambda_grid = NULL,                 # default: data-driven log-spaced grid
     maxit       = 100L,
     tol         = 1e-6,
     nthreads    = getOption("roadrunner.threads"))
# formula interface
plda(species ~ ., data = iris)
```

- S3 class `"plda"`.
- `autotune = TRUE` (default): k-fold CV over `lambda_grid` × candidate K by misclassification error; refit at optimum. Disabled when `lambda` (and K) supplied explicitly and `autotune = FALSE`.
- `penalty = "fused"` treats column order as the natural feature ordering.
- Determinism: `nthreads = 1 ≡ nthreads = N` byte-identical. CV folds seeded.

## 5. Parallelism

Target the same parallel model as `ares`/`krls` (TBB via RcppParallel):
- **CV fold loop** — embarrassingly parallel; each fold writes to a fixed result slot.
- **Per-feature sweeps** — within-class moments and the L1 soft-threshold sweep parallelise across columns.
- The MM iteration itself is sequential (each step depends on the previous).

Determinism: parallel regions only over independent indices with fixed-slot writes → byte-identical across thread counts.

**Milestone fallback:** if the parallel engine risks slipping schedule, a **sequential-correct** engine ships first (Milestone 1), parallelism added in Milestone 2. Correctness and determinism are never traded; only the speed objective is staged.

## 6. Testing

- **Parity** — numeric agreement with `penalizedLDA::PenalizedLDA` / `PenalizedLDA.cv` on fixed seeds; Suggests-gated, `skip_if_not_installed("penalizedLDA")`.
- **Determinism** — `nthreads = 1` vs `N` byte-identical, every path.
- **Multi-class** — G = 2 and G > 2; K up to G−1.
- **Penalties** — L1 and fused both; fused recovers piecewise-constant discriminants on ordered-feature simulation.
- **Edge** — p ≫ n; λ large → all-zero discriminant; single informative feature; constant feature dropped; singleton class handled or rejected with clear message; K > G−1 rejected.
- **CV** — autotune recovers known λ on simulated data with planted signal.
- **Interfaces** — formula ≡ matrix; OOV factor level in `predict` errors clearly.
- Full suite + `lintr::lint_package()` clean + all 4 CI workflows green before merge.

## 7. Scope boundary (YAGNI)

**In v1:** L1 + fused penalties, G ≥ 2 classes, dense numeric X, built-in CV autotune, S3 predict/print/summary/plot, formula + matrix interfaces.

**Punted (explicit non-goals):** sparse-matrix X input; non-diagonal Σ_w; missing-data handling; bagging (discriminant directions do not bag like trees); penalties beyond L1/fused.

## 8. Invariants honored

1. Determinism: nthreads-invariant, byte-identical.
2. Hard deps unchanged — Rcpp, RcppParallel, RcppArmadillo. `penalizedLDA` is Suggests-only (parity tests skip if absent).
3. C++17, cross-platform, no tidyverse/rlang/S4/caret.
4. Two-file engine layout consistent with `krls`/`ares`.
5. pkgdown reference index updated in the same change.
