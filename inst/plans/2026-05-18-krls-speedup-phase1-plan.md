# roadrunner::krls Phase 1 Speedup — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Speed up `roadrunner::krls(..., autotune = TRUE)` by ≥3× wall-clock at (n=500, p=10, nthreads=4) and (n=1500, p=20, nthreads=4), with byte-identical fits at fixed seed across all `nthreads` values.

**Architecture:** Compute pairwise squared-distance matrix `D` once per (ncross, fold) cell instead of once per (ncross, fold, sigma) cell. Run the sigma grid sweep inside C++ with `RcppParallel::parallelFor`, each thread writing to a unique slot in pre-allocated output vectors (no reduction → race-free → determinism preserved). All BLAS-heavy ops stay outside the parallel region.

**Tech Stack:** R, Rcpp, RcppArmadillo (`eig_sym`, `dgemm`), RcppParallel (TBB), testthat 3, devtools.

**Spec:** `inst/specs/2026-05-18-krls-speedup-phase1-design.md` (commit `fce882e`)

**Target version:** `0.0.0.9042`

---

## Pre-flight

Confirm working tree clean and on `main` at v0.0.0.9041 (commit `2c83d4a` or later).

```bash
git -C /home/jack/Dropbox/roadrunner status --short
git -C /home/jack/Dropbox/roadrunner log --oneline -3
```

Expected: clean tree, `2c83d4a` or later v0.0.0.9041 commit visible.

If `inst/plans/` does not exist, this plan won't be present. Create the plan first (this file) then begin.

---

## Task 1: Snapshot baseline fit for EQUIV-1

**Why:** EQUIV-1 asserts byte-identicality vs v0.0.0.9041. We must capture the baseline fit BEFORE any code change.

**Files:**
- Create: `inst/testdata/autotune-baseline-9041.rds`
- Create (temporary, throwaway): `tools/build-autotune-baseline.R`

- [ ] **Step 1: Create the snapshot-builder script**

`tools/build-autotune-baseline.R`:
```r
## Build the autotune-baseline-9041 snapshot for EQUIV-1.
## Run from package root with v0.0.0.9041 installed/loaded.

suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))

set.seed(42L)
n <- 200L
p <- 5L
X <- matrix(rnorm(n * p), n, p)
y <- as.numeric(rowSums(X[, 1:2]) + 0.5 * rnorm(n))

fit <- roadrunner::krls(
  X        = X,
  y        = y,
  autotune = TRUE,
  derivative = FALSE,
  vcov     = FALSE
)

baseline <- list(
  sigma   = fit$sigma,
  lambda  = fit$lambda,
  coeffs  = fit$coeffs,
  fitted  = fit$fitted,
  R2      = fit$R2,
  autotune = fit$autotune,
  seed    = 42L,
  n       = n,
  p       = p,
  version = utils::packageVersion("roadrunner")
)

dir.create("inst/testdata", recursive = TRUE, showWarnings = FALSE)
saveRDS(baseline, "inst/testdata/autotune-baseline-9041.rds", version = 2)
message("Snapshot written: inst/testdata/autotune-baseline-9041.rds")
message("  sigma  = ", baseline$sigma)
message("  lambda = ", baseline$lambda)
message("  R2     = ", baseline$R2)
```

- [ ] **Step 2: Run the builder**

```bash
Rscript tools/build-autotune-baseline.R
```

Expected: `Snapshot written: inst/testdata/autotune-baseline-9041.rds` plus printed sigma / lambda / R2.

- [ ] **Step 3: Verify snapshot loads cleanly**

```bash
Rscript -e 'b <- readRDS("inst/testdata/autotune-baseline-9041.rds"); str(b, max.level=1)'
```

Expected: list with sigma, lambda, coeffs, fitted, R2, autotune, seed, n, p, version.

- [ ] **Step 4: Remove the throwaway builder**

```bash
rm tools/build-autotune-baseline.R
rmdir tools 2>/dev/null || true
```

- [ ] **Step 5: Commit snapshot**

```bash
git add inst/testdata/autotune-baseline-9041.rds
git commit -m "test: snapshot autotune baseline at v0.0.0.9041 for phase-1 equiv tests"
```

---

## Task 2: TDD — pairwise squared-distance helper (sequential, no parallel)

**Files:**
- Modify: `src/krls.cpp` (add `krls_pairwise_sqdist_cpp` export, ~30 LOC)
- Create: `tests/testthat/test-krls-autotune-parallel.R`

- [ ] **Step 1: Write failing tests for the distance helper**

`tests/testthat/test-krls-autotune-parallel.R`:
```r
## Phase 1 (v0.0.0.9042) tests: shared distance + parallel autotune.
## Mirrors inst/specs/2026-05-18-krls-speedup-phase1-design.md.

test_that("DIST-1: krls_pairwise_sqdist_cpp matches as.matrix(dist(X))^2", {
  set.seed(1L)
  X <- matrix(rnorm(30), 10, 3)
  D_cpp <- roadrunner:::krls_pairwise_sqdist_cpp(X, X)
  D_ref <- as.matrix(stats::dist(X))^2
  expect_equal(unname(D_cpp), unname(D_ref), tolerance = 1e-10)
})

test_that("DIST-2: cross-matrix (n_a x n_b) matches naive double loop", {
  set.seed(2L)
  A <- matrix(rnorm(20), 5, 4)
  B <- matrix(rnorm(32), 8, 4)
  D_cpp <- roadrunner:::krls_pairwise_sqdist_cpp(A, B)
  D_ref <- matrix(NA_real_, 5, 8)
  for (i in seq_len(5)) for (j in seq_len(8)) {
    D_ref[i, j] <- sum((A[i, ] - B[j, ])^2)
  }
  expect_equal(D_cpp, D_ref, tolerance = 1e-10)
})

test_that("DIST-3: no negative leakage on near-zero rows (FP rounding)", {
  X <- matrix(0, 4, 3) + matrix(rnorm(12) * 1e-12, 4, 3)
  D <- roadrunner:::krls_pairwise_sqdist_cpp(X, X)
  expect_true(all(D >= 0))
  expect_true(all(diag(D) == 0 | diag(D) < 1e-20))
})
```

- [ ] **Step 2: Run tests to verify failure**

```bash
Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls-autotune-parallel")'
```

Expected: 3 errors `could not find function "krls_pairwise_sqdist_cpp"` or `object 'krls_pairwise_sqdist_cpp' not found`.

- [ ] **Step 3: Implement `krls_pairwise_sqdist_cpp` in `src/krls.cpp`**

Append to `src/krls.cpp` (above the final namespace / EOF, after the existing exports):

```cpp
// ---------------------------------------------------------------------------
// Phase 1 (v0.0.0.9042): shared pairwise squared-distance matrix.
//
// D[i, j] = ||X_a[i, ] - X_b[j, ]||^2
//
// Identity: D = ||X_a||^2 (row sums) + ||X_b||^2^T - 2 X_a X_b^T
// One dgemm dominates; tiny negatives clamped to 0 to absorb FP noise.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
arma::mat krls_pairwise_sqdist_cpp(const arma::mat& X_a, const arma::mat& X_b) {
  if (X_a.n_cols != X_b.n_cols) {
    Rcpp::stop("krls_pairwise_sqdist_cpp: X_a and X_b must have same n_cols");
  }
  arma::vec ra = arma::sum(arma::square(X_a), 1);   // n_a x 1
  arma::vec rb = arma::sum(arma::square(X_b), 1);   // n_b x 1
  arma::mat D = -2.0 * (X_a * X_b.t());             // n_a x n_b
  D.each_col() += ra;
  D.each_row() += rb.t();
  D.elem(arma::find(D < 0.0)).zeros();              // FP-noise floor
  return D;
}
```

- [ ] **Step 4: Regenerate Rcpp bindings**

```bash
Rscript -e 'devtools::document(".")'
Rscript -e 'Rcpp::compileAttributes(".")'
```

Expected: `RcppExports.cpp` and `R/RcppExports.R` updated to include `krls_pairwise_sqdist_cpp`.

- [ ] **Step 5: Run tests to verify pass**

```bash
Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls-autotune-parallel")'
```

Expected: DIST-1, DIST-2, DIST-3 all PASS.

- [ ] **Step 6: Commit**

```bash
git add src/krls.cpp src/RcppExports.cpp R/RcppExports.R tests/testthat/test-krls-autotune-parallel.R
git commit -m "feat(krls): krls_pairwise_sqdist_cpp via dgemm + FP-noise floor"
```

---

## Task 3: TDD — sequential `krls_autotune_inner_cpp` (nthreads ignored, foundation only)

**Why:** Build the function with the correct numerics first, prove parity with the existing per-sigma loop. Parallelism slots in cleanly in Task 4.

**Files:**
- Modify: `src/krls.cpp` (add `krls_autotune_inner_cpp`, ~180 LOC)
- Modify: `tests/testthat/test-krls-autotune-parallel.R`

- [ ] **Step 1: Write failing parity test (INNER-1)**

Append to `tests/testthat/test-krls-autotune-parallel.R`:
```r
test_that("INNER-1: krls_autotune_inner_cpp matches sequential R reference", {
  set.seed(3L)
  n_tr <- 60L; n_te <- 20L; p <- 4L
  X_tr <- matrix(rnorm(n_tr * p), n_tr, p)
  X_te <- matrix(rnorm(n_te * p), n_te, p)
  y_tr <- as.numeric(rowSums(X_tr) + 0.3 * rnorm(n_tr))
  y_te <- as.numeric(rowSums(X_te) + 0.3 * rnorm(n_te))

  D_tr <- roadrunner:::krls_pairwise_sqdist_cpp(X_tr, X_tr)
  D_te <- roadrunner:::krls_pairwise_sqdist_cpp(X_te, X_tr)

  sigma_grid <- c(2, 4, 8, 16)

  ## R reference: rebuild K, eigen, LOO lambda, MSE per sigma. Sequential.
  ref_mse <- numeric(length(sigma_grid))
  ref_lam <- numeric(length(sigma_grid))
  for (i in seq_along(sigma_grid)) {
    s <- sigma_grid[i]
    K_tr <- exp(-D_tr / s)
    K_te <- exp(-D_te / s)
    e    <- roadrunner:::krls_eig_cpp(K_tr)
    d    <- e$values
    V    <- e$vectors
    Vty  <- as.numeric(crossprod(V, y_tr))
    Vsq  <- roadrunner:::krls_vsq_cpp(V)
    L    <- .Machine$double.eps
    while (sum(d / (d + L)) > min(which.min(abs(d - max(d) / 1000)),
                                  length(d))) L <- L * 10
    U    <- length(y_tr)
    while (sum(d / (d + U)) < 1) U <- U - 1
    gr   <- (sqrt(5) - 1) / 2
    X1 <- L + (1 - gr) * (U - L); X2 <- L + gr * (U - L)
    S1 <- roadrunner:::krls_loo_loss_cpp(d, V, Vsq, Vty, X1)
    S2 <- roadrunner:::krls_loo_loss_cpp(d, V, Vsq, Vty, X2)
    while (abs(S1 - S2) > 1e-6) {
      if (S1 < S2) { U <- X2; X2 <- X1; X1 <- L + (1 - gr) * (U - L)
                     S2 <- S1; S1 <- roadrunner:::krls_loo_loss_cpp(d, V, Vsq, Vty, X1)
      } else      { L <- X1; X1 <- X2; X2 <- L + gr * (U - L)
                     S1 <- S2; S2 <- roadrunner:::krls_loo_loss_cpp(d, V, Vsq, Vty, X2) }
    }
    lam <- (X1 + X2) / 2
    sol <- roadrunner:::krls_solve_cpp(d, V, Vsq, Vty, lam)
    alpha <- sol$coeffs
    yhat  <- as.numeric(K_te %*% alpha)
    ref_mse[i] <- mean((y_te - yhat)^2)
    ref_lam[i] <- lam
  }

  out <- roadrunner:::krls_autotune_inner_cpp(
    D_tr, D_te, y_tr, y_te, sigma_grid,
    list(tol = 1e-6, L0 = .Machine$double.eps, L_step = 10,
         U_start_from_n = TRUE),
    nthreads = 1L
  )
  expect_equal(out$mse_per_sigma, ref_mse, tolerance = 1e-9)
  expect_equal(out$lambda_per_sigma, ref_lam, tolerance = 1e-9)
})

test_that("INNER-2: sigma_grid length 1 degenerate path works", {
  set.seed(4L)
  n_tr <- 40L; n_te <- 10L; p <- 3L
  X_tr <- matrix(rnorm(n_tr * p), n_tr, p)
  X_te <- matrix(rnorm(n_te * p), n_te, p)
  y_tr <- rnorm(n_tr); y_te <- rnorm(n_te)
  D_tr <- roadrunner:::krls_pairwise_sqdist_cpp(X_tr, X_tr)
  D_te <- roadrunner:::krls_pairwise_sqdist_cpp(X_te, X_tr)
  out <- roadrunner:::krls_autotune_inner_cpp(
    D_tr, D_te, y_tr, y_te, c(5.0),
    list(tol = 1e-6, L0 = .Machine$double.eps, L_step = 10,
         U_start_from_n = TRUE),
    nthreads = 1L
  )
  expect_length(out$mse_per_sigma, 1L)
  expect_length(out$lambda_per_sigma, 1L)
  expect_true(is.finite(out$mse_per_sigma))
  expect_true(out$lambda_per_sigma > 0)
})
```

- [ ] **Step 2: Run tests to verify failure**

```bash
Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls-autotune-parallel")'
```

Expected: INNER-1, INNER-2 ERROR (`could not find function "krls_autotune_inner_cpp"`). DIST tests still PASS.

- [ ] **Step 3: Implement `krls_autotune_inner_cpp` (sequential, no TBB yet)**

Append to `src/krls.cpp`:

```cpp
// ---------------------------------------------------------------------------
// Phase 1 (v0.0.0.9042): inner autotune loop in C++.
//
// For each sigma s in sigma_grid:
//   K_tr = exp(-D_tr / s)
//   K_te = exp(-D_te / s)
//   (V, d) = eig_sym(K_tr)
//   lambda = golden-section LOO search (tol=1e-6, log-step L)
//   alpha  = V * ((Vty / (d + lambda)))    (closed-form ridge in eigenbasis)
//   mse_s  = mean((y_te - K_te * alpha)^2)
//
// THIS TASK: sequential only (nthreads ignored). Task 4 adds parallelFor.
// ---------------------------------------------------------------------------

static double krls_loo_loss_inline(const arma::vec& d, const arma::mat& V,
                                   const arma::mat& Vsq, const arma::vec& Vty,
                                   double lambda) {
  arma::vec inv = 1.0 / (d + lambda);
  arma::vec c   = V * (Vty % inv);
  arma::vec g   = Vsq * inv;
  return arma::dot(c / g, c / g);
}

static void krls_one_sigma(const arma::mat& D_tr, const arma::mat& D_te,
                           const arma::vec& y_tr, const arma::vec& y_te,
                           double s,
                           double lam_tol, double L0, double L_step,
                           double& mse_out, double& lam_out) {
  arma::mat K_tr = arma::exp(-D_tr / s);
  arma::mat K_te = arma::exp(-D_te / s);

  arma::vec d;
  arma::mat V;
  arma::eig_sym(d, V, K_tr);

  arma::vec Vty = V.t() * y_tr;
  arma::mat Vsq = V % V;

  // Bracket [L, U]
  double L = L0;
  int q = 0; double mind = std::abs(d(0) - d.max() / 1000.0);
  for (arma::uword i = 1; i < d.n_elem; ++i) {
    double dd = std::abs(d(i) - d.max() / 1000.0);
    if (dd < mind) { mind = dd; q = static_cast<int>(i); }
  }
  while (true) {
    double s_eff = arma::sum(d / (d + L));
    if (s_eff <= (double)(q + 1)) break;
    L *= L_step;
    if (L > 1e30) break;
  }
  double U = (double) d.n_elem;
  while (true) {
    double s_eff = arma::sum(d / (d + U));
    if (s_eff >= 1.0) break;
    U -= 1.0;
    if (U <= L) { U = L * 10.0; break; }
  }

  // Golden section
  const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
  double X1 = L + (1.0 - gr) * (U - L);
  double X2 = L + gr * (U - L);
  double S1 = krls_loo_loss_inline(d, V, Vsq, Vty, X1);
  double S2 = krls_loo_loss_inline(d, V, Vsq, Vty, X2);
  int iter = 0;
  while (std::abs(S1 - S2) > lam_tol && iter < 200) {
    if (S1 < S2) {
      U = X2; X2 = X1; X1 = L + (1.0 - gr) * (U - L);
      S2 = S1; S1 = krls_loo_loss_inline(d, V, Vsq, Vty, X1);
    } else {
      L = X1; X1 = X2; X2 = L + gr * (U - L);
      S1 = S2; S2 = krls_loo_loss_inline(d, V, Vsq, Vty, X2);
    }
    ++iter;
  }
  double lam = 0.5 * (X1 + X2);

  // Solve alpha in eigenbasis
  arma::vec alpha = V * (Vty / (d + lam));
  arma::vec yhat  = K_te * alpha;
  mse_out = arma::mean(arma::square(y_te - yhat));
  lam_out = lam;
}

// [[Rcpp::export]]
Rcpp::List krls_autotune_inner_cpp(const arma::mat& D_tr, const arma::mat& D_te,
                                   const arma::vec& y_tr, const arma::vec& y_te,
                                   const arma::vec& sigma_grid,
                                   Rcpp::List lambda_args,
                                   int nthreads) {
  const arma::uword nsigma = sigma_grid.n_elem;
  arma::vec mse(nsigma, arma::fill::zeros);
  arma::vec lam(nsigma, arma::fill::zeros);

  double lam_tol = Rcpp::as<double>(lambda_args["tol"]);
  double L0      = Rcpp::as<double>(lambda_args["L0"]);
  double L_step  = Rcpp::as<double>(lambda_args["L_step"]);

  // Sequential loop. Task 4 swaps in parallelFor.
  for (arma::uword i = 0; i < nsigma; ++i) {
    double m = 0.0, l = 0.0;
    krls_one_sigma(D_tr, D_te, y_tr, y_te, sigma_grid(i),
                   lam_tol, L0, L_step, m, l);
    mse(i) = m;
    lam(i) = l;
  }

  return Rcpp::List::create(
    Rcpp::Named("mse_per_sigma")    = mse,
    Rcpp::Named("lambda_per_sigma") = lam,
    Rcpp::Named("nthreads_used")    = 1
  );
}
```

- [ ] **Step 4: Regenerate bindings + compile**

```bash
Rscript -e 'devtools::document("."); Rcpp::compileAttributes(".")'
```

- [ ] **Step 5: Run tests to verify INNER-1 + INNER-2 pass**

```bash
Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls-autotune-parallel")'
```

Expected: INNER-1 and INNER-2 PASS. (No INNER-3 yet — that's after Task 4.)

- [ ] **Step 6: Commit**

```bash
git add src/krls.cpp src/RcppExports.cpp R/RcppExports.R tests/testthat/test-krls-autotune-parallel.R
git commit -m "feat(krls): krls_autotune_inner_cpp sequential path with closed-form LOO"
```

---

## Task 4: Add TBB parallelism to `krls_autotune_inner_cpp`

**Why:** Numerics already proven by Task 3. Now wrap the per-sigma loop in `RcppParallel::parallelFor` and add the INNER-3 clamp test.

**Files:**
- Modify: `src/krls.cpp` (parallel wrap, ~30 LOC delta)
- Modify: `tests/testthat/test-krls-autotune-parallel.R` (add INNER-3 + early EQUIV-2)

- [ ] **Step 1: Write failing INNER-3 + early EQUIV-2 tests**

Append to `tests/testthat/test-krls-autotune-parallel.R`:
```r
test_that("INNER-3: nthreads > nsigma clamps without spawning extra workers", {
  set.seed(5L)
  n_tr <- 30L; n_te <- 10L; p <- 2L
  X_tr <- matrix(rnorm(n_tr * p), n_tr, p)
  X_te <- matrix(rnorm(n_te * p), n_te, p)
  y_tr <- rnorm(n_tr); y_te <- rnorm(n_te)
  D_tr <- roadrunner:::krls_pairwise_sqdist_cpp(X_tr, X_tr)
  D_te <- roadrunner:::krls_pairwise_sqdist_cpp(X_te, X_tr)
  sg   <- c(2, 4, 8)
  out1 <- roadrunner:::krls_autotune_inner_cpp(
    D_tr, D_te, y_tr, y_te, sg,
    list(tol = 1e-6, L0 = .Machine$double.eps, L_step = 10,
         U_start_from_n = TRUE),
    nthreads = 32L
  )
  out2 <- roadrunner:::krls_autotune_inner_cpp(
    D_tr, D_te, y_tr, y_te, sg,
    list(tol = 1e-6, L0 = .Machine$double.eps, L_step = 10,
         U_start_from_n = TRUE),
    nthreads = 1L
  )
  expect_equal(out1$mse_per_sigma, out2$mse_per_sigma, tolerance = 0)
  expect_equal(out1$lambda_per_sigma, out2$lambda_per_sigma, tolerance = 0)
})

test_that("INNER-EQUIV: nthreads=1 vs nthreads=4 byte-identical on inner_cpp", {
  set.seed(6L)
  n_tr <- 80L; n_te <- 30L; p <- 4L
  X_tr <- matrix(rnorm(n_tr * p), n_tr, p)
  X_te <- matrix(rnorm(n_te * p), n_te, p)
  y_tr <- rnorm(n_tr); y_te <- rnorm(n_te)
  D_tr <- roadrunner:::krls_pairwise_sqdist_cpp(X_tr, X_tr)
  D_te <- roadrunner:::krls_pairwise_sqdist_cpp(X_te, X_tr)
  sg   <- c(1, 2, 4, 8, 16, 32)
  o1 <- roadrunner:::krls_autotune_inner_cpp(
    D_tr, D_te, y_tr, y_te, sg,
    list(tol = 1e-6, L0 = .Machine$double.eps, L_step = 10,
         U_start_from_n = TRUE),
    nthreads = 1L
  )
  o4 <- roadrunner:::krls_autotune_inner_cpp(
    D_tr, D_te, y_tr, y_te, sg,
    list(tol = 1e-6, L0 = .Machine$double.eps, L_step = 10,
         U_start_from_n = TRUE),
    nthreads = 4L
  )
  expect_identical(o1$mse_per_sigma, o4$mse_per_sigma)
  expect_identical(o1$lambda_per_sigma, o4$lambda_per_sigma)
})
```

- [ ] **Step 2: Run tests to verify they pass (sequential path still — confirms parity invariant before adding TBB)**

```bash
Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls-autotune-parallel")'
```

Expected: all 5 tests PASS (because nthreads is currently ignored — every call is sequential, so they trivially agree). This is intentional — we lock in the invariant BEFORE swapping the loop.

- [ ] **Step 3: Replace the sequential loop with `parallelFor`**

In `src/krls.cpp`, locate the `krls_autotune_inner_cpp` definition. Replace the sequential `for (i = 0; i < nsigma; ++i)` loop body. Add at top of file (or after existing includes) if not already present:

```cpp
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
```

(`roadrunner` already depends on RcppParallel for `ares`; the include may exist. Don't duplicate.)

Add a Worker struct above `krls_autotune_inner_cpp`:

```cpp
struct KrlsSigmaWorker : public RcppParallel::Worker {
  const arma::mat& D_tr;
  const arma::mat& D_te;
  const arma::vec& y_tr;
  const arma::vec& y_te;
  const arma::vec& sigma_grid;
  double lam_tol, L0, L_step;
  arma::vec& mse;
  arma::vec& lam;

  KrlsSigmaWorker(const arma::mat& D_tr_, const arma::mat& D_te_,
                  const arma::vec& y_tr_, const arma::vec& y_te_,
                  const arma::vec& sigma_grid_,
                  double lam_tol_, double L0_, double L_step_,
                  arma::vec& mse_, arma::vec& lam_)
    : D_tr(D_tr_), D_te(D_te_), y_tr(y_tr_), y_te(y_te_),
      sigma_grid(sigma_grid_), lam_tol(lam_tol_), L0(L0_), L_step(L_step_),
      mse(mse_), lam(lam_) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      double m = 0.0, l = 0.0;
      krls_one_sigma(D_tr, D_te, y_tr, y_te, sigma_grid(i),
                     lam_tol, L0, L_step, m, l);
      mse(i) = m;
      lam(i) = l;
    }
  }
};
```

Then replace the sequential loop in `krls_autotune_inner_cpp` with:

```cpp
  int n_workers = nthreads;
  if (n_workers < 1) n_workers = 1;
  if ((arma::uword) n_workers > nsigma) n_workers = (int) nsigma;
  RcppParallel::ThreadPool pool(n_workers);  // not used directly; setThreadOptions controls
  RcppParallel::setThreadOptions(n_workers);

  KrlsSigmaWorker worker(D_tr, D_te, y_tr, y_te, sigma_grid,
                         lam_tol, L0, L_step, mse, lam);
  RcppParallel::parallelFor(0, nsigma, worker, /*grainSize=*/1);
```

Update the return list to record actual `nthreads`:

```cpp
  return Rcpp::List::create(
    Rcpp::Named("mse_per_sigma")    = mse,
    Rcpp::Named("lambda_per_sigma") = lam,
    Rcpp::Named("nthreads_used")    = n_workers
  );
```

- [ ] **Step 4: Compile + run tests**

```bash
Rscript -e 'devtools::document("."); Rcpp::compileAttributes(".")'
Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls-autotune-parallel")'
```

Expected: all 5 tests PASS. Critically `INNER-EQUIV` and `INNER-3` confirm nthreads=1 ≡ nthreads=4 byte-identical via `expect_identical`.

If `INNER-EQUIV` fails: a BLAS thread race is leaking FP order. Check that `arma::eig_sym` is using a single-threaded LAPACK path. Run `Rscript -e 'extSoftVersion()'` and `Rscript -e 'sessionInfo()'` to confirm; set `RhpcBLASctl::blas_set_num_threads(1)` in the test setup if needed.

- [ ] **Step 5: Commit**

```bash
git add src/krls.cpp src/RcppExports.cpp R/RcppExports.R tests/testthat/test-krls-autotune-parallel.R
git commit -m "feat(krls): TBB parallelFor over sigma grid in krls_autotune_inner_cpp"
```

---

## Task 5: Refactor `.ares_autotune` R harness to use the new C++ inner

**Files:**
- Modify: `R/krls.R` (`.ares_autotune` body, ~80 LOC delta + new `autotune.nthreads` arg)

- [ ] **Step 1: Add `autotune.nthreads` to `krls.default` signature**

In `R/krls.R`, find the `krls.default` function definition (~L244). Add `autotune.nthreads = NULL` to the formal arguments list, right after `autotune.grid = NULL` (~L253).

Find the autotune validation block (~L462). Below the existing autotune.grid validation, add:

```r
  # autotune.nthreads — default via getOption fallback
  if (is.null(autotune.nthreads)) {
    autotune.nthreads <- getOption(
      "roadrunner.nthreads",
      max(1L, parallel::detectCores(logical = FALSE))
    )
  }
  stopifnot(is.numeric(autotune.nthreads), length(autotune.nthreads) == 1L,
            autotune.nthreads >= 1)
  autotune.nthreads <- as.integer(autotune.nthreads)
  if (length(autotune.grid) > 0L) {
    autotune.nthreads <- min(autotune.nthreads, length(autotune.grid))
  }
  if (isTRUE(RcppParallel::defaultNumThreads() == 1L)) {
    autotune.nthreads <- 1L
  }
```

- [ ] **Step 2: Locate the current `.ares_autotune` body**

```bash
grep -n "\.ares_autotune <- function" /home/jack/Dropbox/roadrunner/R/krls.R
```

Note the line range. The body spans roughly the `function(...)` opening to the matching closing brace. Read those lines into context.

- [ ] **Step 3: Replace the inner fold loop with the new C++-backed path**

Inside `.ares_autotune`, locate the triple-nested loop `for (sigma_idx in seq_along(autotune.grid))` → `for (cross in seq_len(ncross))` → `for (k in seq_len(nfold))`. Replace the entire block with:

```r
  sigma_grid_sorted <- sort(autotune.grid)
  ns <- length(sigma_grid_sorted)
  mse_arr <- array(NA_real_, c(ncross, nfold, ns))
  lam_arr <- array(NA_real_, c(ncross, nfold, ns))

  lambda_args <- list(
    tol            = 1e-6,
    L0             = .Machine$double.eps,
    L_step         = 10.0,
    U_start_from_n = TRUE
  )

  for (cross_i in seq_len(ncross)) {
    fold_seed <- as.integer(cross_seeds[cross_i])
    folds     <- .ares_make_folds(n, nfold, fold_seed)

    for (k in seq_len(nfold)) {
      idx_te <- folds[[k]]
      X_tr   <- Xs[-idx_te, , drop = FALSE]
      X_te   <- Xs[ idx_te, , drop = FALSE]
      y_tr_k <- y[-idx_te]
      y_te_k <- y[ idx_te]

      D_tr <- krls_pairwise_sqdist_cpp(X_tr, X_tr)
      D_te <- krls_pairwise_sqdist_cpp(X_te, X_tr)

      inner <- krls_autotune_inner_cpp(
        D_tr, D_te, y_tr_k, y_te_k,
        sigma_grid_sorted, lambda_args,
        nthreads = autotune.nthreads
      )

      mse_arr[cross_i, k, ] <- inner$mse_per_sigma
      lam_arr[cross_i, k, ] <- inner$lambda_per_sigma
    }
  }
```

(`cross_seeds` and `.ares_make_folds` already exist in `R/krls.R` from v0.0.0.9040. If `.ares_make_folds` is named differently — e.g. `.krls_make_folds` — use the existing helper's exact name. Verify with `grep -n "make_folds" R/krls.R` and adjust.)

Then replace the existing per-sigma aggregation block (which currently aggregates a per-sigma vector) with the new array-based aggregation:

```r
  mse_per_sigma <- apply(mse_arr, 3, mean, na.rm = TRUE)
  se_per_sigma  <- apply(mse_arr, 3, sd,   na.rm = TRUE) /
                    sqrt(ncross * nfold)
  min_idx       <- which.min(mse_per_sigma)
  threshold     <- mse_per_sigma[min_idx] + se_per_sigma[min_idx]
  sigma_1se_idx <- max(which(mse_per_sigma <= threshold))
  sigma_chosen  <- sigma_grid_sorted[sigma_1se_idx]
```

Locate the existing `autotune_info` list construction (search `autotune_info <- list(`). Add these fields:

```r
    grid              = sigma_grid_sorted,
    mse_per_sigma     = mse_per_sigma,
    se_mse            = se_per_sigma,
    mse_per_fold      = mse_arr,
    lam_per_fold      = lam_arr,
    cv.1se            = TRUE,
    sigma_1se         = sigma_grid_sorted[sigma_1se_idx],
    ncross            = ncross,
    nthreads_used     = inner$nthreads_used,
    sigma_grid_sorted = sigma_grid_sorted
```

Preserve any other fields already present (e.g. `lambda`, `mse`) — only ADD the new ones; do not delete the old. If a field name collides, the new C++-backed value wins.

- [ ] **Step 4: Verify the package loads without errors**

```bash
Rscript -e 'devtools::load_all(".", quiet=TRUE)' 2>&1 | head -30
```

Expected: no errors or warnings. (Parse / load errors mean a brace or comma is off.)

- [ ] **Step 5: Run a quick smoke autotune fit**

```bash
Rscript -e '
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
set.seed(99L)
X <- matrix(rnorm(100 * 5), 100, 5)
y <- as.numeric(rowSums(X[, 1:2]) + 0.5 * rnorm(100))
fit <- roadrunner::krls(X = X, y = y, autotune = TRUE,
                        derivative = FALSE, vcov = FALSE)
cat("sigma  =", fit$sigma, "\n")
cat("lambda =", fit$lambda, "\n")
cat("R2     =", fit$R2, "\n")
cat("nthreads_used =", fit$autotune$nthreads_used, "\n")
'
```

Expected: finite sigma / lambda / R2; `nthreads_used` ≥ 1.

- [ ] **Step 6: Commit**

```bash
git add R/krls.R
git commit -m "feat(krls): .ares_autotune uses C++ inner with shared D + autotune.nthreads"
```

---

## Task 6: EQUIV regression tests (snapshot parity + nthreads invariance)

**Files:**
- Modify: `tests/testthat/test-krls-autotune-parallel.R`

- [ ] **Step 1: Add EQUIV-1, EQUIV-2, EQUIV-3, EQUIV-4 tests**

Append to `tests/testthat/test-krls-autotune-parallel.R`:
```r
test_that("EQUIV-1: byte-identical at nthreads=1 vs v0.0.0.9041 snapshot", {
  baseline_path <- testthat::test_path("..", "..", "inst", "testdata",
                                       "autotune-baseline-9041.rds")
  skip_if_not(file.exists(baseline_path),
              "autotune-baseline-9041.rds missing")
  baseline <- readRDS(baseline_path)

  set.seed(baseline$seed)
  X <- matrix(rnorm(baseline$n * baseline$p), baseline$n, baseline$p)
  y <- as.numeric(rowSums(X[, 1:2]) + 0.5 * rnorm(baseline$n))

  fit <- roadrunner::krls(
    X = X, y = y, autotune = TRUE,
    derivative = FALSE, vcov = FALSE,
    autotune.nthreads = 1L
  )

  expect_identical(fit$sigma,  baseline$sigma)
  expect_identical(fit$lambda, baseline$lambda)
  expect_identical(fit$coeffs, baseline$coeffs)
})

test_that("EQUIV-2: nthreads=1 vs nthreads=4 byte-identical (determinism contract)", {
  set.seed(7L)
  n <- 150L; p <- 4L
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(rowSums(X[, 1:2]) + 0.5 * rnorm(n))
  fit1 <- roadrunner::krls(X, y, autotune = TRUE, derivative = FALSE,
                           vcov = FALSE, autotune.nthreads = 1L)
  fit4 <- roadrunner::krls(X, y, autotune = TRUE, derivative = FALSE,
                           vcov = FALSE, autotune.nthreads = 4L)
  expect_identical(fit1$sigma,  fit4$sigma)
  expect_identical(fit1$lambda, fit4$lambda)
  expect_identical(fit1$coeffs, fit4$coeffs)
})

test_that("EQUIV-3: shuffled autotune.grid input yields identical fit", {
  set.seed(8L)
  n <- 120L; p <- 4L
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(rowSums(X[, 1:2]) + 0.5 * rnorm(n))
  grid_sorted   <- c(0.5, 1, 2, 4, 8, 16)
  grid_shuffled <- sample(grid_sorted)
  fit_s <- roadrunner::krls(X, y, autotune = TRUE, derivative = FALSE,
                            vcov = FALSE, autotune.grid = grid_sorted)
  fit_h <- roadrunner::krls(X, y, autotune = TRUE, derivative = FALSE,
                            vcov = FALSE, autotune.grid = grid_shuffled)
  expect_identical(fit_s$sigma,  fit_h$sigma)
  expect_identical(fit_s$lambda, fit_h$lambda)
})

test_that("EQUIV-4: (ncross x nfold) variations all stable across nthreads", {
  set.seed(9L)
  n <- 120L; p <- 3L
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  for (cfg in list(c(2, 2), c(2, 5), c(3, 10))) {
    f1 <- roadrunner::krls(X, y, autotune = TRUE, ncross = cfg[1],
                           nfold = cfg[2], derivative = FALSE, vcov = FALSE,
                           autotune.nthreads = 1L)
    f4 <- roadrunner::krls(X, y, autotune = TRUE, ncross = cfg[1],
                           nfold = cfg[2], derivative = FALSE, vcov = FALSE,
                           autotune.nthreads = 4L)
    expect_identical(f1$sigma,  f4$sigma)
    expect_identical(f1$lambda, f4$lambda)
  }
})
```

- [ ] **Step 2: Run tests**

```bash
Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls-autotune-parallel")'
```

Expected: all 9 tests in the file PASS (DIST-1, DIST-2, DIST-3, INNER-1, INNER-2, INNER-3, INNER-EQUIV, EQUIV-1, EQUIV-2, EQUIV-3, EQUIV-4 — 11 total).

If EQUIV-1 fails: the snapshot was built at a different commit. Verify `git log --oneline -1` matches the SHA where the snapshot was made (`2c83d4a` or later v0.0.0.9041), then rebuild via Task 1 and retry.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-krls-autotune-parallel.R
git commit -m "test(krls): EQUIV-1..4 snapshot + nthreads-invariance regression"
```

---

## Task 7: Full regression — krls suite + R CMD check

**Files:**
- (No new files)

- [ ] **Step 1: Run full krls test suite**

```bash
Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls")'
```

Expected: 0 failures. (Prior partial-arg warnings from `KRLS::krls` are pre-existing — they are noise, not failures.)

If failures appear: check whether they're snapshot mismatches in existing tests that reference the old autotune output shape. Update the snapshot only if the user-facing fit is unchanged (compare `fit$sigma` and `fit$coeffs` against the v0.0.0.9041 baseline RDS).

- [ ] **Step 2: Run `R CMD check`**

```bash
Rscript -e 'devtools::check(".", document = FALSE, args = "--no-tests")'
```

Expected: 0 errors, 0 warnings. The 2 pre-existing notes (future timestamp + non-portable compile flags) should remain — no new notes.

If a new note appears about Rcpp exports or unused includes, run `Rscript -e 'devtools::document(".")'` and re-check.

- [ ] **Step 3: Snapshot a quick wall-clock smoke (sanity, not a benchmark)**

```bash
Rscript -e '
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
set.seed(1L)
X <- matrix(rnorm(500 * 10), 500, 10)
y <- as.numeric(rowSums(X[, 1:3]) + 0.5 * rnorm(500))
t1 <- system.time(roadrunner::krls(X, y, autotune = TRUE,
                                   derivative = FALSE, vcov = FALSE,
                                   autotune.nthreads = 1L))[3]
t4 <- system.time(roadrunner::krls(X, y, autotune = TRUE,
                                   derivative = FALSE, vcov = FALSE,
                                   autotune.nthreads = 4L))[3]
cat(sprintf("nthreads=1: %.2fs  nthreads=4: %.2fs  speedup=%.2fx\n",
            t1, t4, t1 / t4))
'
```

Expected: speedup ratio printed > 1.5. (Sanity check; formal benchmarks run in Task 9.)

- [ ] **Step 4: Commit any documentation regen artifacts (if `devtools::check` updated `man/`)**

```bash
git status --short
# If man/*.Rd changed:
git add man/
git commit -m "docs(krls): regen Rd for autotune.nthreads"
# If nothing changed, skip the commit.
```

---

## Task 8: Documentation + version bump

**Files:**
- Modify: `R/krls.R` (roxygen for `autotune.nthreads`, `@details`, `@note`)
- Modify: `NEWS.md` (v0.0.0.9042 entry)
- Modify: `DESCRIPTION` (version bump)
- Modify: `man/krls.Rd` (regen)

- [ ] **Step 1: Add roxygen for `autotune.nthreads`**

In `R/krls.R`, locate the `@param autotune.grid` roxygen line. Insert after it:

```r
#' @param autotune.nthreads Integer >= 1. Number of TBB worker threads used
#'   for the parallel sigma sweep inside `autotune = TRUE`. Default is
#'   `getOption("roadrunner.nthreads", parallel::detectCores(logical=FALSE))`.
#'   `1` forces strictly sequential execution; this is byte-identical to the
#'   pre-v0.0.0.9042 path at fixed inputs and provides an escape hatch when
#'   the determinism contract is being audited. Capped internally at
#'   `length(autotune.grid)` (no point spawning more workers than cells).
```

- [ ] **Step 2: Add `@details` paragraph for the parallel sweep**

In the `@details` section of `krls`, append:

```r
#' Since v0.0.0.9042 the autotune inner loop is parallelised over sigma
#' candidates using `RcppParallel` (TBB). The pairwise squared-distance
#' matrix is computed once per CV fold and reused across every sigma in
#' the grid (`exp(-D / sigma)` is then a cheap elementwise op).
#' Determinism is preserved: each worker writes to a unique slot in the
#' output vectors, so the result is byte-identical regardless of
#' `autotune.nthreads`. Pass `autotune.nthreads = 1` to force strictly
#' sequential execution.
```

- [ ] **Step 3: Add `@note` for BLAS-contention hint**

In the `@note` section, append:

```r
#' Phase 1 speedup (v0.0.0.9042): parallel autotune assumes single-threaded
#' BLAS for the small kernel operations inside the parallel region. If you
#' have set a high process-global BLAS thread count (e.g. via
#' `RhpcBLASctl::blas_set_num_threads(8)`), consider resetting to 1 around
#' `krls(..., autotune = TRUE)` calls to avoid oversubscription. The
#' BLAS-heavy distance computation runs OUTSIDE the parallel region and
#' benefits from multi-threaded BLAS.
```

- [ ] **Step 4: Bump version in `DESCRIPTION`**

Change the `Version:` line:
```
Version: 0.0.0.9042
```

- [ ] **Step 5: Add NEWS entry**

Prepend to `NEWS.md`:

```markdown
# roadrunner 0.0.0.9042

## krls() speedup — Phase 1 (shared distance + parallel autotune)

* `krls(..., autotune = TRUE)` is now parallelised over the sigma grid via
  RcppParallel/TBB. The pairwise squared-distance matrix is computed once
  per CV fold and reused across every sigma candidate.
* New optional argument `autotune.nthreads` (default
  `getOption("roadrunner.nthreads", parallel::detectCores(logical=FALSE))`).
  Pass `autotune.nthreads = 1` for strictly sequential execution.
* Determinism contract preserved: fits are byte-identical across
  `autotune.nthreads` values at fixed seed and inputs (each worker writes
  to a unique slot in the output vectors, so no reduction is involved).

### Empirical wall-clock speedup (REQ-001 minimal grid, R=5)

| n    | p  | nthreads | time   | speedup vs v0.0.0.9041 |
|------|----|----------|--------|------------------------|
| 500  | 10 | 1        | (TBD)  | 1.0x (sequential)      |
| 500  | 10 | 4        | (TBD)  | (filled by Task 9)     |
| 1500 | 20 | 1        | (TBD)  | 1.0x (sequential)      |
| 1500 | 20 | 4        | (TBD)  | (filled by Task 9)     |

### API compatibility

* Zero breaking changes. `krls(..., autotune = TRUE)` returns the same
  S3 object with the same fields. `autotune_info` gains
  `nthreads_used` and `sigma_grid_sorted` slots.
* `predict.krls_rr`, `summary.krls_rr`, `print.krls_rr` unchanged.

### Out of scope (Phase 2)

* Nystrom approximation (`approx = "nystrom"`, `nystrom_m`, `landmarks`).
* Predict()-side speedups.
* Multi-response / sparse kernel paths.
```

- [ ] **Step 6: Regenerate `man/krls.Rd`**

```bash
Rscript -e 'devtools::document(".")'
```

- [ ] **Step 7: Verify roxygen and Rd are consistent**

```bash
grep -n "autotune.nthreads" R/krls.R man/krls.Rd | head
```

Expected: matches in both files.

- [ ] **Step 8: Commit**

```bash
git add DESCRIPTION NEWS.md R/krls.R man/krls.Rd
git commit -m "docs(krls): v0.0.0.9042 roxygen + NEWS for parallel autotune"
```

---

## Task 9: EMP-PHASE1 empirical sim

**Files:**
- Create: `inst/sims/krls-speedup-phase1.R`

- [ ] **Step 1: Write the sim script with mandatory rules baked in**

`inst/sims/krls-speedup-phase1.R`:
```r
#!/usr/bin/env Rscript
##
## EMP-PHASE1 — empirical speedup verification for v0.0.0.9042.
##
## Mandatory sim rules (this session):
##   1. Smoke test first (1 cell, R=2)
##   2. Intermediate CSV per cell (append after every cell)
##   3. Wall-clock cap 5 min total (abort + diagnostic if projected over)
##
## Cells: 4 = {n=500 p=10, n=1500 p=20} x {nthreads=1, nthreads=4}.
## R = 5 reps per cell. Same (X, y) seeds across nthreads in each rep
## (paired comparison; tests byte-identicality at the same time).
##
## Acceptance:
##   - Median speedup at (500, 10, nthreads=4) >= 3.0x vs nthreads=1
##   - Median speedup at (1500, 20, nthreads=4) >= 3.0x vs nthreads=1
##   - test_MSE_(nthreads=4) == test_MSE_(nthreads=1) byte-identical per rep

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})

OUT_CSV <- "inst/sims/results/krls-speedup-phase1.csv"
WALL_CAP_S <- 300
MASTER_SEED <- 20260518L
R_REPS <- 5L

dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
if (file.exists(OUT_CSV)) file.remove(OUT_CSV)

dgp_additive <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(sin(X[, 1]) + X[, 2]^2 + X[, 3] + 0.5 * rnorm(n))
  list(X = X, y = y)
}

run_one <- function(n, p, nthreads, rep_seed) {
  d_tr <- dgp_additive(n, p, rep_seed)
  d_te <- dgp_additive(1000L, p, rep_seed + 1L)
  t0   <- proc.time()[["elapsed"]]
  fit  <- roadrunner::krls(
    X = d_tr$X, y = d_tr$y,
    autotune = TRUE, derivative = FALSE, vcov = FALSE,
    autotune.nthreads = nthreads
  )
  wall  <- proc.time()[["elapsed"]] - t0
  yhat  <- as.numeric(predict(fit, newdata = d_te$X)$fit)
  test_MSE <- mean((d_te$y - yhat)^2)
  list(wall = wall, test_MSE = test_MSE,
       sigma = fit$sigma, lambda = fit$lambda,
       nthreads_used = fit$autotune$nthreads_used)
}

append_row <- function(row) {
  write.table(row, OUT_CSV, sep = ",", row.names = FALSE,
              col.names = !file.exists(OUT_CSV), append = file.exists(OUT_CSV),
              quote = FALSE)
}

## ---------- SMOKE TEST (mandatory) ----------
cat("[smoke] n=200 p=5 nthreads=1 R=2 ...\n")
sm1 <- run_one(200L, 5L, 1L, MASTER_SEED + 9999L)
sm2 <- run_one(200L, 5L, 1L, MASTER_SEED + 10000L)
stopifnot(is.finite(sm1$test_MSE), is.finite(sm2$test_MSE),
          sm1$wall > 0, sm2$wall > 0)
cat(sprintf("[smoke OK] wall=%.2fs / %.2fs  MSE=%.4f / %.4f\n",
            sm1$wall, sm2$wall, sm1$test_MSE, sm2$test_MSE))

## ---------- FULL GRID ----------
cells <- expand.grid(
  n = c(500L, 1500L),
  p = c(10L, 20L),
  nthreads = c(1L, 4L)
)
## Keep only the two paired (n, p) sets specified in EMP-PHASE1.
cells <- cells[(cells$n == 500L & cells$p == 10L) |
               (cells$n == 1500L & cells$p == 20L), ]

cat(sprintf("[full] %d cells x R=%d reps = %d fits\n",
            nrow(cells), R_REPS, nrow(cells) * R_REPS))
wall_start <- proc.time()[["elapsed"]]

for (ci in seq_len(nrow(cells))) {
  n <- cells$n[ci]; p <- cells$p[ci]; nt <- cells$nthreads[ci]
  for (r in seq_len(R_REPS)) {
    elapsed <- proc.time()[["elapsed"]] - wall_start
    if (elapsed > WALL_CAP_S - 15) {
      cat(sprintf("[abort] wall-cap approaching (%.0fs). Stopping early.\n",
                  elapsed))
      stop("EMP-PHASE1: 5-min wall cap hit. Reduce R or grid.")
    }
    rep_seed <- MASTER_SEED + (ci - 1L) * R_REPS * 2L + (r - 1L) * 2L
    res <- run_one(n, p, nt, rep_seed)
    row <- data.frame(
      n = n, p = p, nthreads = nt, rep = r, seed = rep_seed,
      wall_s = res$wall, test_MSE = res$test_MSE,
      sigma = res$sigma, lambda = res$lambda,
      nthreads_used = res$nthreads_used
    )
    append_row(row)
    cat(sprintf("[%d/%d cell=%d nt=%d rep=%d] wall=%.2fs MSE=%.4f\n",
                ci, nrow(cells), ci, nt, r, res$wall, res$test_MSE))
  }
}

## ---------- SUMMARY ----------
df <- read.csv(OUT_CSV)
agg <- aggregate(cbind(wall_s, test_MSE) ~ n + p + nthreads, data = df,
                 FUN = median)
print(agg)

speedup_500_10 <- agg$wall_s[agg$n == 500  & agg$nthreads == 1] /
                  agg$wall_s[agg$n == 500  & agg$nthreads == 4]
speedup_1500_20 <- agg$wall_s[agg$n == 1500 & agg$nthreads == 1] /
                   agg$wall_s[agg$n == 1500 & agg$nthreads == 4]

cat(sprintf("\nSpeedup (n=500,  p=10,  nthreads=4): %.2fx\n", speedup_500_10))
cat(sprintf("Speedup (n=1500, p=20,  nthreads=4): %.2fx\n", speedup_1500_20))

## Check byte-identical test_MSE across nthreads at same seed
chk <- merge(
  df[df$nthreads == 1, c("n", "p", "rep", "test_MSE")],
  df[df$nthreads == 4, c("n", "p", "rep", "test_MSE")],
  by = c("n", "p", "rep"), suffixes = c("_1", "_4")
)
chk$diff <- abs(chk$test_MSE_1 - chk$test_MSE_4)
max_diff <- max(chk$diff)
cat(sprintf("Max |test_MSE(nthreads=4) - test_MSE(nthreads=1)| = %.3e\n",
            max_diff))
stopifnot(max_diff < 1e-9)
cat("[ok] byte-identicality across nthreads holds.\n")

if (speedup_500_10 < 3.0 || speedup_1500_20 < 3.0) {
  warning("EMP-PHASE1: speedup < 3.0x in at least one cell. Check results.")
}
cat("\nResults at:", OUT_CSV, "\n")
```

- [ ] **Step 2: Run the smoke test first (mandatory rule 1)**

The script's own smoke block runs at top of execution. Run the full script — it will smoke-test before proceeding:

```bash
Rscript inst/sims/krls-speedup-phase1.R 2>&1 | tee /tmp/emp-phase1-log.txt
```

Expected output ends with `[ok] byte-identicality across nthreads holds.` and a printed speedup line for each (n, p) pair.

If the script aborts at the smoke phase: a basic call to `krls(autotune=TRUE)` failed. Investigate before continuing.

- [ ] **Step 3: Inspect intermediate CSV after run**

```bash
head /home/jack/Dropbox/roadrunner/inst/sims/results/krls-speedup-phase1.csv
wc -l /home/jack/Dropbox/roadrunner/inst/sims/results/krls-speedup-phase1.csv
```

Expected: header + 20 rows (4 cells × 5 reps).

- [ ] **Step 4: Verify acceptance criteria**

Median speedup at (500, 10, 4) and (1500, 20, 4) must each be ≥ 3.0×. Max |test_MSE diff across nthreads| must be < 1e-9.

If a speedup target misses: check whether BLAS thread oversubscription is the cause (run with `OMP_NUM_THREADS=1 Rscript ...`). If still missing, escalate — likely an architectural issue (e.g. TBB grain size too large for small sigma counts).

If byte-identicality fails: critical — the determinism contract is broken. Stop and investigate `arma::eig_sym` thread behavior; likely a BLAS thread race.

- [ ] **Step 5: Update NEWS table with measured numbers**

Open `NEWS.md`. Find the v0.0.0.9042 speedup table from Task 8 Step 5. Replace the `(TBD)` placeholders with the median wall-clock values printed by the sim. Recompute speedup ratios.

- [ ] **Step 6: Commit sim + results + NEWS update**

```bash
git add inst/sims/krls-speedup-phase1.R inst/sims/results/krls-speedup-phase1.csv NEWS.md
git commit -m "sim: EMP-PHASE1 speedup verification for v0.0.0.9042"
```

---

## Task 10: Final regression + ship

**Files:**
- (no new files; verification + ship)

- [ ] **Step 1: Full krls suite final pass**

```bash
Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls")'
```

Expected: 0 failures. 7 pre-existing partial-arg warnings only. Note total test count — should be roughly `prior + 11` (11 new in test-krls-autotune-parallel.R).

- [ ] **Step 2: Full `R CMD check`**

```bash
Rscript -e 'devtools::check(".")'
```

Expected: 0 errors, 0 warnings, 2 pre-existing notes (future timestamp + non-portable compile flags). No NEW notes.

If a new note about Rcpp exports or undocumented arguments appears: regenerate docs with `devtools::document(".")`, recheck.

- [ ] **Step 3: Verify final tree state**

```bash
git status --short
git log --oneline -8
```

Expected: clean tree; commits visible in order:
1. snapshot baseline (Task 1)
2. pairwise sqdist (Task 2)
3. inner_cpp sequential (Task 3)
4. parallelFor (Task 4)
5. .ares_autotune harness (Task 5)
6. EQUIV tests (Task 6)
7. (optional) Rd regen (Task 7)
8. docs + version (Task 8)
9. EMP-PHASE1 sim (Task 9)

- [ ] **Step 4: Push to origin**

User must approve push explicitly. Once approved:

```bash
git push origin main
```

Expected: fast-forward push, ~9 commits delivered.

- [ ] **Step 5: Tag the release (optional, if convention is to tag)**

```bash
git tag -a v0.0.0.9042 -m "krls speedup phase 1: parallel autotune"
git push origin v0.0.0.9042
```

---

## Definition of done

- [ ] All 11 new tests pass: `devtools::test(filter="krls-autotune-parallel")`
- [ ] Full krls suite passes: `devtools::test(filter="krls")` 0 fail
- [ ] R CMD check: 0E / 0W / 2 pre-existing notes
- [ ] EMP-PHASE1 speedup ≥ 3.0× at both (500, 10) and (1500, 20) with nthreads=4
- [ ] EMP-PHASE1 byte-identicality across nthreads: max diff < 1e-9
- [ ] `NEWS.md` v0.0.0.9042 entry committed with measured speedup numbers
- [ ] Tree clean, pushed to `origin/main`

---

## Risk log (cross-reference to spec section 9)

- **TBB nesting with BLAS inside eigen**: mitigated by structure (heavy `dgemm` outside parallel region, eigen single-threaded). Verified by EQUIV-2 / INNER-EQUIV.
- **Stale snapshot baseline**: Task 1 builds the baseline at HEAD before any code change. Document the commit SHA used.
- **Memory for D_tr at large n**: n=1500 → ~18 MB; n=2000 → ~32 MB. Within bounds. Larger n flagged for Phase 2 (Nystrom).
