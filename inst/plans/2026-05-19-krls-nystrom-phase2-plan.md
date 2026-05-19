# roadrunner::krls Phase 2 (Nystrom) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add opt-in `approx="nystrom"` to `roadrunner::krls()` for ≥5× single-fit speedup at n=2000 and ≥10× at n=5000, with RMSE within +10% of exact at the default m = ⌈√n × 2⌉.

**Architecture:** Replace O(n³) eigendecomp on n×n K with O(m³) eigendecomp on m×m anchor W plus O(nm²) feature map. Landmarks drawn from training X (random or k-means). Per-fold landmarks in autotune. Three new C++ exports: `krls_nystrom_fit_cpp`, `krls_nystrom_autotune_inner_cpp`, `krls_nystrom_predict_cpp`. R harness gains five optional args (`approx`, `nystrom_m`, `landmarks`, `landmark_method`, `landmark_seed`, `nystrom_eps`).

**Tech Stack:** R, Rcpp, RcppArmadillo (`eig_sym`, `svd`, `dgemm`), RcppParallel (TBB — for autotune inner only), testthat 3, devtools.

**Spec:** `inst/specs/2026-05-19-krls-nystrom-design.md` (commit `c220a15`)

**Target version:** `0.0.0.9043`

**Reference impl** (read-only): `/tmp/krls-reference/R/nystrom.R` (380 LOC pure-R).

---

## Pre-flight

Verify clean tree on `main` at HEAD `c220a15` (the spec commit).

```bash
git -C /home/jack/Dropbox/roadrunner status --short
git -C /home/jack/Dropbox/roadrunner log --oneline -3
```

Expected: tree clean; `c220a15` visible.

---

## Task 1: `.resolve_landmarks` R helper + landmark-resolution tests

**Files:**
- Create: `R/krls-nystrom.R` (new file dedicated to Nystrom internals)
- Create: `tests/testthat/test-krls-nystrom.R`

- [ ] **Step 1: Write failing tests**

Create `tests/testthat/test-krls-nystrom.R`:
```r
## Phase 2 (v0.0.0.9043) tests: Nystrom approximation.
## Mirrors inst/specs/2026-05-19-krls-nystrom-design.md.

test_that("LAND-RAND-1: NULL landmarks with random method returns valid indices", {
  set.seed(1L)
  n <- 100L; d <- 3L
  X <- matrix(rnorm(n * d), n, d)
  out <- roadrunner:::.resolve_landmarks(
    landmarks = NULL, landmark_method = "random", nystrom_m = 25L,
    X_std = X, X_centers = rep(0, d), X_scales = rep(1, d),
    landmark_seed = NULL
  )
  expect_equal(nrow(out$matrix), 25L)
  expect_equal(ncol(out$matrix), d)
  expect_length(out$indices, 25L)
  expect_true(all(out$indices >= 1L & out$indices <= n))
  expect_false(anyDuplicated(out$indices) > 0)
  expect_identical(out$method_used, "random")
})

test_that("LAND-RAND-SEED-1: same landmark_seed produces identical indices", {
  set.seed(99L)
  n <- 100L; d <- 3L
  X <- matrix(rnorm(n * d), n, d)
  out1 <- roadrunner:::.resolve_landmarks(NULL, "random", 25L, X,
                                          rep(0, d), rep(1, d),
                                          landmark_seed = 42L)
  out2 <- roadrunner:::.resolve_landmarks(NULL, "random", 25L, X,
                                          rep(0, d), rep(1, d),
                                          landmark_seed = 42L)
  expect_identical(out1$indices, out2$indices)
  expect_identical(out1$matrix,  out2$matrix)
})

test_that("LAND-IDX-1: integer vector input returns matching X rows", {
  set.seed(2L)
  n <- 50L; d <- 4L
  X <- matrix(rnorm(n * d), n, d)
  idx <- c(3L, 7L, 12L, 25L, 40L)
  out <- roadrunner:::.resolve_landmarks(idx, "random", 5L, X,
                                         rep(0, d), rep(1, d), NULL)
  expect_identical(out$indices, idx)
  expect_equal(out$matrix, X[idx, , drop = FALSE])
  expect_identical(out$method_used, "user_indices")

  # Out-of-range indices error
  expect_error(
    roadrunner:::.resolve_landmarks(c(0L, 5L), "random", 2L, X,
                                    rep(0, d), rep(1, d), NULL),
    "unique integers in 1"
  )
})

test_that("LAND-MAT-1: numeric matrix input is standardized", {
  set.seed(3L)
  n <- 50L; d <- 3L
  X <- matrix(rnorm(n * d), n, d)
  X_centers <- colMeans(X)
  X_scales  <- apply(X, 2, sd)
  X_std     <- scale(X, center = X_centers, scale = X_scales)
  attr(X_std, "scaled:center") <- NULL
  attr(X_std, "scaled:scale")  <- NULL

  # Pass landmarks in original scale; helper should standardize internally
  Z_orig <- X[c(1, 5, 10), , drop = FALSE]
  out <- roadrunner:::.resolve_landmarks(Z_orig, "random", 3L, X_std,
                                         X_centers, X_scales, NULL)
  expect_null(out$indices)
  expect_equal(out$matrix, X_std[c(1, 5, 10), , drop = FALSE],
               tolerance = 1e-12)
  expect_identical(out$method_used, "user_matrix")

  # ncol mismatch error
  expect_error(
    roadrunner:::.resolve_landmarks(matrix(0, 3, 2), "random", 3L, X_std,
                                    X_centers, X_scales, NULL),
    "ncol\\(landmarks\\) must equal"
  )
})

test_that("LAND-VALIDATE-1: nystrom_m and nystrom_eps validators reject bad input", {
  expect_error(roadrunner:::.validate_nystrom_m(NA, 100), "finite integer")
  expect_error(roadrunner:::.validate_nystrom_m(0L, 100), "1 <= nystrom_m")
  expect_error(roadrunner:::.validate_nystrom_m(101L, 100), "1 <= nystrom_m")
  expect_error(roadrunner:::.validate_nystrom_eps(0), "positive scalar")
  expect_error(roadrunner:::.validate_nystrom_eps(-1), "positive scalar")
  expect_identical(roadrunner:::.validate_nystrom_eps(1e-9), 1e-9)
})
```

- [ ] **Step 2: Run tests to verify failure**

```bash
Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls-nystrom")'
```

Expected: 5 ERRORS, `could not find function "..."`.

- [ ] **Step 3: Create `R/krls-nystrom.R` with the resolve_landmarks helpers**

```r
# Internal helpers for krls(..., approx = "nystrom").
# Phase 2 (v0.0.0.9043). See inst/specs/2026-05-19-krls-nystrom-design.md.
#
# Nystrom replaces the full n x n kernel by a low-rank approximation
# anchored at m << n landmarks. The ridge problem then lives in
# m-dimensional space:
#
#   Phi      = C %*% U %*% diag(D_reg^{-1/2})       (n x m feature map)
#   beta_hat = (Phi' Phi + lambda I)^{-1} Phi' y
#   f_hat(x) = K(x, Z) %*% alpha,   alpha = U %*% (D_reg^{-1/2} * beta)
#
# where C = K(X, Z), W = K(Z, Z), and (U, D) is the eigendecomposition
# of W with relative-ridge stabilization D_reg = pmax(D, eps * max(D)).

.select_landmarks_random <- function(n, m) {
  sort.int(sample.int(n, m))
}

# Run `expr` with the RNG temporarily seeded at `seed`, restoring the
# caller's .Random.seed (or its absence) on exit. Used so that
# passing `landmark_seed` to krls() does NOT perturb the user's
# downstream RNG state.
.with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv,
                                inherits = FALSE) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv,
                      inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  })
  set.seed(seed)
  expr
}

.validate_nystrom_m <- function(nystrom_m, n) {
  if (!is.numeric(nystrom_m) || length(nystrom_m) != 1L ||
      is.na(nystrom_m) || !is.finite(nystrom_m) ||
      nystrom_m != as.integer(nystrom_m)) {
    stop("nystrom_m must be a finite integer")
  }
  nystrom_m <- as.integer(nystrom_m)
  if (nystrom_m < 1L || nystrom_m > n) {
    stop("nystrom_m must satisfy 1 <= nystrom_m <= nrow(X)")
  }
  nystrom_m
}

.validate_nystrom_eps <- function(nystrom_eps) {
  if (!is.numeric(nystrom_eps) || length(nystrom_eps) != 1L ||
      is.na(nystrom_eps) || !is.finite(nystrom_eps) ||
      nystrom_eps <= 0) {
    stop("nystrom_eps must be a finite positive scalar")
  }
  nystrom_eps
}

# X_std is the already-standardized training matrix; X_centers / X_scales
# are the column means / SDs of the original-scale training X.
# Index- and NULL-form landmarks resolve via X_std rows (already standardized).
# Matrix-form landmarks are taken as user-supplied original-scale
# coordinates and standardized with the same centers/scales.
.resolve_landmarks <- function(landmarks, landmark_method, nystrom_m,
                               X_std, X_centers, X_scales,
                               landmark_seed = NULL) {
  n <- nrow(X_std)
  d <- ncol(X_std)

  if (is.null(landmarks)) {
    if (is.null(nystrom_m)) {
      nystrom_m <- ceiling(sqrt(n) * 2)
    }
    nystrom_m <- .validate_nystrom_m(nystrom_m, n)
    return(.with_seed(landmark_seed, {
      if (identical(landmark_method, "kmeans")) {
        if (nystrom_m == n) {
          list(indices = seq_len(n),
               matrix  = X_std,
               method_used = "kmeans")
        } else {
          km <- stats::kmeans(X_std, centers = nystrom_m,
                              nstart = 10L, iter.max = 50L)
          Z  <- unname(km$centers)
          list(indices = NULL, matrix = Z, method_used = "kmeans")
        }
      } else {
        idx <- .select_landmarks_random(n, nystrom_m)
        list(indices = idx,
             matrix  = X_std[idx, , drop = FALSE],
             method_used = "random")
      }
    }))
  }

  if (is.numeric(landmarks) && is.null(dim(landmarks))) {
    if (anyNA(landmarks) || any(!is.finite(landmarks)) ||
        any(landmarks != as.integer(landmarks))) {
      stop("landmarks (as indices) must be integer-valued")
    }
    idx <- as.integer(landmarks)
    if (anyNA(idx) || any(idx < 1L) || any(idx > n) || anyDuplicated(idx)) {
      stop("landmarks (as indices) must be unique integers in 1:nrow(X)")
    }
    return(list(indices = idx,
                matrix  = X_std[idx, , drop = FALSE],
                method_used = "user_indices"))
  }

  if (is.matrix(landmarks) || is.data.frame(landmarks)) {
    Z <- as.matrix(landmarks)
    if (!is.numeric(Z))       stop("landmarks matrix must be numeric")
    if (ncol(Z) != d)         stop("ncol(landmarks) must equal ncol(X)")
    if (anyNA(Z) || any(!is.finite(Z)))
      stop("landmarks matrix must contain only finite, non-NA values")
    if (anyDuplicated(Z))
      stop("landmarks matrix must not contain duplicate rows ",
           "(W would be rank-deficient)")
    Z_std <- scale(Z, center = X_centers, scale = X_scales)
    attr(Z_std, "scaled:center") <- NULL
    attr(Z_std, "scaled:scale")  <- NULL
    return(list(indices = NULL, matrix = Z_std, method_used = "user_matrix"))
  }

  stop("`landmarks` must be NULL, an integer vector of row indices, ",
       "or an m x d numeric matrix (in original X-scale)")
}
```

- [ ] **Step 4: Run tests**

```bash
cd /home/jack/Dropbox/roadrunner && Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls-nystrom")'
```

Expected: all 5 tests PASS.

- [ ] **Step 5: Commit**

```bash
cd /home/jack/Dropbox/roadrunner && git add R/krls-nystrom.R tests/testthat/test-krls-nystrom.R && git commit -m "feat(krls): .resolve_landmarks helper for Nystrom approximation"
```

---

## Task 2: `krls_nystrom_predict_cpp` + `get_landmarks` accessor

**Files:**
- Modify: `src/krls.cpp` (~40 LOC)
- Modify: `R/krls-nystrom.R` (add `get_landmarks` exported function)
- Modify: `NAMESPACE` (export `get_landmarks`)
- Modify: `tests/testthat/test-krls-nystrom.R` (append NYS-XK-1, NYS-GETLM-1)

- [ ] **Step 1: Write failing tests**

Append to `tests/testthat/test-krls-nystrom.R`:
```r
test_that("NYS-XK-1: krls_nystrom_predict_cpp matches naive reference", {
  set.seed(10L)
  n_new <- 20L; m <- 8L; d <- 3L
  Xn <- matrix(rnorm(n_new * d), n_new, d)
  Z  <- matrix(rnorm(m * d), m, d)
  alpha <- rnorm(m)
  sigma <- 5.0

  yhat_cpp <- roadrunner:::krls_nystrom_predict_cpp(Xn, Z, alpha, sigma)

  # Naive reference
  K_new <- matrix(NA_real_, n_new, m)
  for (i in seq_len(n_new)) for (j in seq_len(m)) {
    K_new[i, j] <- exp(-sum((Xn[i, ] - Z[j, ])^2) / sigma)
  }
  yhat_ref <- as.numeric(K_new %*% alpha)
  expect_equal(yhat_cpp, yhat_ref, tolerance = 1e-10)
})

test_that("NYS-GETLM-1: get_landmarks returns landmarks in requested scale", {
  # Construct a minimal fake Nystrom fit
  set.seed(11L)
  n <- 50L; d <- 3L
  X <- matrix(rnorm(n * d), n, d)
  X_centers <- colMeans(X)
  X_scales  <- apply(X, 2, sd)
  Z_std <- X[c(1, 5, 10), , drop = FALSE]
  Z_std <- scale(Z_std, center = X_centers, scale = X_scales)
  attr(Z_std, "scaled:center") <- NULL
  attr(Z_std, "scaled:scale")  <- NULL

  fit <- list(approx = "nystrom", landmarks = Z_std,
              X_means = X_centers, X_sds = X_scales)
  class(fit) <- "krls_rr"

  z_std <- get_landmarks(fit, scale = "standardized")
  z_orig <- get_landmarks(fit, scale = "original")
  expect_equal(z_std, Z_std)
  expect_equal(z_orig,
               sweep(sweep(Z_std, 2, X_scales, `*`), 2, X_centers, `+`),
               tolerance = 1e-12)

  # Errors on exact fit
  fit2 <- list(approx = "exact"); class(fit2) <- "krls_rr"
  expect_error(get_landmarks(fit2), "not built with approx")
})
```

- [ ] **Step 2: Run to verify failure**

```bash
cd /home/jack/Dropbox/roadrunner && Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls-nystrom")'
```

Expected: NYS-XK-1, NYS-GETLM-1 ERROR.

- [ ] **Step 3: Implement `krls_nystrom_predict_cpp` in `src/krls.cpp`**

Append to `src/krls.cpp` (after Phase 1's exports):

```cpp
// ---------------------------------------------------------------------------
// Phase 2 (v0.0.0.9043): Nystrom predict.
//
//   yhat = exp(-||X_new - Z||^2 / sigma) @ alpha
//
// Uses the dgemm identity in krls_pairwise_sqdist_cpp for the n_new x m
// cross-distance matrix, then elementwise exp + matrix-vector product.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
arma::vec krls_nystrom_predict_cpp(const arma::mat& X_new, const arma::mat& Z,
                                   const arma::vec& alpha, double sigma) {
  if (X_new.n_cols != Z.n_cols) {
    Rcpp::stop("krls_nystrom_predict_cpp: X_new and Z must have same n_cols");
  }
  if (alpha.n_elem != Z.n_rows) {
    Rcpp::stop("krls_nystrom_predict_cpp: length(alpha) must equal nrow(Z)");
  }
  arma::mat D = krls_pairwise_sqdist_cpp(X_new, Z);   // n_new x m
  arma::mat K = arma::exp(-D / sigma);
  return K * alpha;
}
```

- [ ] **Step 4: Add `get_landmarks` to `R/krls-nystrom.R`**

Append to `R/krls-nystrom.R`:
```r
#' Extract landmark coordinates from a Nystrom krls fit
#'
#' Returns the m x d matrix of landmark coordinates used by a fit built
#' with `approx = "nystrom"`. Internally landmarks live in standardized
#' X-space; the default `scale = "original"` undoes the standardization
#' using the training X's centers and SDs.
#'
#' @param fit A `krls_rr` fit built with `approx = "nystrom"`.
#' @param scale Either `"original"` (default) or `"standardized"`.
#' @return Numeric matrix of dimension `nystrom_m` x `ncol(X_train)`.
#' @export
get_landmarks <- function(fit, scale = c("original", "standardized")) {
  if (!inherits(fit, "krls_rr")) {
    stop("fit is not of class 'krls_rr'")
  }
  if (is.null(fit$approx) || !identical(fit$approx, "nystrom") ||
      is.null(fit$landmarks)) {
    stop("fit was not built with approx = 'nystrom'; no landmarks to return")
  }
  scale <- match.arg(scale)
  if (scale == "standardized") return(fit$landmarks)
  out <- sweep(sweep(fit$landmarks, 2, fit$X_sds, `*`),
               2, fit$X_means, `+`)
  out
}
```

- [ ] **Step 5: Export `get_landmarks`**

Add to `NAMESPACE` (or rely on roxygen `@export` regen):
```r
export(get_landmarks)
```

Run `devtools::document()` to regen NAMESPACE if using roxygen.

- [ ] **Step 6: Compile + run**

```bash
cd /home/jack/Dropbox/roadrunner && Rscript -e 'devtools::document("."); Rcpp::compileAttributes(".")'
cd /home/jack/Dropbox/roadrunner && Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls-nystrom")'
```

Expected: all 7 tests PASS (5 + 2 new).

- [ ] **Step 7: Commit**

```bash
cd /home/jack/Dropbox/roadrunner && git add src/krls.cpp src/RcppExports.cpp R/RcppExports.R R/krls-nystrom.R NAMESPACE && git commit -m "feat(krls): krls_nystrom_predict_cpp + get_landmarks accessor"
```

(Stage `man/` if any new Rd was generated for `get_landmarks`.)

---

## Task 3: `krls_nystrom_fit_cpp` — single-fit Nystrom path

**Files:**
- Modify: `src/krls.cpp` (~200 LOC for fit_cpp)
- Modify: `tests/testthat/test-krls-nystrom.R` (append NYS-FIT-1, NYS-FIT-2)

- [ ] **Step 1: Write failing tests**

Append to `tests/testthat/test-krls-nystrom.R`:
```r
test_that("NYS-FIT-1: krls_nystrom_fit_cpp matches reference R impl", {
  # Reference impl lives at /tmp/krls-reference/R/nystrom.R
  ref_path <- "/tmp/krls-reference/R/nystrom.R"
  skip_if_not(file.exists(ref_path), "reference Nystrom impl missing")
  ref_env <- new.env()
  sys.source(ref_path, envir = ref_env)

  set.seed(12L)
  n <- 80L; d <- 4L; m <- 16L
  X <- matrix(rnorm(n * d), n, d)
  y <- as.numeric(rowSums(X[, 1:2]) + 0.3 * rnorm(n))
  Z <- X[seq(1L, n, length.out = m), , drop = FALSE]
  sigma <- 6.0; eps <- 1e-9

  lambda_args <- list(tol = 1e-6,
                      L0  = .Machine$double.eps,
                      L_step = 10.0,
                      U_start_from_n = TRUE)

  out_cpp <- roadrunner:::krls_nystrom_fit_cpp(
    X, Z, y, sigma, lambda_args, eps, compute_vcov = FALSE
  )

  # Reference
  ref_resolved <- list(matrix = Z, indices = NULL, method_used = "user")
  out_ref <- ref_env$.fit_krls_nystrom(
    X, y, sigma, ref_resolved,
    lambda = NULL, nystrom_eps = eps,
    L = NULL, U = NULL, tol = 1e-6, noisy = FALSE,
    compute_vcov = FALSE, lambda_method = "loo"
  )

  expect_equal(as.numeric(out_cpp$coeffs), as.numeric(out_ref$coeffs),
               tolerance = 1e-7)
  expect_equal(out_cpp$lambda, out_ref$lambda, tolerance = 1e-6)
})

test_that("NYS-FIT-2: byte-identical fits at fixed sigma+landmarks", {
  set.seed(13L)
  n <- 60L; d <- 3L; m <- 12L
  X <- matrix(rnorm(n * d), n, d)
  y <- rnorm(n)
  Z <- X[1:m, , drop = FALSE]
  sigma <- 4.0

  lambda_args <- list(tol = 1e-6, L0 = .Machine$double.eps,
                      L_step = 10.0, U_start_from_n = TRUE)

  o1 <- roadrunner:::krls_nystrom_fit_cpp(X, Z, y, sigma, lambda_args,
                                          1e-9, compute_vcov = FALSE)
  o2 <- roadrunner:::krls_nystrom_fit_cpp(X, Z, y, sigma, lambda_args,
                                          1e-9, compute_vcov = FALSE)
  expect_identical(o1$coeffs, o2$coeffs)
  expect_identical(o1$lambda, o2$lambda)
})
```

- [ ] **Step 2: Run to verify failure**

```bash
cd /home/jack/Dropbox/roadrunner && Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls-nystrom")'
```

Expected: NYS-FIT-1, NYS-FIT-2 ERROR (`could not find function "krls_nystrom_fit_cpp"`).

- [ ] **Step 3: Implement `krls_nystrom_fit_cpp` in `src/krls.cpp`**

Append:

```cpp
// ---------------------------------------------------------------------------
// Phase 2 (v0.0.0.9043): Nystrom single-fit at fixed sigma.
//
// Builds C = K(X, Z), W = K(Z, Z), eig_sym(W) -> (U, D), regularize via
// nystrom_eps relative ridge, build Phi = C @ U @ diag(D_reg^{-1/2}),
// SVD of Phi for the LOO lambda objective, golden-section search,
// recover alpha = U @ (Dinvsqrt * beta).
//
// Reference: /tmp/krls-reference/R/nystrom.R::.fit_krls_nystrom
// ---------------------------------------------------------------------------

static double nystrom_loo_loss(const arma::mat& U_phi,
                               const arma::vec& Sigma2,
                               const arma::vec& y,
                               double lambda) {
  arma::vec w     = Sigma2 / (Sigma2 + lambda);
  arma::vec Uty   = U_phi.t() * y;
  arma::vec yfit  = U_phi * (w % Uty);
  arma::vec diagS = arma::sum(arma::square(U_phi).each_row() % w.t(), 1);
  // ^ row-i sum of (U_phi(i,:)^2 elementwise * w)
  arma::vec denom = 1.0 - diagS;
  if (arma::any(denom <= std::numeric_limits<double>::epsilon())) {
    return std::numeric_limits<double>::max();
  }
  arma::vec resid = (y - yfit) / denom;
  return arma::dot(resid, resid);
}

// [[Rcpp::export]]
Rcpp::List krls_nystrom_fit_cpp(const arma::mat& X_tr, const arma::mat& Z,
                                const arma::vec& y_tr, double sigma,
                                Rcpp::List lambda_args, double eps,
                                bool compute_vcov) {
  const arma::uword n = X_tr.n_rows;
  const arma::uword m = Z.n_rows;

  // C = K(X_tr, Z), W = K(Z, Z)
  arma::mat D_tr_Z = krls_pairwise_sqdist_cpp(X_tr, Z);
  arma::mat D_Z_Z  = krls_pairwise_sqdist_cpp(Z,    Z);
  arma::mat C = arma::exp(-D_tr_Z / sigma);
  arma::mat W = arma::exp(-D_Z_Z  / sigma);

  // eig_sym(W) -> (U, D)
  arma::vec D;
  arma::mat U;
  arma::eig_sym(D, U, W);
  D = arma::clamp(D, 0.0, arma::datum::inf);
  double Dmax = D.max();
  if (Dmax <= 0.0) {
    Rcpp::stop("krls_nystrom_fit_cpp: anchor W numerically zero; try larger sigma");
  }
  arma::vec D_reg   = arma::clamp(D, eps * Dmax, arma::datum::inf);
  arma::vec Dinv2   = 1.0 / arma::sqrt(D_reg);
  arma::uword floored_count = arma::sum(D < eps * Dmax);
  double D_min_raw = D.min();

  // Phi = C @ U @ diag(Dinv2)   (n x m)
  arma::mat Phi = C * U;
  Phi.each_row() %= Dinv2.t();

  // SVD of Phi
  arma::mat U_phi, V_phi;
  arma::vec sigma_phi;
  bool ok = arma::svd(U_phi, sigma_phi, V_phi, Phi);
  if (!ok) {
    Rcpp::stop("krls_nystrom_fit_cpp: SVD of Phi failed");
  }
  arma::vec Sigma2 = arma::square(sigma_phi);
  if (Sigma2.max() <= 0.0) {
    Rcpp::stop("krls_nystrom_fit_cpp: Nystrom feature spectrum is zero; "
               "increase sigma or change landmarks");
  }

  // Golden-section lambda search bounds (anchored at Sigma2)
  double L0     = Rcpp::as<double>(lambda_args["L0"]);
  double L_step = Rcpp::as<double>(lambda_args["L_step"]);
  double tol    = Rcpp::as<double>(lambda_args["tol"]);

  double Smax = Sigma2.max();
  double U_lam = Smax;
  int iter = 0;
  while (arma::sum(Sigma2 / (Sigma2 + U_lam)) >= 1.0 && iter < 200) {
    U_lam *= 2.0;
    ++iter;
  }
  // L: index q = which.min(|Sigma2 - Smax/1000|), 1-based equivalent
  double tgt = Smax / 1000.0;
  arma::uword q_idx = 0;
  double mind = std::abs(Sigma2(0) - tgt);
  for (arma::uword i = 1; i < Sigma2.n_elem; ++i) {
    double dd = std::abs(Sigma2(i) - tgt);
    if (dd < mind) { mind = dd; q_idx = i; }
  }
  double q_R = (double)(q_idx + 1);
  double L_lam = L0;
  while (arma::sum(Sigma2 / (Sigma2 + L_lam)) > q_R && L_lam < Smax) {
    L_lam *= L_step;
  }
  if (!(L_lam < U_lam)) {
    Rcpp::stop("krls_nystrom_fit_cpp: lambda bounds collapsed; "
               "supply L and U explicitly");
  }

  // Golden-section
  const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
  double X1 = L_lam + (1.0 - gr) * (U_lam - L_lam);
  double X2 = L_lam + gr * (U_lam - L_lam);
  double S1 = nystrom_loo_loss(U_phi, Sigma2, y_tr, X1);
  double S2 = nystrom_loo_loss(U_phi, Sigma2, y_tr, X2);
  int it = 0;
  while (std::abs(S1 - S2) > tol && it < 200) {
    if (S1 < S2) {
      U_lam = X2; X2 = X1; X1 = L_lam + (1.0 - gr) * (U_lam - L_lam);
      S2 = S1; S1 = nystrom_loo_loss(U_phi, Sigma2, y_tr, X1);
    } else {
      L_lam = X1; X1 = X2; X2 = L_lam + gr * (U_lam - L_lam);
      S1 = S2; S2 = nystrom_loo_loss(U_phi, Sigma2, y_tr, X2);
    }
    ++it;
  }
  double lambda = 0.5 * (X1 + X2);

  // Solve beta + alpha
  arma::vec Uty   = U_phi.t() * y_tr;
  arma::vec w     = Sigma2 / (Sigma2 + lambda);
  arma::vec yfit  = U_phi * (w % Uty);
  arma::vec beta_coef = V_phi * ((sigma_phi % Uty) / (Sigma2 + lambda));
  arma::vec alpha = U * (Dinv2 % beta_coef);

  // Optional vcov
  arma::mat vcov_alpha;
  if (compute_vcov) {
    double sigmasq = arma::dot(y_tr - yfit, y_tr - yfit) / (double) n;
    arma::vec beta_weights = sigmasq * Sigma2 / arma::square(Sigma2 + lambda);
    arma::mat Vw = V_phi;
    Vw.each_row() %= beta_weights.t();
    arma::mat vcov_beta = Vw * V_phi.t();
    arma::mat alpha_map = U;
    alpha_map.each_row() %= Dinv2.t();
    vcov_alpha = alpha_map * vcov_beta * alpha_map.t();
  }

  Rcpp::List W_eigen = Rcpp::List::create(
    Rcpp::Named("values")  = D,
    Rcpp::Named("vectors") = U
  );

  return Rcpp::List::create(
    Rcpp::Named("coeffs")           = Rcpp::NumericVector(alpha.begin(),
                                                          alpha.end()),
    Rcpp::Named("fitted_std")       = yfit,
    Rcpp::Named("lambda")           = lambda,
    Rcpp::Named("landmarks")        = Z,
    Rcpp::Named("W_eigen")          = W_eigen,
    Rcpp::Named("Dinvsqrt")         = Dinv2,
    Rcpp::Named("Sigma2")           = Sigma2,
    Rcpp::Named("vcov_alpha")       = vcov_alpha,
    Rcpp::Named("floored_count")    = (int) floored_count,
    Rcpp::Named("D_min_raw")        = D_min_raw,
    Rcpp::Named("D_max_raw")        = Dmax,
    Rcpp::Named("nystrom_m")        = (int) m
  );
}
```

- [ ] **Step 4: Compile + test**

```bash
cd /home/jack/Dropbox/roadrunner && Rscript -e 'devtools::document("."); Rcpp::compileAttributes(".")'
cd /home/jack/Dropbox/roadrunner && Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls-nystrom")'
```

Expected: NYS-FIT-1 PASS at tol=1e-7 vs reference; NYS-FIT-2 PASS.

If NYS-FIT-1 fails: likely L-bracket increment mismatch (reference uses `L + 0.05` linear; C++ uses `L * 10` log-step). Either:
(a) Use the SAME linear `+0.05` step in C++ to match reference exactly.
(b) Loosen tol to 1e-6 if both converge to same lambda within tol.

Investigate before loosening.

- [ ] **Step 5: Commit**

```bash
cd /home/jack/Dropbox/roadrunner && git add src/krls.cpp src/RcppExports.cpp R/RcppExports.R tests/testthat/test-krls-nystrom.R && git commit -m "feat(krls): krls_nystrom_fit_cpp single-fit Nystrom path"
```

---

## Task 4: R harness — route non-autotune `approx="nystrom"` path

**Files:**
- Modify: `R/krls.R` (~80 LOC: add 5 new args + routing)
- Modify: `tests/testthat/test-krls-nystrom.R` (append NYS-DET-1, NYS-DET-2, NYS-PRED-1)

- [ ] **Step 1: Write failing tests**

Append to `tests/testthat/test-krls-nystrom.R`:
```r
test_that("NYS-DET-1: same landmark_seed yields byte-identical fits", {
  set.seed(20L)
  n <- 120L; p <- 4L
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(rowSums(X[, 1:2]) + 0.3 * rnorm(n))
  f1 <- roadrunner::krls(X, y, approx = "nystrom",
                         nystrom_m = 30L, landmark_seed = 42L,
                         derivative = FALSE, vcov = FALSE)
  f2 <- roadrunner::krls(X, y, approx = "nystrom",
                         nystrom_m = 30L, landmark_seed = 42L,
                         derivative = FALSE, vcov = FALSE)
  expect_identical(f1$sigma,  f2$sigma)
  expect_identical(f1$lambda, f2$lambda)
  expect_identical(f1$coeffs, f2$coeffs)
})

test_that("NYS-DET-2: landmark_seed does not disturb caller's RNG state", {
  set.seed(21L)
  n <- 80L; p <- 3L
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  set.seed(999L)
  v_before <- rnorm(5L)
  set.seed(999L)
  fit <- roadrunner::krls(X, y, approx = "nystrom", nystrom_m = 20L,
                          landmark_seed = 42L, derivative = FALSE,
                          vcov = FALSE)
  v_after <- rnorm(5L)

  expect_identical(v_before, v_after)
})

test_that("NYS-PRED-1: predict on a Nystrom fit returns length-n_new numeric", {
  set.seed(22L)
  n <- 100L; p <- 3L
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(rowSums(X) + 0.3 * rnorm(n))
  fit <- roadrunner::krls(X, y, approx = "nystrom", nystrom_m = 25L,
                          landmark_seed = 42L, derivative = FALSE,
                          vcov = FALSE)
  Xn <- matrix(rnorm(15 * p), 15, p)
  yhat <- predict(fit, newdata = Xn)$fit
  expect_length(yhat, 15L)
  expect_true(all(is.finite(yhat)))
})
```

- [ ] **Step 2: Run to verify failure**

```bash
cd /home/jack/Dropbox/roadrunner && Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls-nystrom")'
```

Expected: NYS-DET-1/2 ERROR (`approx` arg not supported), NYS-PRED-1 ERROR.

- [ ] **Step 3: Add five new args to `krls.default` signature**

In `R/krls.R`, locate the `krls.default <- function(` signature. After the last formal arg (before `...` or before the closing paren), add:

```r
                         approx = c("exact", "nystrom"),
                         nystrom_m = NULL,
                         landmarks = NULL,
                         landmark_method = c("random", "kmeans"),
                         landmark_seed = NULL,
                         nystrom_eps = 1e-9,
```

- [ ] **Step 4: Add validation block at top of `krls.default` body**

After the existing `match.arg` calls in `krls.default`, add:
```r
  approx <- match.arg(approx)
  landmark_method <- match.arg(landmark_method)
  nystrom_eps <- .validate_nystrom_eps(nystrom_eps)
```

- [ ] **Step 5: Add Nystrom routing**

Locate where the non-autotune, non-Nystrom fit happens (the existing exact path inside `krls.default`). Add a branch BEFORE that block:

```r
  if (identical(approx, "nystrom") && !isTRUE(autotune)) {
    # ---- Single-fit Nystrom path ----
    lm <- .resolve_landmarks(landmarks, landmark_method, nystrom_m,
                             Xs, X_centers = X_means, X_scales = X_sds,
                             landmark_seed = landmark_seed)
    Z_std <- lm$matrix

    lambda_args <- list(
      tol = 1e-6, L0 = .Machine$double.eps,
      L_step = 10.0, U_start_from_n = TRUE
    )

    if (is.null(sigma)) {
      sigma <- .krls_sigma_anchor(Xs)
    }

    nys <- krls_nystrom_fit_cpp(
      Xs, Z_std, y_std, sigma, lambda_args,
      nystrom_eps, compute_vcov = isTRUE(vcov)
    )

    # Assemble fit object
    fit <- list(
      sigma     = sigma,
      lambda    = nys$lambda,
      coeffs    = nys$coeffs,
      fitted    = nys$fitted_std * y_sd + y_mean,
      X         = X,
      y         = y,
      X_means   = X_means,
      X_sds     = X_sds,
      y_mean    = y_mean,
      y_sd      = y_sd,
      approx    = "nystrom",
      landmarks = Z_std,
      landmark_indices = lm$indices,
      landmark_method  = lm$method_used,
      nystrom_m        = nys$nystrom_m,
      nystrom_eps      = nystrom_eps,
      W_eigen          = nys$W_eigen,
      Dinvsqrt         = nys$Dinvsqrt,
      Sigma2           = nys$Sigma2,
      vcov_alpha       = if (isTRUE(vcov)) nys$vcov_alpha else NULL,
      nystrom_diagnostics = list(
        floored_count = nys$floored_count,
        D_min_raw     = nys$D_min_raw,
        D_max_raw     = nys$D_max_raw
      )
    )
    class(fit) <- "krls_rr"

    # Add R2 like the exact path
    fit$R2 <- 1 - sum((y - fit$fitted)^2) / sum((y - mean(y))^2)

    return(fit)
  }
```

Note: the exact variable names (`Xs`, `y_std`, `X_means`, `X_sds`, `y_mean`, `y_sd`, `.krls_sigma_anchor`) must match what's already in `krls.default`. Verify via `grep -n` before pasting. Adapt names if different.

- [ ] **Step 6: Extend `predict.krls_rr` for the Nystrom path**

Locate `predict.krls_rr` in `R/predict.R`. At the TOP of the function body, after argument validation, add:

```r
  if (identical(object$approx, "nystrom")) {
    # Standardize newdata using training centers / scales
    newdata <- as.matrix(newdata)
    Xn_std <- scale(newdata, center = object$X_means, scale = object$X_sds)
    attr(Xn_std, "scaled:center") <- NULL
    attr(Xn_std, "scaled:scale")  <- NULL
    yhat_std <- as.numeric(krls_nystrom_predict_cpp(
      Xn_std, object$landmarks, as.numeric(object$coeffs), object$sigma
    ))
    yhat <- yhat_std * object$y_sd + object$y_mean
    return(list(fit = yhat))
  }
```

(Adapt to whatever existing return shape `predict.krls_rr` uses. If it returns a numeric vector when `interval=NULL`, return a numeric vector here too. Check existing code.)

- [ ] **Step 7: Compile + test**

```bash
cd /home/jack/Dropbox/roadrunner && Rscript -e 'devtools::document("."); Rcpp::compileAttributes(".")'
cd /home/jack/Dropbox/roadrunner && Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls-nystrom")'
```

Expected: NYS-DET-1, NYS-DET-2, NYS-PRED-1 PASS. Total 10 tests pass in the suite.

- [ ] **Step 8: Smoke fit**

```bash
cd /home/jack/Dropbox/roadrunner && Rscript -e '
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
set.seed(1L)
X <- matrix(rnorm(2000 * 10), 2000, 10)
y <- as.numeric(rowSums(X[, 1:3]) + 0.5 * rnorm(2000))
t_exact <- system.time(fit_e <- roadrunner::krls(X, y, derivative=FALSE,
                                                 vcov=FALSE))[3]
t_nys   <- system.time(fit_n <- roadrunner::krls(X, y, approx="nystrom",
                                                 landmark_seed=42L,
                                                 derivative=FALSE,
                                                 vcov=FALSE))[3]
cat(sprintf("n=2000 p=10: exact=%.2fs nystrom=%.2fs speedup=%.2fx\n",
            t_exact, t_nys, t_exact/t_nys))
cat(sprintf("m used = %d\n", fit_n$nystrom_m))
'
```

Expected: speedup ≥ 3x (informal sanity; formal target is ≥5x in Task 7 sim).

- [ ] **Step 9: Commit**

```bash
cd /home/jack/Dropbox/roadrunner && git add R/krls.R R/predict.R tests/testthat/test-krls-nystrom.R && git commit -m "feat(krls): route non-autotune approx='nystrom' through C++ fit_cpp"
```

(Stage `man/*.Rd` if any regen happened.)

---

## Task 5: `krls_nystrom_autotune_inner_cpp` + autotune harness wiring

**Files:**
- Modify: `src/krls.cpp` (~200 LOC)
- Modify: `R/krls.R` (~50 LOC: autotune branch routes through Nystrom inner)
- Modify: `tests/testthat/test-krls-nystrom.R` (append NYS-AT-1, NYS-AT-EQUIV, NYS-AT-FALL)

- [ ] **Step 1: Write failing tests**

Append:
```r
test_that("NYS-AT-1: krls(approx='nystrom', autotune=TRUE) produces finite fit", {
  set.seed(30L)
  n <- 200L; p <- 5L
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(rowSums(X[, 1:2]) + 0.3 * rnorm(n))
  fit <- roadrunner::krls(X, y, approx = "nystrom", autotune = TRUE,
                          landmark_seed = 42L, seed.cv = 42L,
                          derivative = FALSE, vcov = FALSE,
                          autotune.nthreads = 1L)
  expect_true(is.finite(fit$sigma))
  expect_true(is.finite(fit$lambda))
  expect_identical(fit$approx, "nystrom")
  expect_true(fit$nystrom_m >= 1L)
})

test_that("NYS-AT-EQUIV: autotune Nystrom nthreads=1 vs 4 byte-identical", {
  set.seed(31L)
  n <- 180L; p <- 4L
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(rowSums(X[, 1:2]) + 0.3 * rnorm(n))
  f1 <- roadrunner::krls(X, y, approx="nystrom", autotune=TRUE,
                         landmark_seed=42L, seed.cv=42L,
                         derivative=FALSE, vcov=FALSE,
                         autotune.nthreads=1L)
  f4 <- roadrunner::krls(X, y, approx="nystrom", autotune=TRUE,
                         landmark_seed=42L, seed.cv=42L,
                         derivative=FALSE, vcov=FALSE,
                         autotune.nthreads=4L)
  expect_identical(f1$sigma,  f4$sigma)
  expect_identical(f1$lambda, f4$lambda)
  expect_identical(f1$coeffs, f4$coeffs)
})

test_that("NYS-AT-FALL: autotune Nystrom with weights errors clearly", {
  set.seed(32L)
  n <- 100L; p <- 3L
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  w <- runif(n, 0.5, 2)
  expect_error(
    roadrunner::krls(X, y, weights = w, approx = "nystrom",
                     autotune = TRUE, landmark_seed = 42L),
    "weights.*not supported.*nystrom"
  )
})
```

- [ ] **Step 2: Run to verify failure**

Expected: all 3 ERROR.

- [ ] **Step 3: Implement `krls_nystrom_autotune_inner_cpp` in `src/krls.cpp`**

Append:

```cpp
struct KrlsNystromSigmaWorker : public RcppParallel::Worker {
  const arma::mat& D_tr_Z;
  const arma::mat& D_Z_Z;
  const arma::mat& D_te_Z;
  const arma::vec& y_tr;
  const arma::vec& y_te;
  const arma::vec& sigma_grid;
  double L0, L_step, tol, eps;
  arma::vec& mse;
  arma::vec& lam;

  KrlsNystromSigmaWorker(const arma::mat& D_tr_Z_, const arma::mat& D_Z_Z_,
                         const arma::mat& D_te_Z_,
                         const arma::vec& y_tr_, const arma::vec& y_te_,
                         const arma::vec& sigma_grid_,
                         double L0_, double L_step_, double tol_, double eps_,
                         arma::vec& mse_, arma::vec& lam_)
    : D_tr_Z(D_tr_Z_), D_Z_Z(D_Z_Z_), D_te_Z(D_te_Z_),
      y_tr(y_tr_), y_te(y_te_), sigma_grid(sigma_grid_),
      L0(L0_), L_step(L_step_), tol(tol_), eps(eps_),
      mse(mse_), lam(lam_) {}

  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      double s = sigma_grid(i);
      arma::mat C    = arma::exp(-D_tr_Z / s);
      arma::mat W    = arma::exp(-D_Z_Z  / s);
      arma::mat K_te = arma::exp(-D_te_Z / s);

      arma::vec D; arma::mat U;
      arma::eig_sym(D, U, W);
      D = arma::clamp(D, 0.0, arma::datum::inf);
      double Dmax = D.max();
      if (Dmax <= 0.0) { mse(i) = std::numeric_limits<double>::max(); lam(i) = 0.0; continue; }
      arma::vec D_reg = arma::clamp(D, eps * Dmax, arma::datum::inf);
      arma::vec Dinv2 = 1.0 / arma::sqrt(D_reg);

      arma::mat Phi = C * U;
      Phi.each_row() %= Dinv2.t();

      arma::mat U_phi, V_phi; arma::vec sigma_phi;
      bool ok = arma::svd(U_phi, sigma_phi, V_phi, Phi);
      if (!ok) { mse(i) = std::numeric_limits<double>::max(); lam(i) = 0.0; continue; }
      arma::vec Sigma2 = arma::square(sigma_phi);

      // Lambda bracket
      double Smax = Sigma2.max();
      double U_lam = Smax;
      int it = 0;
      while (arma::sum(Sigma2 / (Sigma2 + U_lam)) >= 1.0 && it < 200) {
        U_lam *= 2.0; ++it;
      }
      double tgt = Smax / 1000.0;
      arma::uword q_idx = 0;
      double mind = std::abs(Sigma2(0) - tgt);
      for (arma::uword j = 1; j < Sigma2.n_elem; ++j) {
        double dd = std::abs(Sigma2(j) - tgt);
        if (dd < mind) { mind = dd; q_idx = j; }
      }
      double q_R = (double)(q_idx + 1);
      double L_lam = L0;
      while (arma::sum(Sigma2 / (Sigma2 + L_lam)) > q_R && L_lam < Smax) {
        L_lam *= L_step;
      }
      if (!(L_lam < U_lam)) { mse(i) = std::numeric_limits<double>::max(); lam(i) = 0.0; continue; }

      // Golden section
      const double gr = (std::sqrt(5.0) - 1.0) / 2.0;
      double X1 = L_lam + (1.0 - gr) * (U_lam - L_lam);
      double X2 = L_lam + gr * (U_lam - L_lam);
      double S1 = nystrom_loo_loss(U_phi, Sigma2, y_tr, X1);
      double S2 = nystrom_loo_loss(U_phi, Sigma2, y_tr, X2);
      int it2 = 0;
      while (std::abs(S1 - S2) > tol && it2 < 200) {
        if (S1 < S2) {
          U_lam = X2; X2 = X1; X1 = L_lam + (1.0 - gr) * (U_lam - L_lam);
          S2 = S1; S1 = nystrom_loo_loss(U_phi, Sigma2, y_tr, X1);
        } else {
          L_lam = X1; X1 = X2; X2 = L_lam + gr * (U_lam - L_lam);
          S1 = S2; S2 = nystrom_loo_loss(U_phi, Sigma2, y_tr, X2);
        }
        ++it2;
      }
      double lambda = 0.5 * (X1 + X2);

      // alpha + test MSE
      arma::vec Uty   = U_phi.t() * y_tr;
      arma::vec beta_coef = V_phi * ((sigma_phi % Uty) / (Sigma2 + lambda));
      arma::vec alpha = U * (Dinv2 % beta_coef);
      arma::vec yhat_te = K_te * alpha;
      mse(i) = arma::mean(arma::square(y_te - yhat_te));
      lam(i) = lambda;
    }
  }
};

// [[Rcpp::export]]
Rcpp::List krls_nystrom_autotune_inner_cpp(
    const arma::mat& X_tr, const arma::mat& Z, const arma::mat& X_te,
    const arma::vec& y_tr, const arma::vec& y_te,
    const arma::vec& sigma_grid, Rcpp::List lambda_args, double eps,
    int nthreads) {
  const arma::uword nsigma = sigma_grid.n_elem;
  arma::vec mse(nsigma, arma::fill::zeros);
  arma::vec lam(nsigma, arma::fill::zeros);

  arma::mat D_tr_Z = krls_pairwise_sqdist_cpp(X_tr, Z);
  arma::mat D_Z_Z  = krls_pairwise_sqdist_cpp(Z,    Z);
  arma::mat D_te_Z = krls_pairwise_sqdist_cpp(X_te, Z);

  double L0     = Rcpp::as<double>(lambda_args["L0"]);
  double L_step = Rcpp::as<double>(lambda_args["L_step"]);
  double tol    = Rcpp::as<double>(lambda_args["tol"]);

  int n_workers = nthreads;
  if (n_workers < 1) n_workers = 1;
  if ((arma::uword) n_workers > nsigma) n_workers = (int) nsigma;

  KrlsNystromSigmaWorker worker(D_tr_Z, D_Z_Z, D_te_Z, y_tr, y_te, sigma_grid,
                                L0, L_step, tol, eps, mse, lam);
  RcppParallel::parallelFor(0, nsigma, worker, /*grainSize=*/1,
                            /*numThreads=*/n_workers);

  return Rcpp::List::create(
    Rcpp::Named("mse_per_sigma")    = mse,
    Rcpp::Named("lambda_per_sigma") = lam,
    Rcpp::Named("nthreads_used")    = n_workers
  );
}
```

- [ ] **Step 4: Wire autotune harness in `R/krls.R`**

In the autotune block of `krls.default`, after the existing weights-check, add:
```r
  if (identical(approx, "nystrom") && isTRUE(autotune)) {
    if (!is.null(weights)) {
      stop("weights are not supported with approx = 'nystrom'; use approx = 'exact'")
    }
    if (!is.null(L) || !is.null(U) || !is.null(tol)) {
      stop("user-supplied L/U/tol not supported with approx = 'nystrom' (yet)")
    }

    sigma_grid_sorted <- sort(autotune.grid)
    ns <- length(sigma_grid_sorted)
    mse_arr <- array(NA_real_, c(ncross, nfold, ns))
    lam_arr <- array(NA_real_, c(ncross, nfold, ns))

    lambda_args <- list(tol = 1e-6, L0 = .Machine$double.eps,
                        L_step = 10.0, U_start_from_n = TRUE)

    for (cross_i in seq_len(ncross)) {
      fold_seed <- as.integer(cross_seeds[cross_i])
      folds <- .ares_make_folds(n, nfold, fold_seed)

      for (k in seq_len(nfold)) {
        idx_te <- folds[[k]]
        X_tr   <- Xs[-idx_te, , drop = FALSE]
        X_te   <- Xs[ idx_te, , drop = FALSE]
        y_tr_k <- y_std[-idx_te]
        y_te_k <- y_std[ idx_te]
        n_tr   <- nrow(X_tr)

        nystrom_m_eff <- if (is.null(nystrom_m)) ceiling(sqrt(n_tr) * 2) else nystrom_m
        fold_landmark_seed <- if (is.null(landmark_seed)) NULL
                              else as.integer(landmark_seed + cross_i * 1000L + k)
        lm <- .resolve_landmarks(landmarks, landmark_method, nystrom_m_eff,
                                 X_tr, X_means, X_sds, fold_landmark_seed)
        Z_std <- lm$matrix

        inner <- krls_nystrom_autotune_inner_cpp(
          X_tr, Z_std, X_te, y_tr_k, y_te_k,
          sigma_grid_sorted, lambda_args, nystrom_eps,
          nthreads = autotune.nthreads
        )

        mse_arr[cross_i, k, ] <- inner$mse_per_sigma
        lam_arr[cross_i, k, ] <- inner$lambda_per_sigma
      }
    }

    mse_per_sigma <- apply(mse_arr, 3, mean, na.rm = TRUE)
    se_per_sigma  <- apply(mse_arr, 3, sd,   na.rm = TRUE) / sqrt(ncross * nfold)
    min_idx       <- which.min(mse_per_sigma)
    threshold     <- mse_per_sigma[min_idx] + se_per_sigma[min_idx]
    sigma_1se_idx <- max(which(mse_per_sigma <= threshold))
    sigma_chosen  <- sigma_grid_sorted[sigma_1se_idx]

    # Refit at chosen sigma using full Xs
    nystrom_m_eff <- if (is.null(nystrom_m)) ceiling(sqrt(n) * 2) else nystrom_m
    lm_full <- .resolve_landmarks(landmarks, landmark_method, nystrom_m_eff,
                                  Xs, X_means, X_sds, landmark_seed)
    nys <- krls_nystrom_fit_cpp(Xs, lm_full$matrix, y_std, sigma_chosen,
                                lambda_args, nystrom_eps,
                                compute_vcov = isTRUE(vcov))

    fit <- list(
      sigma = sigma_chosen, lambda = nys$lambda,
      coeffs = nys$coeffs,
      fitted = nys$fitted_std * y_sd + y_mean,
      X = X, y = y, X_means = X_means, X_sds = X_sds,
      y_mean = y_mean, y_sd = y_sd,
      approx = "nystrom",
      landmarks = lm_full$matrix,
      landmark_indices = lm_full$indices,
      landmark_method = lm_full$method_used,
      nystrom_m = nys$nystrom_m, nystrom_eps = nystrom_eps,
      W_eigen = nys$W_eigen, Dinvsqrt = nys$Dinvsqrt, Sigma2 = nys$Sigma2,
      vcov_alpha = if (isTRUE(vcov)) nys$vcov_alpha else NULL,
      autotune = list(
        grid = sigma_grid_sorted,
        mse_per_sigma = mse_per_sigma,
        se_mse = se_per_sigma,
        mse_per_fold = mse_arr,
        lam_per_fold = lam_arr,
        sigma_1se = sigma_chosen,
        ncross = ncross, nfold = nfold,
        nthreads_used = inner$nthreads_used
      ),
      nystrom_diagnostics = list(
        floored_count = nys$floored_count,
        D_min_raw     = nys$D_min_raw,
        D_max_raw     = nys$D_max_raw
      )
    )
    class(fit) <- "krls_rr"
    fit$R2 <- 1 - sum((y - fit$fitted)^2) / sum((y - mean(y))^2)
    return(fit)
  }
```

- [ ] **Step 5: Compile + test**

```bash
cd /home/jack/Dropbox/roadrunner && Rscript -e 'devtools::document("."); Rcpp::compileAttributes(".")'
cd /home/jack/Dropbox/roadrunner && Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls-nystrom")'
```

Expected: NYS-AT-1, NYS-AT-EQUIV, NYS-AT-FALL all PASS. 13 tests total.

- [ ] **Step 6: Commit**

```bash
cd /home/jack/Dropbox/roadrunner && git add src/krls.cpp src/RcppExports.cpp R/RcppExports.R R/krls.R tests/testthat/test-krls-nystrom.R && git commit -m "feat(krls): krls_nystrom_autotune_inner_cpp + autotune harness"
```

---

## Task 6: Full regression + R CMD check

- [ ] **Step 1: Full krls test suite**

```bash
cd /home/jack/Dropbox/roadrunner && Rscript -e 'devtools::load_all(".", quiet=TRUE); devtools::test(filter="krls")'
```

Expected: 0 failures. Phase 1 tests (krls-autotune-parallel) still pass.

- [ ] **Step 2: R CMD check**

```bash
cd /home/jack/Dropbox/roadrunner && Rscript -e 'devtools::check(".", document = FALSE, args = "--no-tests")'
```

Expected: 0E/0W/no new notes.

- [ ] **Step 3: If `devtools::check` triggers a doc regen, commit man/**

```bash
cd /home/jack/Dropbox/roadrunner && git status --short
# If only man/*.Rd changed:
cd /home/jack/Dropbox/roadrunner && git add man/ && git commit -m "docs: regen Rd for Nystrom args"
```

---

## Task 7: Docs + version bump + EMP-PHASE2 sim + final ship

**Files:**
- Modify: `R/krls.R` (roxygen `@param` for 5 new args + `@details` + `@note`)
- Modify: `NEWS.md` (v0.0.0.9043 entry)
- Modify: `DESCRIPTION` (Version)
- Modify: `man/krls.Rd` (regen)
- Create: `inst/sims/krls-nystrom-phase2.R`

- [ ] **Step 1: Add roxygen for the 5 new args + details + note**

In `R/krls.R`, in the `@param` block, add (in alphabetical or logical order):

```r
#' @param approx Either `"exact"` (default) or `"nystrom"`. The exact
#'   path runs the original O(n^3) eigendecomposition on the full kernel.
#'   `"nystrom"` replaces it with a low-rank Nystrom approximation
#'   anchored at `nystrom_m` landmarks, giving O(m^3) + O(n m^2) cost.
#'   At default `nystrom_m = ceiling(sqrt(n) * 2)` test RMSE is typically
#'   within +10% of exact on smooth DGPs.
#' @param nystrom_m Integer. Number of landmarks when `approx = "nystrom"`.
#'   Default `NULL` resolves to `ceiling(sqrt(n) * 2)` (e.g. n=2000 -> m=90).
#'   Ignored when `landmarks` is supplied explicitly.
#' @param landmarks NULL (default; auto-draw via `landmark_method`), an
#'   integer vector of row indices into X, or an m x d numeric matrix of
#'   landmark coordinates in the original X scale (auto-standardized
#'   internally using the training centers/scales).
#' @param landmark_method Either `"random"` (default; uniform subsample of
#'   training rows) or `"kmeans"` (Hartigan-Wong centers on standardized X;
#'   not bit-stable across R versions). Ignored when `landmarks` is supplied.
#' @param landmark_seed Optional integer seed for the local landmark draw.
#'   When supplied, uses `.with_seed()` to avoid disturbing the caller's
#'   global RNG state. Required for byte-identical fits across consecutive
#'   `krls(approx="nystrom")` calls.
#' @param nystrom_eps Numeric. Relative ridge floor applied to the
#'   landmark-kernel eigenvalues: `D_reg = max(D, nystrom_eps * max(D))`.
#'   Defaults to `1e-9`. Stabilizes the m x m eigendecomposition when
#'   landmarks are near-collinear; the fit object's
#'   `nystrom_diagnostics$floored_count` reports how many eigenvalues hit
#'   this floor (useful for tuning m).
```

Add a `@details` paragraph:

```r
#' Since v0.0.0.9043 `krls()` supports an opt-in Nystrom low-rank
#' approximation via `approx = "nystrom"`. Replacing the full n x n
#' eigendecomposition with an m x m one (m = `nystrom_m`, default
#' `ceiling(sqrt(n) * 2)`) gives ~5x speedup at n=2000 and ~10x at
#' n=5000 with RMSE typically within +10% of the exact fit on smooth
#' DGPs. Landmarks default to a uniform random subsample of training
#' rows; pass `landmark_method = "kmeans"` for centroidal landmarks.
#' Fits are byte-identical at fixed `landmark_seed` and `seed.cv`.
```

Add a `@note` paragraph:

```r
#' `approx = "nystrom"` is currently incompatible with observation
#' `weights` and with user-supplied lambda-search bracket `L`/`U`/`tol`.
#' Both error with a clear message at fit time. `vcov = TRUE` is
#' supported in the single-fit Nystrom path but the resulting `vcov`
#' is for the m-length dual coefficients, not the n-length kernel
#' coefficients. Phase 3 may extend these.
```

- [ ] **Step 2: Bump DESCRIPTION**

```
Version: 0.0.0.9043
```

- [ ] **Step 3: Prepend NEWS.md entry**

```markdown
# roadrunner 0.0.0.9043

## krls() speedup — Phase 2 (Nystrom low-rank approximation)

* New optional argument `approx = "nystrom"` enables the Nystrom
  low-rank approximation. Replaces O(n^3) eigendecomposition on the
  full kernel with O(m^3) + O(n m^2) where m = `nystrom_m`
  (default `ceiling(sqrt(n) * 2)`, e.g. n=2000 -> m=90).
* Five new optional args: `approx`, `nystrom_m`, `landmarks`,
  `landmark_method`, `landmark_seed`, `nystrom_eps`.
* Default exact path is byte-identical to v0.0.0.9042 — opt-in only.
* New exported helper `get_landmarks(fit)` returns landmark
  coordinates (original or standardized X scale).
* Determinism: at fixed `landmark_seed` (+ `seed.cv` when autotune),
  Nystrom fits are byte-identical across `autotune.nthreads`.

### Empirical wall-clock speedup (EMP-PHASE2, R=5 reps)

| n    | p  | approx   | autotune | wall (s) | speedup vs exact |
|------|----|----------|----------|----------|------------------|
| 2000 | 10 | exact    | FALSE    | (TBD)    | 1.0x             |
| 2000 | 10 | nystrom  | FALSE    | (TBD)    | (TBD)            |
| 5000 | 10 | exact    | FALSE    | (TBD)    | 1.0x             |
| 5000 | 10 | nystrom  | FALSE    | (TBD)    | (TBD)            |

(TBD entries filled by the EMP-PHASE2 simulation.)

### API compatibility

* Zero breaking changes. `krls()` without `approx` defaults to `"exact"`
  and behaves identically to v0.0.0.9042.
* `predict.krls_rr` auto-detects Nystrom fits via `fit$approx` and uses
  the cheap cross-kernel path.

### Out of scope (Phase 3)

* Leverage-score landmark selection.
* Auto-engage Nystrom at large n.
* Sparse kernel truncation.
* Bagging + Nystrom integration.
```

- [ ] **Step 4: Regen man/**

```bash
cd /home/jack/Dropbox/roadrunner && Rscript -e 'devtools::document(".")'
```

- [ ] **Step 5: Write EMP-PHASE2 sim**

`inst/sims/krls-nystrom-phase2.R`:
```r
#!/usr/bin/env Rscript
##
## EMP-PHASE2 - Nystrom speedup verification for v0.0.0.9043.
##
## Mandatory sim rules (this session):
##   1. Smoke test first (1 cell, R=2)
##   2. Intermediate CSV per cell (append after every fit)
##   3. Wall-clock cap 5 min total (abort if projected over)

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})

OUT_CSV     <- "inst/sims/results/krls-nystrom-phase2.csv"
WALL_CAP_S  <- 300
MASTER_SEED <- 20260519L
R_REPS      <- 5L

dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
if (file.exists(OUT_CSV)) file.remove(OUT_CSV)

dgp_additive <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(sin(X[, 1]) + X[, 2]^2 + X[, 3] + 0.5 * rnorm(n))
  list(X = X, y = y)
}

run_one <- function(n, p, approx, rep_seed) {
  d_tr <- dgp_additive(n, p, rep_seed)
  d_te <- dgp_additive(1000L, p, rep_seed + 1L)
  t0 <- proc.time()[["elapsed"]]
  if (approx == "exact") {
    fit <- roadrunner::krls(X = d_tr$X, y = d_tr$y,
                            derivative = FALSE, vcov = FALSE)
  } else {
    fit <- roadrunner::krls(X = d_tr$X, y = d_tr$y,
                            approx = "nystrom", landmark_seed = 42L,
                            derivative = FALSE, vcov = FALSE)
  }
  wall <- proc.time()[["elapsed"]] - t0
  yhat <- as.numeric(predict(fit, newdata = d_te$X)$fit)
  test_RMSE <- sqrt(mean((d_te$y - yhat)^2))
  list(wall = wall, test_RMSE = test_RMSE,
       sigma = fit$sigma, lambda = fit$lambda,
       nystrom_m = if (approx == "nystrom") fit$nystrom_m else NA_integer_)
}

append_row <- function(row) {
  write.table(row, OUT_CSV, sep = ",", row.names = FALSE,
              col.names = !file.exists(OUT_CSV),
              append = file.exists(OUT_CSV), quote = FALSE)
}

## ---------- SMOKE ----------
cat("[smoke] n=500 p=5 approx=nystrom R=2 ...\n")
sm1 <- run_one(500L, 5L, "nystrom", MASTER_SEED + 9999L)
sm2 <- run_one(500L, 5L, "nystrom", MASTER_SEED + 10000L)
stopifnot(is.finite(sm1$wall), is.finite(sm2$wall),
          is.finite(sm1$test_RMSE), is.finite(sm2$test_RMSE))
cat(sprintf("[smoke OK] wall %.2fs/%.2fs  RMSE %.4f/%.4f\n",
            sm1$wall, sm2$wall, sm1$test_RMSE, sm2$test_RMSE))

## ---------- FULL ----------
cells <- expand.grid(
  n = c(2000L, 5000L),
  p = c(10L),
  approx = c("exact", "nystrom"),
  stringsAsFactors = FALSE
)

wall_start <- proc.time()[["elapsed"]]
for (ci in seq_len(nrow(cells))) {
  n <- cells$n[ci]; p <- cells$p[ci]; approx <- cells$approx[ci]
  for (r in seq_len(R_REPS)) {
    elapsed <- proc.time()[["elapsed"]] - wall_start
    if (elapsed > WALL_CAP_S - 30) {
      cat(sprintf("[abort] wall-cap approaching (%.0fs).\n", elapsed))
      stop("EMP-PHASE2: 5-min wall cap hit.")
    }
    rep_seed <- MASTER_SEED + (ci - 1L) * R_REPS * 2L + (r - 1L) * 2L
    res <- run_one(n, p, approx, rep_seed)
    row <- data.frame(n = n, p = p, approx = approx, rep = r,
                      seed = rep_seed, wall_s = res$wall,
                      test_RMSE = res$test_RMSE,
                      sigma = res$sigma, lambda = res$lambda,
                      nystrom_m = res$nystrom_m)
    append_row(row)
    cat(sprintf("[%d/%d %s n=%d rep=%d] wall=%.2fs RMSE=%.4f\n",
                ci, nrow(cells), approx, n, r, res$wall, res$test_RMSE))
  }
}

## ---------- SUMMARY ----------
df <- read.csv(OUT_CSV)
agg <- aggregate(cbind(wall_s, test_RMSE) ~ n + p + approx, data = df,
                 FUN = median)
print(agg)

for (n_val in c(2000L, 5000L)) {
  w_e <- agg$wall_s[agg$n == n_val & agg$approx == "exact"]
  w_n <- agg$wall_s[agg$n == n_val & agg$approx == "nystrom"]
  r_e <- agg$test_RMSE[agg$n == n_val & agg$approx == "exact"]
  r_n <- agg$test_RMSE[agg$n == n_val & agg$approx == "nystrom"]
  cat(sprintf("\nn=%d: speedup=%.2fx  RMSE exact=%.4f  RMSE nystrom=%.4f  RMSE delta=%.1f%%\n",
              n_val, w_e/w_n, r_e, r_n, 100 * (r_n - r_e) / r_e))
}

cat("\nResults at:", OUT_CSV, "\n")
```

- [ ] **Step 6: Run EMP-PHASE2**

```bash
cd /home/jack/Dropbox/roadrunner && Rscript inst/sims/krls-nystrom-phase2.R 2>&1 | tee /tmp/emp-phase2-log.txt
```

Expected: speedup ≥ 5x at n=2000, ≥ 10x at n=5000, RMSE delta ≤ +10% at both.

If acceptance misses, BLOCK and report. Do not adjust thresholds.

- [ ] **Step 7: Update NEWS table with measured numbers**

Replace `(TBD)` in NEWS.md v0.0.0.9043 entry with median wall + speedup from the sim output.

- [ ] **Step 8: Final tests + R CMD check**

```bash
cd /home/jack/Dropbox/roadrunner && Rscript -e 'devtools::test(filter="krls")' 2>&1 | tail -8
cd /home/jack/Dropbox/roadrunner && Rscript -e 'devtools::check(".")' 2>&1 | tail -10
```

Expected: 0 fail, 0E/0W.

- [ ] **Step 9: Commit final**

```bash
cd /home/jack/Dropbox/roadrunner && git add DESCRIPTION NEWS.md R/krls.R man/krls.Rd inst/sims/krls-nystrom-phase2.R inst/sims/results/krls-nystrom-phase2.csv && git commit -m "docs+sim: v0.0.0.9043 Nystrom docs + EMP-PHASE2 verification"
```

- [ ] **Step 10: Push**

```bash
git -C /home/jack/Dropbox/roadrunner pull --ff-only origin main && git -C /home/jack/Dropbox/roadrunner push origin main
```

---

## Definition of done

- All 13 Nystrom unit tests pass.
- EMP-PHASE2 acceptance: speedup ≥ 5x at n=2000, ≥ 10x at n=5000; RMSE within +10% of exact.
- Phase 1 tests still pass (exact path unaffected).
- R CMD check 0E/0W.
- NEWS v0.0.0.9043 entry with measured speedup table.
- Pushed to origin/main.
