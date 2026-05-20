# pLDA — Penalized LDA Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `plda()`, a fourth roadrunner fitter implementing Witten & Tibshirani (2011) penalized Fisher's LDA with L1 and fused-lasso penalties, multi-class support, and built-in CV autotune.

**Architecture:** Two-file engine (`R/plda.R` thin wrapper + `src/plda.cpp` C++17/Armadillo engine), identical layout to `krls`/`ares`. The engine runs a minorize-maximize loop whose inner step is a penalized-prox-then-normalize update; multi-class via deflation. Sequential-correct first, TBB parallelism added last.

**Tech Stack:** Rcpp, RcppArmadillo, RcppParallel (TBB), testthat 3e, roxygen2. Parity reference: the `penalizedLDA` CRAN package (Suggests-only).

**Design spec:** `dev/plans/2026-05-20-plda-design.md`

---

## Algorithm reference (read before Task 6)

Standardized data: center each feature by its global mean, divide by within-class sd
`sd_w[j] = sqrt( sum_g sum_{i in g} (x_ij - mean_gj)^2 / n )` (note divisor `n`, matching `penalizedLDA::wcsd.matrix`).

Let `M` be the G×p matrix of standardized class means and `w_g = n_g / n`. Between-class
scatter times a vector is matrix-free: `Sigma_b %*% beta = t(M) %*% (w * (M %*% beta))` — O(Gp), never form p×p.

Each discriminant vector solves `max beta' Sigma_b beta - P(beta)` s.t. `||beta||_2 <= 1`.
Minorize-maximize: given `beta_old`, set `u = Sigma_b %*% beta_old`, then
`beta_new = prox(u) / ||prox(u)||_2` (zero vector if `prox(u) == 0`). Iterate to `tol`/`maxit`.
- L1: `prox(u) = soft_threshold(u, lambda)`.
- Fused: `prox(u) = soft_threshold( tv1d(u, lambda2), lambda )` — Friedman et al. (2007)
  two-stage property: 1-D total-variation prox, then soft-threshold.

Multi-class: after `beta_1..beta_{k-1}` found, deflate the class-mean matrix
`M_k = M %*% (I - B B')` with `B = [beta_1..beta_{k-1}]` (orthonormal), then repeat. Up to `K <= G-1`.

Prediction: project standardized X onto the K discriminant vectors; project class centroids
likewise; assign nearest centroid (Euclidean); posterior = softmax of `-0.5 * squared distance`.

**Correctness is enforced by parity tests against `penalizedLDA`.** When a parity test fails,
the deflation order and the standardization divisor are the first things to check.

---

## File structure

- Create `src/plda.cpp` — C++ engine (within-class moments, MM solver, soft-threshold, TV/FLSA, deflation, projection).
- Create `R/plda.R` — `plda()` generic, `plda.formula`, `plda.default`, `.plda_cv`, validation.
- Modify `R/predict.R` — append `predict.plda`.
- Modify `R/print.R` — append `print.plda`, `summary.plda`, `plot.plda`.
- Create `tests/testthat/test-plda-engine.R`, `test-plda-api.R`, `test-plda-fused.R`, `test-plda-cv.R`, `test-plda-determinism.R`, `test-plda-parity.R`.
- Modify `DESCRIPTION`, `_pkgdown.yml`, `NEWS.md`, `README.md`, `CLAUDE.md`.
- Auto-generated: `R/RcppExports.R`, `src/RcppExports.cpp`, `man/*.Rd` — never hand-edit.

---

## Milestone 0 — Scaffold

### Task 1: DESCRIPTION

**Files:** Modify `DESCRIPTION`

- [ ] **Step 1: Bump version and add Suggests dependency**

In `DESCRIPTION`, bump `Version:` to `0.0.0.9052`. In the `Suggests:` field add `penalizedLDA` (keep alphabetical-ish ordering with the other suggests).

- [ ] **Step 2: Commit**

```bash
git add DESCRIPTION
git commit -m "chore(plda): bump version, add penalizedLDA to Suggests"
```

### Task 2: Engine + wrapper skeleton

**Files:**
- Create: `src/plda.cpp`
- Create: `R/plda.R`
- Test: `tests/testthat/test-plda-engine.R`

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-plda-engine.R
test_that("plda_wcsd_cpp returns within-class sd with divisor n", {
  set.seed(1)
  x <- matrix(rnorm(40), 20, 2)
  g <- rep(1:2, each = 10)
  got <- roadrunner:::plda_wcsd_cpp(x, g, 2L)
  expect_length(got, 2L)
  ref <- vapply(1:2, function(j) {
    ss <- sum(tapply(seq_len(20), g, function(ix) sum((x[ix, j] - mean(x[ix, j]))^2)))
    sqrt(ss / 20)
  }, numeric(1))
  expect_equal(got, ref, tolerance = 1e-10)
})
```

- [ ] **Step 2: Run test to verify it fails**

Run: `Rscript -e 'devtools::load_all(); devtools::test(filter = "plda-engine")'`
Expected: FAIL — `plda_wcsd_cpp` not found.

- [ ] **Step 3: Create the C++ engine skeleton with `plda_wcsd_cpp`**

```cpp
// src/plda.cpp
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Within-class standard deviation per feature, divisor n (matches penalizedLDA::wcsd.matrix).
// [[Rcpp::export]]
arma::vec plda_wcsd_cpp(const arma::mat& x, const arma::ivec& y, int G) {
  const arma::uword n = x.n_rows, p = x.n_cols;
  arma::vec ss(p, arma::fill::zeros);
  for (int g = 1; g <= G; ++g) {
    arma::uvec idx = arma::find(y == g);
    if (idx.n_elem == 0) continue;
    arma::mat xg = x.rows(idx);
    arma::rowvec mg = arma::mean(xg, 0);
    xg.each_row() -= mg;
    ss += arma::sum(arma::square(xg), 0).t();
  }
  return arma::sqrt(ss / (double) n);
}
```

- [ ] **Step 4: Create the R wrapper skeleton**

```r
# R/plda.R

#' Penalized linear discriminant analysis
#'
#' Penalized Fisher's linear discriminant (Witten & Tibshirani 2011) with L1 or
#' fused-lasso penalties, multi-class support, and built-in cross-validation.
#'
#' @param x Numeric predictor matrix, or a formula for the formula interface.
#' @param y Factor (or coercible) class label of length `nrow(x)`.
#' @param ... Passed to methods.
#' @return An object of S3 class `"plda"`.
#' @export
plda <- function(x, ...) UseMethod("plda")
```

- [ ] **Step 5: Regenerate bindings and run the test**

Run: `Rscript -e 'devtools::document(); devtools::load_all(); devtools::test(filter = "plda-engine")'`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/plda.cpp R/plda.R R/RcppExports.R src/RcppExports.cpp NAMESPACE man/
git commit -m "feat(plda): engine skeleton + within-class sd"
```

---

## Milestone 1 — Sequential L1 engine

### Task 3: Soft-threshold and normalize helpers

**Files:** Modify `src/plda.cpp`; Test `tests/testthat/test-plda-engine.R`

- [ ] **Step 1: Write the failing test**

```r
test_that("plda_softthresh_cpp soft-thresholds and is exported for testing", {
  expect_equal(roadrunner:::plda_softthresh_cpp(c(-3, -0.5, 0.5, 3), 1),
               c(-2, 0, 0, 2), tolerance = 1e-12)
})
```

- [ ] **Step 2: Run test — Expected FAIL** (`plda_softthresh_cpp` not found).

Run: `Rscript -e 'devtools::load_all(); devtools::test(filter = "plda-engine")'`

- [ ] **Step 3: Add helpers to `src/plda.cpp`**

```cpp
// Element-wise soft-threshold: sign(u) * max(|u| - lam, 0).
static inline arma::vec soft_threshold(const arma::vec& u, double lam) {
  return arma::sign(u) % arma::clamp(arma::abs(u) - lam, 0.0, arma::datum::inf);
}

// L2-normalize; returns the zero vector unchanged.
static inline arma::vec normalize_l2(const arma::vec& v) {
  double nv = arma::norm(v, 2);
  return (nv > 0.0) ? arma::vec(v / nv) : v;
}

// [[Rcpp::export]]
arma::vec plda_softthresh_cpp(const arma::vec& u, double lam) {
  return soft_threshold(u, lam);
}
```

- [ ] **Step 4: Run test — Expected PASS.**

Run: `Rscript -e 'devtools::load_all(); devtools::test(filter = "plda-engine")'`

- [ ] **Step 5: Commit**

```bash
git add src/plda.cpp R/RcppExports.R src/RcppExports.cpp
git commit -m "feat(plda): soft-threshold and normalize helpers"
```

### Task 4: L1 MM solver — single discriminant (`plda_fit_cpp`, K=1)

**Files:** Modify `src/plda.cpp`; Test `tests/testthat/test-plda-parity.R`

- [ ] **Step 1: Write the failing parity test**

```r
# tests/testthat/test-plda-parity.R
test_that("plda_fit_cpp K=1 L1 matches penalizedLDA::PenalizedLDA", {
  skip_if_not_installed("penalizedLDA")
  set.seed(42)
  n <- 60; p <- 20
  x <- matrix(rnorm(n * p), n, p)
  y <- rep(1:3, length.out = n)
  x[y == 1, 1:3] <- x[y == 1, 1:3] + 2
  ref <- penalizedLDA::PenalizedLDA(x, y, lambda = 0.1, K = 1, type = "standard")
  got <- roadrunner:::plda_fit_cpp(x, as.integer(y), 3L, 1L, 0.1, 0.0,
                                   0L, 100L, 1e-6)
  # discriminant direction matches up to sign
  d_got <- got$discrim[, 1]; d_ref <- ref$discrim[, 1]
  if (sum((d_got - d_ref)^2) > sum((d_got + d_ref)^2)) d_got <- -d_got
  expect_equal(d_got, d_ref, tolerance = 1e-4)
})
```

- [ ] **Step 2: Run test — Expected FAIL** (`plda_fit_cpp` not found).

Run: `Rscript -e 'devtools::load_all(); devtools::test(filter = "plda-parity")'`

- [ ] **Step 3: Implement `plda_fit_cpp` (K=1 path; penalty arg 0 = L1)**

Add to `src/plda.cpp`. `penalty`: 0 = L1, 1 = fused. `lambda2` unused for L1 here.

```cpp
// Standardize x: subtract global column means, divide by within-class sd.
// Drops constant features by leaving their sd as 1 (column becomes ~0).
struct Standardized { arma::mat xs; arma::rowvec mu; arma::vec sdw; };

static Standardized standardize(const arma::mat& x, const arma::ivec& y, int G) {
  Standardized S;
  S.mu  = arma::mean(x, 0);
  S.sdw = plda_wcsd_cpp(x, y, G);
  S.sdw.transform([](double v) { return v > 1e-12 ? v : 1.0; });
  S.xs = x;
  S.xs.each_row() -= S.mu;
  S.xs.each_row() /= S.sdw.t();
  return S;
}

// G x p matrix of class means of the standardized data, and class weights n_g/n.
static void class_means(const arma::mat& xs, const arma::ivec& y, int G,
                        arma::mat& M, arma::vec& w) {
  const arma::uword n = xs.n_rows, p = xs.n_cols;
  M.set_size(G, p); w.set_size(G);
  for (int g = 1; g <= G; ++g) {
    arma::uvec idx = arma::find(y == g);
    M.row(g - 1) = arma::mean(xs.rows(idx), 0);
    w(g - 1) = (double) idx.n_elem / (double) n;
  }
}

// One penalized discriminant via minorize-maximize. M,w already deflated for k>1.
// prox: 0 = L1 soft-threshold; 1 = fused (handled in Task 8).
static arma::vec mm_discriminant(const arma::mat& M, const arma::vec& w,
                                 double lambda, double lambda2, int penalty,
                                 int maxit, double tol) {
  const arma::uword p = M.n_cols;
  arma::vec beta(p, arma::fill::ones);
  beta = normalize_l2(beta);
  for (int it = 0; it < maxit; ++it) {
    arma::vec u = M.t() * (w % (M * beta));      // Sigma_b %*% beta, matrix-free
    arma::vec prox;
    if (penalty == 0) prox = soft_threshold(u, lambda);
    else              prox = soft_threshold(tv1d(u, lambda2), lambda); // Task 8
    arma::vec beta_new = normalize_l2(prox);
    if (arma::norm(beta_new - beta, 2) < tol) { beta = beta_new; break; }
    beta = beta_new;
  }
  return beta;
}

// [[Rcpp::export]]
List plda_fit_cpp(const arma::mat& x, const arma::ivec& y, int G, int K,
                  double lambda, double lambda2, int penalty,
                  int maxit, double tol) {
  Standardized S = standardize(x, y, G);
  arma::mat M; arma::vec w;
  class_means(S.xs, y, G, M, w);

  arma::mat discrim(S.xs.n_cols, K, arma::fill::zeros);
  arma::mat Md = M;                              // deflated class means
  for (int k = 0; k < K; ++k) {
    arma::vec beta = mm_discriminant(Md, w, lambda, lambda2, penalty, maxit, tol);
    discrim.col(k) = beta;
    if (k + 1 < K && arma::norm(beta, 2) > 0.0)  // deflate: M (I - beta beta')
      Md = Md - (Md * beta) * beta.t();
  }
  return List::create(_["discrim"] = discrim, _["mu"] = S.mu,
                      _["sdw"] = S.sdw, _["cmeans"] = M, _["cw"] = w);
}
```

NOTE: `tv1d` is referenced for the fused path but only defined in Task 8. For now, add a
temporary forward stub so the file compiles: `static arma::vec tv1d(const arma::vec& u, double){ return u; }`
placed above `mm_discriminant`. Task 8 replaces the stub body.

- [ ] **Step 4: Run test — Expected PASS.**

Run: `Rscript -e 'devtools::load_all(); devtools::test(filter = "plda-parity")'`
If the direction is off, recheck the standardization divisor (`n`) and that `u` uses the matrix-free form.

- [ ] **Step 5: Commit**

```bash
git add src/plda.cpp R/RcppExports.R src/RcppExports.cpp
git commit -m "feat(plda): L1 minorize-maximize solver, single discriminant"
```

### Task 5: Multi-class deflation (`plda_fit_cpp`, K>1)

**Files:** Test `tests/testthat/test-plda-parity.R`

- [ ] **Step 1: Write the failing test**

```r
test_that("plda_fit_cpp K=2 L1 matches penalizedLDA discrim matrix", {
  skip_if_not_installed("penalizedLDA")
  set.seed(7)
  n <- 90; p <- 25
  x <- matrix(rnorm(n * p), n, p)
  y <- rep(1:3, length.out = n)
  x[y == 1, 1:4] <- x[y == 1, 1:4] + 1.5
  x[y == 2, 5:8] <- x[y == 2, 5:8] - 1.5
  ref <- penalizedLDA::PenalizedLDA(x, y, lambda = 0.08, K = 2, type = "standard")
  got <- roadrunner:::plda_fit_cpp(x, as.integer(y), 3L, 2L, 0.08, 0.0, 0L, 200L, 1e-7)
  for (k in 1:2) {
    a <- got$discrim[, k]; b <- ref$discrim[, k]
    if (sum((a - b)^2) > sum((a + b)^2)) a <- -a
    expect_equal(a, b, tolerance = 1e-3, info = paste("discriminant", k))
  }
})
```

- [ ] **Step 2: Run test**

Run: `Rscript -e 'devtools::load_all(); devtools::test(filter = "plda-parity")'`
Expected: PASS if Task 4's deflation is correct. If discriminant 2 mismatches, the deflation
update `Md = Md - (Md*beta)*beta.t()` is the place to fix — compare against the deflation in
`penalizedLDA:::PenalizedLDA` source.

- [ ] **Step 3: Commit (only if a deflation fix was needed)**

```bash
git add src/plda.cpp R/RcppExports.R src/RcppExports.cpp
git commit -m "test(plda): multi-class deflation parity"
```

### Task 6: `plda.default` R wrapper

**Files:** Modify `R/plda.R`; Test `tests/testthat/test-plda-api.R`

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-plda-api.R
test_that("plda.default returns a well-formed plda object", {
  set.seed(3)
  x <- as.matrix(iris[, 1:4]); y <- iris$Species
  fit <- plda(x, y, K = 2, lambda = 0.1, autotune = FALSE)
  expect_s3_class(fit, "plda")
  expect_equal(dim(fit$discrim), c(4L, 2L))
  expect_equal(fit$classes, levels(y))
  expect_equal(fit$penalty, "L1")
})

test_that("plda rejects K greater than G-1", {
  x <- as.matrix(iris[, 1:4]); y <- iris$Species
  expect_error(plda(x, y, K = 3, lambda = 0.1, autotune = FALSE), "K")
})
```

- [ ] **Step 2: Run test — Expected FAIL** (no `plda.default`).

Run: `Rscript -e 'devtools::load_all(); devtools::test(filter = "plda-api")'`

- [ ] **Step 3: Implement `plda.default` in `R/plda.R`**

```r
#' @param K Number of discriminant vectors (`<= G-1`). Default `G-1`.
#' @param lambda Magnitude penalty. Required when `autotune = FALSE`.
#' @param penalty `"L1"` (default) or `"fused"`.
#' @param lambda2 Fused-lasso difference penalty (used when `penalty = "fused"`).
#' @param autotune If `TRUE` (default), cross-validate `lambda` (and `K`).
#' @param nfold CV folds. @param lambda_grid Optional CV grid.
#' @param maxit,tol MM solver controls.
#' @rdname plda
#' @export
plda.default <- function(x, y, K = NULL, lambda = NULL,
                         penalty = c("L1", "fused"), lambda2 = NULL,
                         autotune = TRUE, nfold = 5L, lambda_grid = NULL,
                         maxit = 100L, tol = 1e-6, ...) {
  penalty <- match.arg(penalty)
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("plda: `x` must be numeric.", call. = FALSE)
  y <- as.factor(y)
  if (nrow(x) != length(y)) stop("plda: nrow(x) must equal length(y).", call. = FALSE)
  classes <- levels(y)
  G <- length(classes)
  if (G < 2L) stop("plda: need at least two classes.", call. = FALSE)
  yint <- as.integer(y)
  if (is.null(K)) K <- G - 1L
  if (K < 1L || K > G - 1L)
    stop(sprintf("plda: K must be in 1..%d (G-1).", G - 1L), call. = FALSE)
  lam2 <- if (is.null(lambda2)) 0 else lambda2
  pen_code <- if (penalty == "fused") 1L else 0L

  if (autotune) {
    cv <- .plda_cv(x, yint, G, K, penalty, pen_code, lam2, nfold,
                   lambda_grid, maxit, tol)
    lambda <- cv$lambda; K <- cv$K
  } else if (is.null(lambda)) {
    stop("plda: supply `lambda` or use autotune = TRUE.", call. = FALSE)
  } else {
    cv <- NULL
  }

  eng <- plda_fit_cpp(x, yint, G, K, lambda, lam2, pen_code,
                      as.integer(maxit), tol)
  structure(list(discrim = eng$discrim, mu = eng$mu, sdw = eng$sdw,
                 cmeans = eng$cmeans, cw = eng$cw, classes = classes,
                 K = K, lambda = lambda, lambda2 = lam2, penalty = penalty,
                 cv = cv, call = match.call()),
            class = "plda")
}
```

For this task only, add a temporary `.plda_cv` stub so `autotune = TRUE` does not error
(Task 11 replaces it): `.plda_cv <- function(...) list(lambda = 0.1, K = 1L, grid = NULL)`.

- [ ] **Step 4: Run test — Expected PASS.**

Run: `Rscript -e 'devtools::document(); devtools::load_all(); devtools::test(filter = "plda-api")'`

- [ ] **Step 5: Commit**

```bash
git add R/plda.R NAMESPACE man/
git commit -m "feat(plda): plda.default matrix interface"
```

### Task 7: `plda.formula` interface

**Files:** Modify `R/plda.R`; Test `tests/testthat/test-plda-api.R`

- [ ] **Step 1: Write the failing test**

```r
test_that("plda.formula equals plda.default", {
  f <- plda(Species ~ ., data = iris, K = 2, lambda = 0.1, autotune = FALSE)
  d <- plda(as.matrix(iris[, 1:4]), iris$Species, K = 2, lambda = 0.1, autotune = FALSE)
  expect_equal(unname(f$discrim), unname(d$discrim), tolerance = 1e-10)
  expect_s3_class(f, "plda")
})
```

- [ ] **Step 2: Run test — Expected FAIL.**

Run: `Rscript -e 'devtools::load_all(); devtools::test(filter = "plda-api")'`

- [ ] **Step 3: Implement `plda.formula`**

```r
#' @rdname plda
#' @export
plda.formula <- function(formula, data, ...) {
  mf <- model.frame(formula, data)
  y <- model.response(mf)
  x <- model.matrix(attr(mf, "terms"), mf)
  icpt <- match("(Intercept)", colnames(x), nomatch = 0L)
  if (icpt > 0L) x <- x[, -icpt, drop = FALSE]
  fit <- plda.default(x, y, ...)
  fit$call <- match.call()
  fit$terms <- attr(mf, "terms")
  fit$xlevels <- .getXlevels(attr(mf, "terms"), mf)
  fit
}
```

- [ ] **Step 4: Run test — Expected PASS.**

Run: `Rscript -e 'devtools::document(); devtools::load_all(); devtools::test(filter = "plda-api")'`

- [ ] **Step 5: Commit**

```bash
git add R/plda.R NAMESPACE man/
git commit -m "feat(plda): formula interface"
```

### Task 8: `predict.plda`

**Files:** Modify `src/plda.cpp`, `R/predict.R`; Test `tests/testthat/test-plda-api.R`

- [ ] **Step 1: Write the failing test**

```r
test_that("predict.plda returns class / posterior / projection", {
  set.seed(11)
  x <- as.matrix(iris[, 1:4]); y <- iris$Species
  fit <- plda(x, y, K = 2, lambda = 0.05, autotune = FALSE)
  cl <- predict(fit, x)
  expect_s3_class(cl, "factor")
  expect_equal(levels(cl), levels(y))
  expect_gt(mean(cl == y), 0.9)
  po <- predict(fit, x, type = "posterior")
  expect_equal(dim(po), c(150L, 3L))
  expect_equal(unname(rowSums(po)), rep(1, 150), tolerance = 1e-8)
  pr <- predict(fit, x, type = "projection")
  expect_equal(dim(pr), c(150L, 2L))
})
```

- [ ] **Step 2: Run test — Expected FAIL** (no `predict.plda`).

Run: `Rscript -e 'devtools::load_all(); devtools::test(filter = "plda-api")'`

- [ ] **Step 3: Add the C++ projection export**

```cpp
// Project new data onto stored discriminant vectors after standardizing
// with the training mu/sdw. Returns an m x K score matrix.
// [[Rcpp::export]]
arma::mat plda_project_cpp(const arma::mat& xnew, const arma::rowvec& mu,
                           const arma::vec& sdw, const arma::mat& discrim) {
  arma::mat xs = xnew;
  xs.each_row() -= mu;
  xs.each_row() /= sdw.t();
  return xs * discrim;
}
```

- [ ] **Step 4: Implement `predict.plda` in `R/predict.R`**

```r
#' Predict from a penalized LDA fit
#'
#' @param object A `"plda"` object.
#' @param newdata Numeric matrix or data.frame of predictors.
#' @param type `"class"` (default), `"posterior"`, or `"projection"`.
#' @param ... Unused.
#' @return Factor of class labels, posterior probability matrix, or projection matrix.
#' @export
predict.plda <- function(object, newdata, type = c("class", "posterior", "projection"), ...) {
  type <- match.arg(type)
  if (!is.null(object$terms)) {
    Terms <- delete.response(object$terms)
    mf <- model.frame(Terms, as.data.frame(newdata), xlev = object$xlevels)
    x <- model.matrix(Terms, mf)
    icpt <- match("(Intercept)", colnames(x), nomatch = 0L)
    if (icpt > 0L) x <- x[, -icpt, drop = FALSE]
  } else {
    x <- as.matrix(newdata)
  }
  scores  <- plda_project_cpp(x, object$mu, object$sdw, object$discrim)
  cscores <- object$cmeans %*% object$discrim          # G x K centroid projections
  if (type == "projection") return(scores)
  # squared distance from each row to each class centroid in projected space
  d2 <- outer(rowSums(scores^2), rowSums(cscores^2), `+`) - 2 * scores %*% t(cscores)
  if (type == "posterior") {
    m <- d2 - apply(d2, 1L, min)
    e <- exp(-0.5 * m)
    p <- e / rowSums(e)
    colnames(p) <- object$classes
    return(p)
  }
  factor(object$classes[max.col(-d2, ties.method = "first")], levels = object$classes)
}
```

- [ ] **Step 5: Run test — Expected PASS.**

Run: `Rscript -e 'devtools::document(); devtools::load_all(); devtools::test(filter = "plda-api")'`

- [ ] **Step 6: Commit**

```bash
git add src/plda.cpp R/predict.R R/RcppExports.R src/RcppExports.cpp NAMESPACE man/
git commit -m "feat(plda): predict method (class/posterior/projection)"
```

### Task 9: `print` / `summary` / `plot`

**Files:** Modify `R/print.R`; Test `tests/testthat/test-plda-api.R`

- [ ] **Step 1: Write the failing test**

```r
test_that("plda print and summary work", {
  fit <- plda(Species ~ ., data = iris, K = 2, lambda = 0.1, autotune = FALSE)
  expect_output(print(fit), "Penalized LDA")
  s <- summary(fit)
  expect_s3_class(s, "summary.plda")
  expect_output(print(s), "nonzero")
})
```

- [ ] **Step 2: Run test — Expected FAIL.**

Run: `Rscript -e 'devtools::load_all(); devtools::test(filter = "plda-api")'`

- [ ] **Step 3: Implement the three S3 methods in `R/print.R`**

```r
#' @export
print.plda <- function(x, ...) {
  cat("Penalized LDA (", x$penalty, " penalty)\n", sep = "")
  cat("  Classes:      ", paste(x$classes, collapse = ", "), "\n", sep = "")
  cat("  Discriminants:", x$K, "\n")
  cat("  lambda:       ", format(x$lambda, digits = 4), "\n", sep = "")
  if (x$penalty == "fused")
    cat("  lambda2:      ", format(x$lambda2, digits = 4), "\n", sep = "")
  nz <- colSums(abs(x$discrim) > 0)
  cat("  Nonzero feats per discriminant:", paste(nz, collapse = ", "), "\n")
  invisible(x)
}

#' @export
summary.plda <- function(object, ...) {
  nz <- colSums(abs(object$discrim) > 0)
  structure(list(penalty = object$penalty, classes = object$classes,
                 K = object$K, lambda = object$lambda, lambda2 = object$lambda2,
                 nonzero = nz, npred = nrow(object$discrim),
                 cv = object$cv), class = "summary.plda")
}

#' @export
print.summary.plda <- function(x, ...) {
  cat("Penalized LDA summary\n")
  cat("  penalty:", x$penalty, "| classes:", length(x$classes),
      "| discriminants:", x$K, "\n")
  cat("  lambda:", format(x$lambda, digits = 4), "\n")
  for (k in seq_len(x$K))
    cat(sprintf("  discriminant %d: %d / %d nonzero features\n",
                k, x$nonzero[k], x$npred))
  if (!is.null(x$cv))
    cat("  CV-selected via", length(x$cv$grid), "lambda grid points\n")
  invisible(x)
}

#' @export
plot.plda <- function(x, data = NULL, labels = NULL, ...) {
  if (is.null(data)) stop("plot.plda: supply `data` (predictor matrix used to fit).",
                          call. = FALSE)
  sc <- plda_project_cpp(as.matrix(data), x$mu, x$sdw, x$discrim)
  cols <- if (is.null(labels)) 1L else as.integer(as.factor(labels))
  if (x$K >= 2L) {
    plot(sc[, 1], sc[, 2], col = cols, pch = 19,
         xlab = "Discriminant 1", ylab = "Discriminant 2",
         main = "Penalized LDA projection", ...)
  } else {
    d <- sc[, 1]
    plot(d, jitter(rep(0, length(d))), col = cols, pch = 19, yaxt = "n",
         xlab = "Discriminant 1", ylab = "", main = "Penalized LDA projection", ...)
  }
  invisible(x)
}
```

- [ ] **Step 4: Run test — Expected PASS.**

Run: `Rscript -e 'devtools::document(); devtools::load_all(); devtools::test(filter = "plda-api")'`

- [ ] **Step 5: Commit**

```bash
git add R/print.R NAMESPACE man/
git commit -m "feat(plda): print, summary, plot methods"
```

---

## Milestone 2 — Fused-lasso penalty

### Task 10: 1-D total-variation prox (`tv1d`)

**Files:** Modify `src/plda.cpp`; Test `tests/testthat/test-plda-fused.R`

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-plda-fused.R
test_that("plda_tv1d_cpp solves the 1-D fused-lasso signal approximator", {
  # prox of lambda2 * TV: minimize 0.5||u-b||^2 + lam2 * sum|b_j - b_{j-1}|
  u <- c(1, 1, 1, 5, 5, 5)
  got <- roadrunner:::plda_tv1d_cpp(u, 0.0)
  expect_equal(got, u, tolerance = 1e-10)            # lam2 = 0 is identity
  big <- roadrunner:::plda_tv1d_cpp(u, 100)
  expect_equal(big, rep(mean(u), 6), tolerance = 1e-6)  # huge lam2 -> constant = mean
})
```

- [ ] **Step 2: Run test — Expected FAIL.**

Run: `Rscript -e 'devtools::load_all(); devtools::test(filter = "plda-fused")'`

- [ ] **Step 3: Replace the `tv1d` stub with Condat's (2013) direct 1-D TV algorithm**

Reference: L. Condat, "A direct algorithm for 1-D total variation denoising," IEEE SPL 2013.
Replace the temporary stub from Task 4 with the full algorithm and add the export.

```cpp
// Condat (2013) direct 1-D total-variation denoising:
// minimize 0.5 * sum (b_j - u_j)^2 + lam * sum |b_{j+1} - b_j|.
static arma::vec tv1d(const arma::vec& u, double lam) {
  const int N = (int) u.n_elem;
  arma::vec x(N);
  if (N == 0) return x;
  if (lam <= 0.0) { x = u; return x; }
  int k = 0, k0 = 0, kplus = 0, kminus = 0;
  double vmin = u[0] - lam, vmax = u[0] + lam;
  double umin = lam, umax = -lam;
  while (true) {
    if (k == N - 1) {
      x[k] = vmin + umin;
      break;
    }
    if (u[k + 1] + umin < vmin - lam) {
      for (int i = k0; i <= kminus; ++i) x[i] = vmin;
      k = k0 = kminus = kplus = kminus + 1;
      vmin = u[k]; vmax = u[k] + 2 * lam;
      umin = lam; umax = -lam;
    } else if (u[k + 1] + umax > vmax + lam) {
      for (int i = k0; i <= kplus; ++i) x[i] = vmax;
      k = k0 = kminus = kplus = kplus + 1;
      vmin = u[k] - 2 * lam; vmax = u[k];
      umin = lam; umax = -lam;
    } else {
      ++k;
      umin += u[k] - vmin;
      umax += u[k] - vmax;
      if (umin >= lam) { vmin += (umin - lam) / (k - k0 + 1); umin = lam; kminus = k; }
      if (umax <= -lam) { vmax += (umax + lam) / (k - k0 + 1); umax = -lam; kplus = k; }
    }
    if (k == N - 1) {
      if (umin < 0.0) {
        for (int i = k0; i <= kminus; ++i) x[i] = vmin;
        k = k0 = kminus = kminus + 1;
        vmin = u[k]; umin = lam; umax = u[k] + lam - vmax;
      } else if (umax > 0.0) {
        for (int i = k0; i <= kplus; ++i) x[i] = vmax;
        k = k0 = kplus = kplus + 1;
        vmax = u[k]; umax = -lam; umin = u[k] - lam - vmin;
      } else {
        for (int i = k0; i <= N - 1; ++i) x[i] = vmin + umin / (k - k0 + 1);
        break;
      }
    }
  }
  return x;
}

// [[Rcpp::export]]
arma::vec plda_tv1d_cpp(const arma::vec& u, double lam) { return tv1d(u, lam); }
```

- [ ] **Step 4: Run test — Expected PASS.**

Run: `Rscript -e 'devtools::load_all(); devtools::test(filter = "plda-fused")'`
If the constant-limit case fails, Condat's edge handling is the suspect — cross-check the
endpoint branch against the published pseudocode.

- [ ] **Step 5: Commit**

```bash
git add src/plda.cpp R/RcppExports.R src/RcppExports.cpp
git commit -m "feat(plda): Condat 1-D total-variation prox"
```

### Task 11: Wire fused penalty through `plda()`

**Files:** Test `tests/testthat/test-plda-fused.R`

- [ ] **Step 1: Write the failing parity test**

```r
test_that("plda fused penalty matches penalizedLDA type='ordered'", {
  skip_if_not_installed("penalizedLDA")
  set.seed(99)
  n <- 80; p <- 30
  x <- matrix(rnorm(n * p), n, p)
  y <- rep(1:2, length.out = n)
  x[y == 1, 10:20] <- x[y == 1, 10:20] + 1.2
  ref <- penalizedLDA::PenalizedLDA(x, y, lambda = 0.05, K = 1,
                                    type = "ordered", lambda2 = 0.1)
  fit <- plda(x, factor(y), K = 1, lambda = 0.05, penalty = "fused",
              lambda2 = 0.1, autotune = FALSE)
  a <- fit$discrim[, 1]; b <- ref$discrim[, 1]
  if (sum((a - b)^2) > sum((a + b)^2)) a <- -a
  expect_equal(a, b, tolerance = 1e-3)
})
```

- [ ] **Step 2: Run test**

Run: `Rscript -e 'devtools::load_all(); devtools::test(filter = "plda-fused")'`
Expected: PASS — the fused path (`pen_code = 1`) already routes through
`soft_threshold(tv1d(u, lambda2), lambda)` in `mm_discriminant` (Task 4), and `tv1d` is now
real (Task 10). If it fails, confirm `plda.default` forwards `lambda2` and `penalty = "fused"`.

- [ ] **Step 3: Commit (only if a wiring fix was needed)**

```bash
git add R/plda.R src/plda.cpp R/RcppExports.R src/RcppExports.cpp
git commit -m "test(plda): fused-lasso parity"
```

---

## Milestone 3 — CV autotune

### Task 12: `.plda_cv` cross-validation

**Files:** Modify `R/plda.R`; Test `tests/testthat/test-plda-cv.R`

- [ ] **Step 1: Write the failing test**

```r
# tests/testthat/test-plda-cv.R
test_that(".plda_cv selects a lambda and reports a grid", {
  set.seed(5)
  x <- as.matrix(iris[, 1:4]); yint <- as.integer(iris$Species)
  cv <- roadrunner:::.plda_cv(x, yint, G = 3L, K = 2L, penalty = "L1",
                              pen_code = 0L, lam2 = 0, nfold = 5L,
                              lambda_grid = NULL, maxit = 100L, tol = 1e-6)
  expect_true(is.numeric(cv$lambda) && cv$lambda > 0)
  expect_true(cv$K >= 1L && cv$K <= 2L)
  expect_true(length(cv$grid) >= 5L)
  expect_equal(length(cv$errors), length(cv$grid))
})

test_that("plda autotune fit predicts iris well", {
  set.seed(6)
  fit <- plda(Species ~ ., data = iris)   # autotune = TRUE default
  acc <- mean(predict(fit, iris) == iris$Species)
  expect_gt(acc, 0.9)
})
```

- [ ] **Step 2: Run test — Expected FAIL** (stub `.plda_cv` from Task 6 has no `$grid`).

Run: `Rscript -e 'devtools::load_all(); devtools::test(filter = "plda-cv")'`

- [ ] **Step 3: Replace the `.plda_cv` stub in `R/plda.R`**

```r
# Default data-driven log-spaced lambda grid.
.plda_lambda_grid <- function(x, yint, G, n = 12L) {
  S <- plda_wcsd_cpp(x, yint, G)
  S[S < 1e-12] <- 1
  xs <- scale(x, center = TRUE, scale = FALSE) %*% diag(1 / S)
  hi <- max(abs(colMeans(xs))) + 1e-8
  exp(seq(log(hi * 1e-3), log(hi), length.out = n))
}

# k-fold CV over (lambda grid) x (1..K); picks the lambda/K with lowest mean
# misclassification error. Folds are seeded for determinism.
.plda_cv <- function(x, yint, G, K, penalty, pen_code, lam2, nfold,
                     lambda_grid, maxit, tol) {
  if (is.null(lambda_grid)) lambda_grid <- .plda_lambda_grid(x, yint, G)
  n <- nrow(x)
  set.seed(0L)
  folds <- sample(rep_len(seq_len(nfold), n))
  err <- matrix(0, length(lambda_grid), K)
  for (f in seq_len(nfold)) {
    tr <- folds != f; te <- !tr
    for (li in seq_along(lambda_grid)) {
      eng <- plda_fit_cpp(x[tr, , drop = FALSE], yint[tr], G, K,
                          lambda_grid[li], lam2, pen_code,
                          as.integer(maxit), tol)
      sc  <- plda_project_cpp(x[te, , drop = FALSE], eng$mu, eng$sdw, eng$discrim)
      csc <- eng$cmeans %*% eng$discrim
      for (k in seq_len(K)) {
        d2 <- outer(rowSums(sc[, 1:k, drop = FALSE]^2),
                    rowSums(csc[, 1:k, drop = FALSE]^2), `+`) -
              2 * sc[, 1:k, drop = FALSE] %*% t(csc[, 1:k, drop = FALSE])
        pred <- max.col(-d2, ties.method = "first")
        err[li, k] <- err[li, k] + sum(pred != yint[te])
      }
    }
  }
  err <- err / n
  best <- which(err == min(err), arr.ind = TRUE)[1, ]
  list(lambda = lambda_grid[best[1]], K = as.integer(best[2]),
       grid = lambda_grid, errors = err[, best[2]])
}
```

- [ ] **Step 4: Run test — Expected PASS.**

Run: `Rscript -e 'devtools::document(); devtools::load_all(); devtools::test(filter = "plda-cv")'`

- [ ] **Step 5: Commit**

```bash
git add R/plda.R man/
git commit -m "feat(plda): k-fold CV autotune over lambda and K"
```

---

## Milestone 4 — Parallelism, docs, CI

### Task 13: Parallelize the CV fold loop (TBB) + determinism test

**Files:** Modify `src/plda.cpp`, `R/plda.R`; Test `tests/testthat/test-plda-determinism.R`

- [ ] **Step 1: Write the determinism test**

```r
# tests/testthat/test-plda-determinism.R
test_that("plda fit is byte-identical across thread counts", {
  set.seed(8)
  x <- matrix(rnorm(120 * 15), 120, 15)
  y <- factor(rep(1:3, length.out = 120))
  f1 <- plda(x, y, K = 2, lambda = 0.07, autotune = FALSE, nthreads = 1L)
  f4 <- plda(x, y, K = 2, lambda = 0.07, autotune = FALSE, nthreads = 4L)
  expect_identical(f1$discrim, f4$discrim)
})

test_that("plda autotune is byte-identical across thread counts", {
  set.seed(9)
  x <- matrix(rnorm(150 * 12), 150, 12)
  y <- factor(rep(1:3, length.out = 150))
  f1 <- plda(x, y, nthreads = 1L)
  f4 <- plda(x, y, nthreads = 4L)
  expect_identical(f1$discrim, f4$discrim)
  expect_identical(f1$lambda, f4$lambda)
})
```

- [ ] **Step 2: Run test — Expected FAIL** (`nthreads` arg not yet accepted).

Run: `Rscript -e 'devtools::load_all(); devtools::test(filter = "plda-determinism")'`

- [ ] **Step 3: Add an `nthreads` argument and parallelize CV folds**

Add `nthreads = getOption("roadrunner.threads", 1L)` to the `plda.default` signature and pass
it to `.plda_cv`. Inside `.plda_cv`, parallelize the **fold loop** with `RcppParallel` (or
`parallel::mclapply` guarded for Windows) so each fold computes into its own `err` matrix slice,
then reduce by summation in fixed fold order. Because each fold writes a disjoint slice and the
reduction order is fixed, the result is byte-identical regardless of `nthreads`. Follow the TBB
pattern already used in `src/krls.cpp` (`krls_autotune_inner_cpp`): a `Worker` struct over the
fold index, fixed-slot writes, deterministic serial reduction.

The single fit (`plda_fit_cpp`) MM loop stays sequential; only per-feature sweeps
(`plda_wcsd_cpp`, `soft_threshold`) may be parallelized across columns with fixed-index writes.

- [ ] **Step 4: Run test — Expected PASS.**

Run: `Rscript -e 'devtools::document(); devtools::load_all(); devtools::test(filter = "plda-determinism")'`

- [ ] **Step 5: Commit**

```bash
git add R/plda.R src/plda.cpp R/RcppExports.R src/RcppExports.cpp man/
git commit -m "perf(plda): parallel CV folds, determinism preserved"
```

NOTE: if parallelization risks breaking determinism or schedule, ship Task 13 as the
**sequential** path (skip Step 3 parallelization, keep `nthreads` accepted-but-ignored) and the
determinism tests still pass trivially. Parallelism can be a follow-up commit. Correctness and
determinism are never traded.

### Task 14: Adversarial / edge-case tests

**Files:** Test `tests/testthat/test-plda-api.R`

- [ ] **Step 1: Write the tests**

```r
test_that("plda handles edge cases", {
  set.seed(12)
  # p >> n
  x <- matrix(rnorm(20 * 200), 20, 200); y <- factor(rep(1:2, each = 10))
  expect_s3_class(plda(x, y, K = 1, lambda = 0.2, autotune = FALSE), "plda")
  # huge lambda -> all-zero discriminant
  f0 <- plda(x, y, K = 1, lambda = 1e6, autotune = FALSE)
  expect_true(all(f0$discrim[, 1] == 0))
  # constant feature does not crash
  xc <- x; xc[, 1] <- 3
  expect_s3_class(plda(xc, y, K = 1, lambda = 0.2, autotune = FALSE), "plda")
  # OOV factor level in predict errors clearly
  fit <- plda(Species ~ ., data = iris, K = 2, lambda = 0.1, autotune = FALSE)
  bad <- iris[1:5, ]
  expect_error(predict(fit, bad[, 1:3]))   # missing predictor column
})
```

- [ ] **Step 2: Run test**

Run: `Rscript -e 'devtools::load_all(); devtools::test(filter = "plda-api")'`
Expected: PASS. If the all-zero case fails, confirm `normalize_l2` returns the zero vector
unchanged when `prox(u)` is zero.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-plda-api.R
git commit -m "test(plda): adversarial and edge-case coverage"
```

### Task 15: roxygen docs + pkgdown + NEWS + README + CLAUDE.md

**Files:** Modify `_pkgdown.yml`, `NEWS.md`, `README.md`, `CLAUDE.md`

- [ ] **Step 1: Add a pkgdown reference section**

In `_pkgdown.yml`, add immediately before the `- title: "Package"` section:

```yaml
- title: "Penalized LDA: plda()"
  desc: "Penalized Fisher's linear discriminant (Witten & Tibshirani 2011) with L1 and fused-lasso penalties, multi-class, and CV autotune."
  contents:
  - plda
  - predict.plda
  - print.plda
  - summary.plda
  - plot.plda
```

- [ ] **Step 2: Verify the pkgdown reference index locally**

Run: `Rscript -e 'pkg <- pkgdown::as_pkgdown("."); invisible(pkgdown:::data_reference_index(pkg)); cat("REFERENCE INDEX OK\n")'`
Expected: `REFERENCE INDEX OK`. If a topic is reported missing, every exported `plda` symbol
must appear either in `_pkgdown.yml` or carry `@keywords internal`.

- [ ] **Step 3: Update `NEWS.md`**

Add a top entry: `# roadrunner 0.0.0.9052` with a bullet describing `plda()` — penalized LDA,
L1 + fused-lasso penalties, multi-class, CV autotune.

- [ ] **Step 4: Update `README.md`**

Add a `### Penalized LDA via plda()` usage block alongside the `ares()`/`krls()` examples:

```r
fit <- plda(Species ~ ., data = iris)
predict(fit, iris)
```

Update the intro line to say the package ships three algorithms plus `meep()` — now four core
fitters (`ares`, `krls`, `plda`) plus the `meep()` ensemble.

- [ ] **Step 5: Update `CLAUDE.md`**

Add `plda()` to the "Three algorithms shipped" list (now four), add an `R/plda.R` /
`src/plda.cpp` row to the "Where things live" tables, and add a `## plda()` section
mirroring the `krls()` section (algorithm, penalties, tuning, S3 methods).

- [ ] **Step 6: Regenerate documentation**

Run: `Rscript -e 'devtools::document()'`
Expected: `man/plda.Rd`, `man/predict.plda.Rd`, etc. generated; `NAMESPACE` exports `plda`,
`predict.plda`, `print.plda`, `summary.plda`, `plot.plda`.

- [ ] **Step 7: Commit**

```bash
git add _pkgdown.yml NEWS.md README.md CLAUDE.md man/ NAMESPACE
git commit -m "docs(plda): reference index, NEWS, README, handoff brief"
```

### Task 16: Full check, lint, push, CI green

**Files:** none (verification only)

- [ ] **Step 1: Lint**

Run: `Rscript -e 'print(lintr::lint_package())'`
Expected: no lints. Fix any in `R/plda.R` / `R/predict.R` / `R/print.R` (infix spaces, line
length 120).

- [ ] **Step 2: Full test suite**

Run: `Rscript -e 'Sys.setenv(ARES_FULL_TESTS=1); devtools::test()'`
Expected: 0 failures. All `test-plda-*` plus the pre-existing suite green.

- [ ] **Step 3: R CMD check**

Run: `Rscript -e 'devtools::check()'`
Expected: 0 errors, 0 warnings. Notes only if pre-existing.

- [ ] **Step 4: Push**

```bash
git push origin main
```

- [ ] **Step 5: Confirm CI green**

Run: `gh run list --limit 5`
Expected: `R-CMD-check`, `lint`, `pkgdown`, `test-coverage` all `success` for the head commit.
If any fail, diagnose from `gh run view <id> --log-failed`, fix, recommit, repeat.

---

## Self-review notes

- **Spec coverage:** §2 algorithm → Tasks 2–5, 10–11; §3 architecture → all tasks map to the
  listed files; §4 API → Tasks 6–7, 12; §5 parallelism → Task 13 (with sequential fallback);
  §6 testing → Tasks 4–5, 8, 11, 12, 14 + parity/determinism files; §7 scope → no tasks for
  punted items (correct); §8 invariants → Tasks 1, 13, 15, 16.
- **Type consistency:** the engine object fields (`discrim`, `mu`, `sdw`, `cmeans`, `cw`) are
  fixed in Task 4 and consumed unchanged in Tasks 6, 8, 12, 13. `pen_code` (0=L1, 1=fused) is
  consistent across Tasks 4, 6, 12. `.plda_cv` returns `lambda`/`K`/`grid`/`errors` — the stub
  in Task 6 is explicitly flagged as temporary and replaced in Task 12.
- **Known risk:** multi-class deflation and the standardization divisor are the parity-fragile
  spots; Tasks 4–5 and 11 guard them with `penalizedLDA` parity tests. If parity fails, fix
  against the reference package source — do not weaken the test tolerance.
