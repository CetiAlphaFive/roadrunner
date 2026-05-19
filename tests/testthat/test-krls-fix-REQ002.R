# test-krls-fix-REQ002.R
# Unit and regression tests for REQ-20260518-002: four KRLS fixes.
# Covers: Fix 1 (sigma heuristic), Fix 2 (lambda tol), Fix 3 (autotune CV),
#         Fix 4 (autotune grid), Determinism, Parity, Regression.
#
# This file is the INDEPENDENT validation authored by tester (pipeline-isolated).
# Do NOT modify source code from this file.

library(roadrunner)

# ===========================================================================
# FIX 1: Scale-aware default sigma (median pairwise squared distance)
# ===========================================================================

# F1-1: Default sigma equals sqrt(median pairwise sq dist * p) on scaled X
test_that("F1-1: default sigma equals sqrt(median pairwise sq dist * p) on standardised X", {
  set.seed(101L)
  n <- 80L; p <- 5L
  Xs <- scale(matrix(rnorm(n * p), n, p))   # pre-standardised X

  # Expected: sqrt(median of pairwise squared Euclidean distances on Xs * p)
  sigma_expected <- sqrt(median(as.numeric(dist(Xs))^2) * p)

  fit <- krls(Xs, y = Xs[, 1] + rnorm(n, sd = 0.1),
              derivative = FALSE, vcov = FALSE)

  expect_equal(fit$sigma, sigma_expected, tolerance = 1e-10)
})

# F1-2: Default sigma is NOT equal to ncol(X) on a non-trivial DGP
test_that("F1-2: default sigma is not ncol(X) after Fix 1", {
  set.seed(102L)
  X <- matrix(rnorm(60 * 4), 60, 4)   # p = 4
  y <- X[, 1] + rnorm(60, sd = 0.2)
  fit <- krls(X, y, derivative = FALSE, vcov = FALSE)

  expect_false(
    isTRUE(all.equal(fit$sigma, ncol(X), tolerance = 0.01)),
    info = "Default sigma should not equal ncol(X) after Fix 1"
  )
})

# F1-3: Explicit sigma = ncol(X) is still honoured (backward compat)
test_that("F1-3: explicit sigma = ncol(X) is honoured", {
  set.seed(103L)
  X <- matrix(rnorm(50 * 3), 50, 3)
  y <- rnorm(50)
  p <- ncol(X)
  fit <- krls(X, y, sigma = p, derivative = FALSE, vcov = FALSE)

  expect_equal(fit$sigma, p)
})

# F1-4: Default sigma is positive and finite
test_that("F1-4: default sigma is positive and finite", {
  set.seed(104L)
  X <- matrix(rnorm(40 * 2), 40, 2)
  y <- rnorm(40)
  fit <- krls(X, y, derivative = FALSE, vcov = FALSE)

  expect_true(is.finite(fit$sigma))
  expect_true(fit$sigma > 0)
})

# F1-5: Default sigma reproducibility (same data -> same sigma)
test_that("F1-5: default sigma is deterministic (same data -> same sigma)", {
  set.seed(999L)  # generate X with fixed seed outside
  X <- matrix(rnorm(60 * 3), 60, 3)
  # generate y separately but identically
  set.seed(105L)
  y <- X[, 1] + rnorm(60, sd = 0.1)
  fit1 <- krls(X, y, derivative = FALSE, vcov = FALSE)
  fit2 <- krls(X, y, derivative = FALSE, vcov = FALSE)

  expect_identical(fit1$sigma, fit2$sigma)
})

# ===========================================================================
# FIX 2: Tighter LOO lambda tolerance
# ===========================================================================

# F2-1: Lambda is positive and finite under new default tol
test_that("F2-1: lambda is positive and finite under new and old tol", {
  set.seed(201L)
  n <- 100L; p <- 5L
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1]^2 + rnorm(n, sd = 0.5)
  sigma_fixed <- median(as.numeric(dist(scale(X)))^2)

  fit_new  <- krls(X, y, sigma = sigma_fixed, derivative = FALSE, vcov = FALSE)
  fit_old  <- krls(X, y, sigma = sigma_fixed, derivative = FALSE, vcov = FALSE,
                   tol = 1e-3 * n)

  expect_true(is.finite(fit_new$lambda) && fit_new$lambda > 0)
  expect_true(is.finite(fit_old$lambda)  && fit_old$lambda  > 0)
})

# F2-2: tol=1e-6 and tol=1e-10 yield close lambdas (precision test)
test_that("F2-2: tol=1e-6 and tol=1e-10 yield lambdas within 1e-3", {
  set.seed(202L)
  X <- matrix(rnorm(80 * 4), 80, 4)
  y <- sin(X[, 1]) + rnorm(80, sd = 0.3)
  sig <- 10

  fit_tol6  <- krls(X, y, sigma = sig, tol = 1e-6,  derivative = FALSE, vcov = FALSE)
  fit_tol10 <- krls(X, y, sigma = sig, tol = 1e-10, derivative = FALSE, vcov = FALSE)

  # With tol=1e-6, lambda must be within 1e-3 of the tol=1e-10 result.
  expect_equal(fit_tol6$lambda, fit_tol10$lambda, tolerance = 1e-3)
})

# F2-3: Default tol is tighter than old tol=1e-3*n (closer to tol=1e-10 reference)
test_that("F2-3: default tol yields lambda at least as precise as old tol=1e-3*n", {
  set.seed(203L)
  n <- 200L
  X <- matrix(rnorm(n * 5), n, 5)
  y <- X[, 1] + 0.5 * X[, 2]^2 + rnorm(n, sd = 0.5)
  sig <- 12

  fit_default  <- krls(X, y, sigma = sig, derivative = FALSE, vcov = FALSE)
  fit_old_tol  <- krls(X, y, sigma = sig, tol = 1e-3 * n, derivative = FALSE, vcov = FALSE)
  fit_ref_tol  <- krls(X, y, sigma = sig, tol = 1e-10, derivative = FALSE, vcov = FALSE)

  # new default must be closer (or equal) to the tol=1e-10 reference
  err_new <- abs(fit_default$lambda - fit_ref_tol$lambda)
  err_old <- abs(fit_old_tol$lambda - fit_ref_tol$lambda)
  expect_lte(err_new, err_old + 1e-6,
             label = "new default tol yields lambda at least as precise as old tol")
})

# F2-4: Explicit tol argument is honoured
test_that("F2-4: explicit tol argument produces a valid lambda", {
  set.seed(204L)
  X <- matrix(rnorm(60 * 3), 60, 3)
  y <- rnorm(60)
  fit <- krls(X, y, tol = 0.01, derivative = FALSE, vcov = FALSE)

  expect_true(is.finite(fit$lambda) && fit$lambda > 0)
})

# ===========================================================================
# FIX 3: Stabilised autotune CV (nfold=10, ncross=2, 1-SE rule, new fields)
# ===========================================================================

# F3-1: Default autotune nfold is 10 when autotune=TRUE
test_that("F3-1: default autotune nfold is 10", {
  set.seed(301L)
  X <- matrix(rnorm(80 * 3), 80, 3)
  y <- sin(X[, 1]) + rnorm(80, sd = 0.2)
  fit <- krls(X, y, autotune = TRUE, seed.cv = 1L,
              derivative = FALSE, vcov = FALSE)

  expect_equal(fit$autotune$nfold, 10L)
})

# F3-2: Default autotune ncross is 2
test_that("F3-2: default autotune ncross is 2", {
  set.seed(301L)
  X <- matrix(rnorm(80 * 3), 80, 3)
  y <- sin(X[, 1]) + rnorm(80, sd = 0.2)
  fit <- krls(X, y, autotune = TRUE, seed.cv = 1L,
              derivative = FALSE, vcov = FALSE)

  expect_equal(fit$autotune$ncross, 2L)
})

# F3-3: Explicit nfold and ncross override the defaults
test_that("F3-3: explicit nfold and ncross override defaults", {
  set.seed(302L)
  X <- matrix(rnorm(60 * 2), 60, 2)
  y <- X[, 1] + rnorm(60, sd = 0.1)
  fit <- krls(X, y, autotune = TRUE, nfold = 5L, ncross = 1L, seed.cv = 2L,
              derivative = FALSE, vcov = FALSE)

  expect_equal(fit$autotune$nfold, 5L)
  expect_equal(fit$autotune$ncross, 1L)
})

# F3-4: 1-SE rule causes selected sigma >= argmin sigma
test_that("F3-4: 1-SE rule selects sigma >= raw argmin", {
  set.seed(303L)
  X <- matrix(rnorm(100 * 3), 100, 3)
  y <- 0.5 * X[, 1] + 0.3 * X[, 2] + rnorm(100, sd = 0.2)
  fit <- krls(X, y, autotune = TRUE, nfold = 5L, ncross = 2L, seed.cv = 3L,
              derivative = FALSE, vcov = FALSE)

  # raw argmin sigma (before 1-SE rule)
  best_raw <- fit$autotune$grid[which.min(fit$autotune$mse)]

  expect_gte(fit$autotune$winner, best_raw,
             label = "1-SE rule selects sigma >= argmin sigma")
})

# F3-5: Autotune reproducibility with seed.cv
test_that("F3-5: autotune is reproducible with seed.cv", {
  set.seed(304L)
  X <- matrix(rnorm(80 * 3), 80, 3)
  y <- sin(X[, 1]) + rnorm(80, sd = 0.2)

  f1 <- krls(X, y, autotune = TRUE, seed.cv = 9L, derivative = FALSE, vcov = FALSE)
  f2 <- krls(X, y, autotune = TRUE, seed.cv = 9L, derivative = FALSE, vcov = FALSE)

  expect_identical(f1$autotune$winner, f2$autotune$winner)
  expect_equal(f1$autotune$mse, f2$autotune$mse, tolerance = 1e-14)
})

# F3-6: autotune_info list has new fields (ncross, mse_per_fold, se_mse)
test_that("F3-6: autotune_info contains ncross, mse_per_fold, se_mse", {
  set.seed(301L)
  X <- matrix(rnorm(80 * 3), 80, 3)
  y <- sin(X[, 1]) + rnorm(80, sd = 0.2)
  fit <- krls(X, y, autotune = TRUE, seed.cv = 1L,
              derivative = FALSE, vcov = FALSE)
  at <- fit$autotune

  expect_true(!is.null(at$ncross))
  expect_true(!is.null(at$mse_per_fold))
  expect_true(!is.null(at$se_mse))

  # mse_per_fold dimensions: n_sigma rows x (nfold * ncross) cols
  n_sig <- length(at$grid)
  expect_equal(nrow(at$mse_per_fold), n_sig)
  expect_equal(ncol(at$mse_per_fold), at$nfold * at$ncross)
})

# F3-7: lambda.method='cv' path with ncross=NULL defaults to ncross=1
test_that("F3-7: lambda.method='cv' with ncross=NULL defaults ncross=1", {
  set.seed(305L)
  X <- matrix(rnorm(60 * 3), 60, 3)
  y <- rnorm(60)
  fit <- krls(X, y, lambda.method = "cv", nfold = 5L, seed.cv = 5L,
              derivative = FALSE, vcov = FALSE)

  expect_equal(fit$cv$ncross, 1L)
})

# ===========================================================================
# FIX 4: Autotune grid range and centering
# ===========================================================================

# F4-1: Default autotune grid has 9 points
test_that("F4-1: default autotune grid has 9 points", {
  set.seed(401L)
  X <- matrix(rnorm(80 * 3), 80, 3)
  y <- sin(X[, 1]) + rnorm(80, sd = 0.2)
  fit <- krls(X, y, autotune = TRUE, seed.cv = 1L,
              derivative = FALSE, vcov = FALSE)

  expect_equal(length(fit$autotune$grid), 9L)
})

# F4-2: Default grid is NOT the old d * {0.25..8}
test_that("F4-2: default grid is not the old d * {0.25, 0.5, 1, 2, 4, 8}", {
  set.seed(401L)
  X <- matrix(rnorm(80 * 3), 80, 3)
  y <- sin(X[, 1]) + rnorm(80, sd = 0.2)
  fit <- krls(X, y, autotune = TRUE, seed.cv = 1L,
              derivative = FALSE, vcov = FALSE)

  old_grid <- ncol(X) * c(0.25, 0.5, 1, 2, 4, 8)
  expect_false(
    isTRUE(all.equal(sort(fit$autotune$grid), sort(old_grid))),
    info = "new grid should differ from old d * {0.25..8}"
  )
})

# F4-3: Default grid is centred on sigma_anchor with ratios {0.125..32}
test_that("F4-3: default grid is sigma_anchor * {0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32}", {
  set.seed(402L)
  n <- 80L; p <- 3L
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(X[, 1]) + rnorm(n, sd = 0.2)

  # Compute expected sigma_anchor: sqrt(median pairwise sq dist on scale(X) * p)
  Xs_ref <- scale(X)
  sigma_anchor_ref <- sqrt(median(as.numeric(dist(Xs_ref))^2) * p)

  fit <- krls(X, y, autotune = TRUE, seed.cv = 1L,
              derivative = FALSE, vcov = FALSE)

  expected_grid <- sort(sigma_anchor_ref * c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32))

  expect_equal(sort(fit$autotune$grid), expected_grid, tolerance = 1e-10)
})

# F4-4: Custom autotune.grid overrides the new default (backward compat)
test_that("F4-4: custom autotune.grid overrides the default", {
  set.seed(403L)
  X <- matrix(rnorm(60 * 3), 60, 3)
  y <- rnorm(60)
  custom_grid <- c(1, 5, 10, 20)
  fit <- krls(X, y, autotune = TRUE, autotune.grid = custom_grid,
              seed.cv = 1L, derivative = FALSE, vcov = FALSE)

  expect_equal(fit$autotune$grid, custom_grid)
})

# ===========================================================================
# DETERMINISM TESTS
# ===========================================================================

# DET-1: nthreads=1 vs nthreads=4, default call
test_that("DET-1: nthreads=1 vs nthreads=4 byte-identical (no autotune)", {
  set.seed(999L)
  X <- matrix(rnorm(80 * 4), 80, 4)
  y <- sin(X[, 1]) + 0.3 * X[, 2]^2 + rnorm(80, sd = 0.2)

  f1 <- krls(X, y, nthreads = 1L, derivative = FALSE, vcov = FALSE)
  f4 <- krls(X, y, nthreads = 4L, derivative = FALSE, vcov = FALSE)

  expect_identical(f1$coeffs, f4$coeffs)
  expect_identical(f1$fitted, f4$fitted)
  expect_identical(f1$sigma,  f4$sigma)
  expect_identical(f1$lambda, f4$lambda)
})

# DET-2: nthreads=1 vs nthreads=4, with autotune
test_that("DET-2: nthreads=1 vs nthreads=4 byte-identical (with autotune)", {
  skip_if_quick()
  set.seed(998L)
  X <- matrix(rnorm(80 * 3), 80, 3)
  y <- X[, 1] + rnorm(80, sd = 0.2)

  f1 <- krls(X, y, autotune = TRUE, nthreads = 1L, seed.cv = 42L,
             derivative = FALSE, vcov = FALSE)
  f4 <- krls(X, y, autotune = TRUE, nthreads = 4L, seed.cv = 42L,
             derivative = FALSE, vcov = FALSE)

  expect_identical(f1$sigma,  f4$sigma)
  expect_identical(f1$lambda, f4$lambda)
  expect_identical(f1$coeffs, f4$coeffs)
})

# ===========================================================================
# PARITY TESTS (KRLS package, if installed)
# ===========================================================================

# PAR-1 + PAR-2: Run existing parity test suite via sourcing
# The existing test-krls.R already runs these; we verify here that at fixed
# (sigma, lambda) the coefficients still match KRLS::krls.
test_that("PAR-1: at fixed sigma+lambda, coefficients match KRLS::krls", {
  testthat::skip_if_not_installed("KRLS")
  set.seed(1L)
  n <- 60L; p <- 3L
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  y <- sin(X[, 1L]) + 0.5 * X[, 2L]^2 - 0.3 * X[, 3L] + rnorm(n, sd = 0.1)

  sigma <- ncol(X); lambda_val <- 0.5
  ref <- KRLS::krls(X, y, sigma = sigma, lambda = lambda_val,
                    derivative = FALSE, vcov = TRUE, print.level = 0)
  fit <- roadrunner::krls(X, y, sigma = sigma, lambda = lambda_val,
                          derivative = FALSE, vcov = TRUE, print.level = 0)

  expect_equal(as.numeric(fit$coeffs), as.numeric(ref$coeffs), tolerance = 1e-9)
  expect_equal(as.numeric(fit$fitted), as.numeric(ref$fitted), tolerance = 1e-9)
})

# PAR-2: Auto-lambda with new tol=1e-6 still agrees with KRLS within old tol
test_that("PAR-2: auto-lambda agrees with KRLS::krls within the old tolerance (backward compat)", {
  testthat::skip_if_not_installed("KRLS")
  set.seed(2L)
  n <- 60L; p <- 3L
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  y <- sin(X[, 1L]) + 0.5 * X[, 2L]^2 - 0.3 * X[, 3L] + rnorm(n, sd = 0.1)
  sigma <- ncol(X)

  ref <- KRLS::krls(X, y, sigma = sigma, derivative = FALSE, vcov = TRUE,
                    print.level = 0)
  fit <- roadrunner::krls(X, y, sigma = sigma, derivative = FALSE, vcov = TRUE,
                          print.level = 0)

  # New lambda is more precise; old tolerance 1e-3*n = 0.06 still holds at lambda level.
  # Coefficient differences scale with lambda differences; allow moderate tolerance.
  expect_equal(fit$lambda, ref$lambda, tolerance = 1e-3 * n)
})

# ===========================================================================
# REGRESSION: existing test suite must stay green
# ===========================================================================

# REG-1 is validated by running devtools::test(filter="krls") — see run step.
# This file does not re-implement those tests; they are checked separately.
