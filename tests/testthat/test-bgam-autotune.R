# test-bgam-autotune.R
# T-AUTO-01 through T-AUTO-03: CV mstop selection
# ALL TESTS GATED BY skip_if_quick()
# Derived from test-spec.md Section 3.4

library(testthat)

test_that("T-AUTO-01: autotune returns mstop_opt in valid range", {
  skip_if_quick()
  set.seed(20)
  n <- 300
  X <- matrix(rnorm(n * 5), n, 5)
  y <- sin(2 * pi * X[, 1]) + 0.3 * X[, 2] + rnorm(n, 0, 0.5)
  fit <- bgam(X, y, autotune = TRUE, mstop_max = 100, nfold = 5,
              seed.cv = 42, nu = 0.1, nthreads = 1)
  expect_false(is.null(fit$mstop_opt))
  expect_gte(fit$mstop_opt, 10L)
  expect_lte(fit$mstop_opt, 100L)
  expect_false(is.null(fit$cv))
  expect_true(is.numeric(fit$cv$cv_loss))
  expect_gte(length(fit$cv$cv_loss), 2L)
})

test_that("T-AUTO-02: CV-selected mstop beats overfit mstop on held-out RMSE", {
  skip_if_quick()
  set.seed(20)
  n <- 300
  X <- matrix(rnorm(n * 5), n, 5)
  y <- sin(2 * pi * X[, 1]) + 0.3 * X[, 2] + rnorm(n, 0, 0.5)
  set.seed(21)
  n_test <- 100
  X_test <- matrix(rnorm(n_test * 5), n_test, 5)
  y_test <- sin(2 * pi * X_test[, 1]) + 0.3 * X_test[, 2] + rnorm(n_test, 0, 0.5)
  fit_cv   <- bgam(X, y, autotune = TRUE, mstop_max = 100, seed.cv = 42,
                   nu = 0.1, nthreads = 1)
  fit_full <- bgam(X, y, autotune = FALSE, mstop = 100, nu = 0.1, nthreads = 1)
  rmse_cv   <- sqrt(mean((y_test - predict(fit_cv,   X_test))^2))
  rmse_full <- sqrt(mean((y_test - predict(fit_full, X_test))^2))
  expect_lte(rmse_cv, rmse_full * 1.05)
})

test_that("T-AUTO-03: CV result is reproducible with same seed", {
  skip_if_quick()
  set.seed(20)
  n <- 300
  X <- matrix(rnorm(n * 5), n, 5)
  y <- sin(2 * pi * X[, 1]) + 0.3 * X[, 2] + rnorm(n, 0, 0.5)
  fit_a <- bgam(X, y, autotune = TRUE, mstop_max = 80, seed.cv = 7,
                nu = 0.1, nthreads = 1)
  fit_b <- bgam(X, y, autotune = TRUE, mstop_max = 80, seed.cv = 7,
                nu = 0.1, nthreads = 1)
  expect_identical(fit_a$mstop_opt, fit_b$mstop_opt)
  expect_identical(fit_a$cv$cv_loss, fit_b$cv$cv_loss)
  expect_lt(max(abs(fit_a$fitted.values - fit_b$fitted.values)), 1e-14)
})
