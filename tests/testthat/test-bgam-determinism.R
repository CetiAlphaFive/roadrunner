# test-bgam-determinism.R
# T-DET-01 through T-DET-03: thread-identical results
# Derived from test-spec.md Section 3.5

library(testthat)

test_that("T-DET-01: nthreads=1 identical to nthreads=4 (gaussian)", {
  set.seed(99)
  n <- 100
  X <- matrix(rnorm(n * 5), n, 5)
  y <- X[, 1] + rnorm(n)
  fit1 <- bgam(X, y, autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  fit4 <- bgam(X, y, autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 4)
  expect_identical(fit1$fitted.values, fit4$fitted.values)
  expect_identical(fit1$selection_path, fit4$selection_path)
  expect_identical(fit1$loss_path, fit4$loss_path)
})

test_that("T-DET-02: nthreads=1 identical to nthreads=4 (binomial)", {
  set.seed(100)
  n <- 100
  X <- matrix(rnorm(n * 4), n, 4)
  y <- rbinom(n, 1, plogis(X[, 1]))
  fit1 <- bgam(X, y, family = "binomial", autotune = FALSE, mstop = 50,
               nu = 0.1, nthreads = 1)
  fit4 <- bgam(X, y, family = "binomial", autotune = FALSE, mstop = 50,
               nu = 0.1, nthreads = 4)
  expect_identical(fit1$fitted.values, fit4$fitted.values)
  expect_identical(fit1$selection_path, fit4$selection_path)
})

test_that("T-DET-03: autotune CV is thread-identical", {
  skip_if_quick()
  set.seed(20)
  n <- 300
  X_full <- matrix(rnorm(n * 5), n, 5)
  y_full <- sin(2 * pi * X_full[, 1]) + 0.3 * X_full[, 2] + rnorm(n, 0, 0.5)
  X <- X_full[1:100, ]
  y <- y_full[1:100]
  fit_t1 <- bgam(X, y, autotune = TRUE, mstop_max = 50, seed.cv = 5,
                 nu = 0.1, nthreads = 1)
  fit_t4 <- bgam(X, y, autotune = TRUE, mstop_max = 50, seed.cv = 5,
                 nu = 0.1, nthreads = 4)
  expect_identical(fit_t1$mstop_opt, fit_t4$mstop_opt)
  expect_identical(fit_t1$fitted.values, fit_t4$fitted.values)
})
