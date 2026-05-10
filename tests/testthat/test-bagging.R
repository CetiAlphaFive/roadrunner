# Tests for n.boot bagging (Phase 2 — v0.18+).

test_that("n.boot > 0 produces a $boot slot with the right number of fits", {
  set.seed(20260510)
  x <- matrix(stats::runif(150 * 4), 150, 4)
  y <- 5 * x[, 1] + stats::rnorm(150)
  fit <- ares(x, y, n.boot = 10L, seed.cv = 1, nthreads = 2)
  expect_true("boot" %in% names(fit))
  expect_equal(fit$boot$n.boot, 10L)
  expect_length(fit$boot$fits, 10L)
})

test_that("n.boot = 0 (default) does not create a $boot slot", {
  set.seed(20260510)
  x <- matrix(stats::runif(100 * 4), 100, 4)
  y <- 5 * x[, 1] + stats::rnorm(100)
  fit <- ares(x, y, nthreads = 2)
  expect_false("boot" %in% names(fit))
})

test_that("predict(.., se.fit = TRUE) returns sd attribute matching length", {
  set.seed(20260510)
  x <- matrix(stats::runif(150 * 4), 150, 4)
  y <- 5 * x[, 1] + 3 * x[, 2] + stats::rnorm(150)
  fit <- ares(x, y, n.boot = 8L, seed.cv = 1, nthreads = 2)
  xte <- matrix(stats::runif(40 * 4), 40, 4)
  yhat <- predict(fit, xte, se.fit = TRUE)
  expect_length(yhat, 40L)
  sds <- attr(yhat, "sd")
  expect_length(sds, 40L)
  expect_true(all(is.finite(sds)))
  expect_true(all(sds >= 0))
})

test_that("predict without se.fit returns plain numeric for bag fits", {
  set.seed(20260510)
  x <- matrix(stats::runif(120 * 4), 120, 4)
  y <- 5 * x[, 1] + stats::rnorm(120)
  fit <- ares(x, y, n.boot = 5L, seed.cv = 1, nthreads = 2)
  xte <- matrix(stats::runif(20 * 4), 20, 4)
  yhat <- predict(fit, xte)
  expect_null(attr(yhat, "sd"))
  expect_length(yhat, 20L)
})

test_that("bagging is deterministic at fixed seed.cv across nthreads", {
  set.seed(20260510)
  x <- matrix(stats::runif(120 * 4), 120, 4)
  y <- 5 * x[, 1] + stats::rnorm(120)
  xte <- matrix(stats::runif(30 * 4), 30, 4)
  f1 <- ares(x, y, n.boot = 5L, seed.cv = 42, nthreads = 1)
  f4 <- ares(x, y, n.boot = 5L, seed.cv = 42, nthreads = 4)
  p1 <- predict(f1, xte, se.fit = TRUE)
  p4 <- predict(f4, xte, se.fit = TRUE)
  expect_equal(as.numeric(p1), as.numeric(p4), tolerance = 1e-12)
  expect_equal(attr(p1, "sd"), attr(p4, "sd"), tolerance = 1e-12)
})

test_that("autotune composes with n.boot", {
  set.seed(20260510)
  x <- matrix(stats::runif(180 * 5), 180, 5)
  y <- 10 * sin(pi * x[, 1] * x[, 2]) + stats::rnorm(180)
  fit <- ares(x, y, autotune = TRUE, n.boot = 5L,
              seed.cv = 1, nthreads = 2)
  expect_true("autotune" %in% names(fit))
  expect_true("boot" %in% names(fit))
  # Each replicate should have used the central fit's tuned hyperparams,
  # so its forward-pass nk equals the autotune winner's nk.
  for (b in fit$boot$fits) {
    expect_equal(b$nk, fit$autotune$nk)
  }
})

test_that("seed.cv preserves user RNG with bagging", {
  x <- matrix(stats::runif(100 * 3), 100, 3)
  y <- 5 * x[, 1] + stats::rnorm(100)
  set.seed(123)
  before <- stats::runif(1)
  set.seed(123)
  fit <- ares(x, y, n.boot = 5L, seed.cv = 99, nthreads = 2)
  after <- stats::runif(1)
  expect_equal(before, after, tolerance = 1e-12)
})
