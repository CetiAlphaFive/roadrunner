# Fit-quality tests on Friedman-style DGPs. ares is no longer required to
# match earth; instead we check that fits explain a reasonable share of the
# DGP's signal and that ares does not under-perform earth by a wide margin.
# (ares may legitimately beat earth on rank-deficient designs because it
# returns the true OLS minimum where earth's rank-clamped path does not.)

test_that("ares fits Friedman-1 deg=2 with R^2 >= 0.92 on signal", {
  skip_on_cran()
  set.seed(20260509)
  n <- 300; p <- 10
  x <- matrix(stats::runif(n * p), n, p)
  signal <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
            10 * x[, 4] + 5 * x[, 5]
  y <- signal + stats::rnorm(n)
  fa <- ares(x, y, degree = 2, nthreads = 2)
  pa <- predict(fa, x)
  # R^2 vs DGP signal (population-truth R^2 — a fit-quality test that does
  # not depend on earth at all).
  ss_signal <- sum((signal - mean(signal))^2)
  ss_resid_signal <- sum((signal - pa)^2)
  r2_signal <- 1 - ss_resid_signal / ss_signal
  expect_gt(r2_signal, 0.92)
})

test_that("ares fits Friedman-1 deg=1 with R^2 >= 0.96 on signal", {
  skip_on_cran()
  set.seed(20260509)
  n <- 200; p <- 5
  x <- matrix(stats::runif(n * p), n, p)
  signal <- 10 * x[, 1] + 5 * x[, 2] + 3 * x[, 3]
  y <- signal + stats::rnorm(n)
  fa <- ares(x, y, degree = 1, nthreads = 2)
  pa <- predict(fa, x)
  ss_signal <- sum((signal - mean(signal))^2)
  ss_resid_signal <- sum((signal - pa)^2)
  r2_signal <- 1 - ss_resid_signal / ss_signal
  expect_gt(r2_signal, 0.96)
})

test_that("ares does not lose to earth by a wide margin (informational)", {
  skip_if_not_installed("earth")
  skip_on_cran()
  set.seed(20260509)
  n <- 300; p <- 10
  x <- matrix(stats::runif(n * p), n, p)
  y <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
       10 * x[, 4] + 5 * x[, 5] + stats::rnorm(n)
  fa <- ares(x, y, degree = 2, nthreads = 2)
  fe <- earth::earth(x, y, degree = 2)
  ratio <- fa$rss / fe$rss
  # ares may legitimately match or beat earth; only fail if ares is more
  # than 25% worse (which would indicate a forward-pass regression).
  expect_lt(ratio, 1.25)
})
