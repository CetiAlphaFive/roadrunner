# test-bgam-plot.R
# T-PLOT-01 through T-PLOT-02: plot runs headless without error
# Plus property-based invariants INV-01 through INV-07
# Derived from test-spec.md Sections 3.10 and 5

library(testthat)

test_that("T-PLOT-01: plot.bgam runs headless without error", {
  set.seed(4)
  n <- 300
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n, 0, 0.5)
  fit <- bgam(X, y, family = "gaussian",
              autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  grDevices::pdf(NULL)
  expect_no_error(plot(fit))
  grDevices::dev.off()
})

test_that("T-PLOT-02: summary.bgam prints without error and produces output", {
  set.seed(4)
  n <- 300
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n, 0, 0.5)
  fit <- bgam(X, y, family = "gaussian",
              autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  out <- capture.output(summary(fit))
  expect_no_error(out)
  expect_gt(length(out), 0L)
})

# ---- Property-based invariants (INV-01 through INV-07) ----------------------

test_that("INV-01: selection frequencies sum to 1 (gaussian)", {
  set.seed(4)
  n <- 300; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n, 0, 0.5)
  fit <- bgam(X, y, family = "gaussian",
              autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  expect_lt(abs(sum(fit$selection_frequency) - 1.0), 1e-10)
})

# INV-02: selection_path covers only valid predictor indices (1-based).
# test-spec: all(fit$selection_path >= 1)
# BLOCKED: builder stores 0-based indices (from C++ 0-indexed output).
# Values are 0..(p-1) not 1..p.
test_that("INV-02: selection_path in valid range 1..p [BLOCKED - 0-based indices]", {
  set.seed(4)
  n <- 300; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n, 0, 0.5)
  fit <- bgam(X, y, family = "gaussian",
              autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  # test-spec: all(fit$selection_path >= 1) AND all(fit$selection_path <= fit$p)
  # Actual: 0-based, so min is 0 not 1
  expect_true(all(fit$selection_path >= 1L))
  expect_true(all(fit$selection_path <= fit$p))
})

# INV-03: selection frequency consistent with selection path.
# test-spec: tabulate(fit$selection_path, nbins=fit$p) / fit$mstop
# BLOCKED: 0-based selection_path makes tabulate(path, nbins=p) wrong;
# must use tabulate(path + 1L, nbins=p) to match selection_frequency.
test_that("INV-03: selection_frequency consistent with selection_path [BLOCKED - 0-based]", {
  set.seed(4)
  n <- 300; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n, 0, 0.5)
  fit <- bgam(X, y, family = "gaussian",
              autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  # test-spec formula (1-based): fails with 0-based path
  computed_freq <- tabulate(fit$selection_path, nbins = fit$p) / fit$mstop
  expect_lt(max(abs(computed_freq - fit$selection_frequency)), 1e-10)
})

test_that("INV-04: mstop consistent with path lengths", {
  set.seed(4)
  n <- 300; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n, 0, 0.5)
  fit <- bgam(X, y, family = "gaussian",
              autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  expect_equal(length(fit$selection_path), fit$mstop)
  expect_equal(length(fit$loss_path), fit$mstop)
})

test_that("INV-05: binomial fitted values are probabilities", {
  set.seed(10)
  n <- 200
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(1.5 * x))
  fit <- bgam(matrix(x), y, family = "binomial",
              autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  expect_true(all(fit$fitted.values > 0))
  expect_true(all(fit$fitted.values < 1))
})

test_that("INV-06: gaussian residuals are approximately zero-mean", {
  set.seed(4)
  n <- 300; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n, 0, 0.5)
  fit <- bgam(X, y, family = "gaussian",
              autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  expect_lt(abs(mean(fit$residuals)) / sd(fit$residuals), 0.1)
})

test_that("INV-07: loss path is non-increasing", {
  set.seed(4)
  n <- 300; p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n, 0, 0.5)
  fit <- bgam(X, y, family = "gaussian",
              autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  expect_true(all(diff(fit$loss_path) <= 1e-8))
})
