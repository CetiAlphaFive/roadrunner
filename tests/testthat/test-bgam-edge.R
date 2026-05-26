# test-bgam-edge.R
# T-EDGE-01 through T-EDGE-08: edge case scenarios
# Derived from test-spec.md Section 4

library(testthat)

test_that("T-EDGE-01: single predictor works and is always selected", {
  set.seed(70)
  n <- 50
  X <- matrix(rnorm(n), n, 1)
  y <- X[, 1] + rnorm(n)
  fit <- bgam(X, y, autotune = FALSE, mstop = 10, nu = 0.1)
  expect_no_error(fit)
  expect_equal(length(fit$fitted.values), 50L)
  expect_equal(length(fit$selection_frequency), 1L)
  # Use expect_true to avoid names attribute mismatch with expect_equal
  expect_true(fit$selection_frequency[1] == 1.0)
})

# T-EDGE-02: zero-column predictor matrix rejected with clear error.
# test-spec: stops with message containing "at least 1 column"
# BLOCKED: builder does not guard against 0-column x before setting colnames;
# the internal R error is "length of 'dimnames' [2] not equal to array extent"
# rather than the expected clear user-facing message.
test_that("T-EDGE-02: zero-column predictor matrix rejected", {
  X <- matrix(numeric(0), nrow = 20, ncol = 0)
  y <- rnorm(20)
  err_msg <- tryCatch(
    { bgam(X, y, autotune = FALSE, mstop = 5, nu = 0.1); NA_character_ },
    error = function(e) conditionMessage(e)
  )
  # test-spec: message should contain "at least 1 column"
  # Actual: "length of 'dimnames' [2] not equal to array extent"
  expect_match(err_msg, "at least 1 column")
})

test_that("T-EDGE-03: n=10 (minimum) works", {
  set.seed(71)
  n <- 10
  X <- matrix(rnorm(10), 10, 1)
  y <- rnorm(10)
  fit <- bgam(X, y, autotune = FALSE, mstop = 5, nu = 0.1, nknots = 3)
  expect_no_error(fit)
  expect_equal(length(fit$fitted.values), 10L)
})

test_that("T-EDGE-04: n=9 (below minimum) stopped with clear error", {
  set.seed(72)
  n <- 9
  X <- matrix(rnorm(9), 9, 1)
  y <- rnorm(9)
  expect_error(bgam(X, y, autotune = FALSE, mstop = 5, nu = 0.1),
               "at least 10 observations")
})

test_that("T-EDGE-05: near-constant predictor switches to linear with warning", {
  set.seed(73)
  n <- 50
  X <- cbind(rnorm(n), rep(3.0, n) + rnorm(n) * 1e-12)
  y <- X[, 1] + rnorm(n)
  expect_warning(
    fit <- bgam(X, y, autotune = FALSE, mstop = 20, nu = 0.1),
    regexp = "near-constant|near constant"
  )
  expect_no_error(fit)
  expect_gt(fit$selection_frequency[1], fit$selection_frequency[2])
})

test_that("T-EDGE-06: all-equal y (zero-variance response) does not crash", {
  X <- matrix(rnorm(50), 50, 2)
  y <- rep(1, 50)
  fit <- bgam(X, y, family = "gaussian", autotune = FALSE, mstop = 5, nu = 0.1)
  expect_no_error(fit)
  expect_true(all(abs(fit$fitted.values - 1) < 0.01))
  expect_true(all(fit$loss_path < 1e-6))
})

test_that("T-EDGE-07: factor predictor in formula via model.matrix expansion", {
  set.seed(74)
  n <- 60
  df <- data.frame(
    y = rnorm(n),
    g = factor(rep(c("A", "B", "C"), 20)),
    x = rnorm(n)
  )
  fit <- bgam(y ~ g + x, data = df, autotune = FALSE, mstop = 20, nu = 0.1)
  expect_no_error(fit)
  expect_s3_class(fit, "bgam")
  expect_equal(length(fit$fitted.values), 60L)
})

test_that("T-EDGE-08: p > n (large p) does not crash", {
  set.seed(75)
  n <- 50
  p <- 100
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] + rnorm(n)
  fit <- bgam(X, y, autotune = FALSE, mstop = 20, nu = 0.1)
  expect_no_error(fit)
  expect_equal(length(fit$fitted.values), n)
})
