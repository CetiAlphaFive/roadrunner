# test-bgam-api.R
# T-API-01 through T-API-12: API basics, error handling, formula/default interfaces
# Derived from test-spec.md Section 3.1

library(testthat)

test_that("T-API-01: formula interface works", {
  set.seed(1)
  n <- 50
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- sin(x1) + rnorm(n, 0, 0.3)
  df <- data.frame(y = y, x1 = x1, x2 = x2)
  fit <- bgam(y ~ x1 + x2, data = df, autotune = FALSE, mstop = 20, nu = 0.1)
  expect_s3_class(fit, "bgam")
  expect_equal(length(fit$fitted.values), 50L)
})

test_that("T-API-02: default (matrix) interface works and agrees with formula", {
  set.seed(1)
  n <- 50
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- sin(x1) + rnorm(n, 0, 0.3)
  df <- data.frame(y = y, x1 = x1, x2 = x2)
  X <- cbind(x1, x2)
  fit_formula <- bgam(y ~ x1 + x2, data = df, autotune = FALSE, mstop = 20, nu = 0.1)
  fit_default <- bgam(x = X, y = y, autotune = FALSE, mstop = 20, nu = 0.1)
  expect_s3_class(fit_default, "bgam")
  expect_equal(length(fit_default$fitted.values), 50L)
  expect_lt(max(abs(fit_formula$fitted.values - fit_default$fitted.values)), 1e-10)
})

# T-API-03: unsupported family rejected with message "family must be".
# test-spec: stops with message containing "family must be"
# BLOCKED: builder uses match.arg() which produces:
#   "'arg' should be one of 'gaussian', 'binomial'"
# This does not contain "family must be". Builder must add an explicit
# stop() with the required message before match.arg().
test_that("T-API-03: unsupported family rejected", {
  n <- 20
  x <- matrix(rnorm(n * 2), n, 2)
  y <- rnorm(n)
  err_msg <- tryCatch(
    { bgam(x, y, family = "poisson"); NA_character_ },
    error = function(e) conditionMessage(e)
  )
  # test-spec requires this to contain "family must be"
  expect_match(err_msg, "family must be")
})

test_that("T-API-04: zero weights rejected", {
  set.seed(1)
  n <- 20
  x <- matrix(rnorm(n * 2), n, 2)
  y <- rnorm(n)
  w <- c(0, rep(1, n - 1))
  expect_error(bgam(x, y, weights = w), "weights must be strictly positive")
})

test_that("T-API-05: negative weights rejected", {
  set.seed(1)
  n <- 20
  x <- matrix(rnorm(n * 2), n, 2)
  y <- rnorm(n)
  w <- c(-1, rep(1, n - 1))
  expect_error(bgam(x, y, weights = w), "weights must be strictly positive")
})

test_that("T-API-06: NA weights rejected", {
  set.seed(1)
  n <- 20
  x <- matrix(rnorm(n * 2), n, 2)
  y <- rnorm(n)
  w <- c(NA, rep(1, n - 1))
  expect_error(bgam(x, y, weights = w), "weights must be strictly positive")
})

test_that("T-API-07: mstop < 1 rejected", {
  set.seed(1)
  n <- 20
  x <- matrix(rnorm(n * 2), n, 2)
  y <- rnorm(n)
  expect_error(bgam(x, y, mstop = 0), "mstop must be >= 1")
})

test_that("T-API-08: nu out of range rejected", {
  set.seed(1)
  n <- 20
  x <- matrix(rnorm(n * 2), n, 2)
  y <- rnorm(n)
  expect_error(bgam(x, y, nu = 0),   "nu must be in")
  expect_error(bgam(x, y, nu = 1.5), "nu must be in")
})

test_that("T-API-09: binomial y not binary rejected", {
  set.seed(1)
  n <- 50
  x <- matrix(rnorm(n * 2), n, 2)
  y <- rnorm(n)
  expect_error(bgam(x, y, family = "binomial"), "binary")
})

test_that("T-API-10: OOV factor level in predict errors clearly", {
  set.seed(2)
  n <- 20
  df <- data.frame(y = rnorm(n), g = factor(rep(c("A", "B"), 10)), x = rnorm(n))
  fit <- bgam(y ~ g + x, data = df, autotune = FALSE, mstop = 10)
  newdf <- data.frame(g = factor("C"), x = 0)
  expect_error(predict(fit, newdata = newdf), "not seen in training|factor level")
})

test_that("T-API-11: fit object has all required fields", {
  set.seed(1)
  n <- 50
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- sin(x1) + rnorm(n, 0, 0.3)
  df <- data.frame(y = y, x1 = x1, x2 = x2)
  fit <- bgam(y ~ x1 + x2, data = df, autotune = FALSE, mstop = 20, nu = 0.1)
  required <- c("coefficients", "fitted.values", "linear.predictors",
                "residuals", "selection_path", "selection_frequency",
                "loss_path", "nu", "mstop", "family", "base_learners",
                "predictor_names", "call")
  for (nm in required) {
    expect_true(!is.null(fit[[nm]]), info = paste("field missing:", nm))
  }
  expect_equal(length(fit$selection_path), fit$mstop)
  expect_equal(length(fit$loss_path), fit$mstop)
  expect_equal(length(fit$selection_frequency), fit$p)
  expect_lt(abs(sum(fit$selection_frequency) - 1.0), 1e-10)
  expect_true(all(fit$selection_frequency >= 0))
})

test_that("T-API-12: print and summary run without error", {
  set.seed(1)
  n <- 50
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- sin(x1) + rnorm(n, 0, 0.3)
  df <- data.frame(y = y, x1 = x1, x2 = x2)
  fit <- bgam(y ~ x1 + x2, data = df, autotune = FALSE, mstop = 20, nu = 0.1)
  expect_invisible(print(fit))
  expect_identical(print(fit), fit)
  expect_no_error(summary(fit))
})
