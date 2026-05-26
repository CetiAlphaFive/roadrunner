# test-bgam-predict.R
# T-PRED-01 through T-PRED-05: predict types, se.fit, terms
# Also T-REG-02, T-REG-03 (regression scenarios without bagging)
# Derived from test-spec.md Sections 3.7 and 6

library(testthat)

test_that("T-PRED-01: binomial type='response' and type='link' differ by plogis", {
  set.seed(10)
  n <- 200
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(1.5 * x))
  fit <- bgam(matrix(x), y, family = "binomial",
              autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  set.seed(777)
  X_new <- matrix(rnorm(20), 20, 1)
  r <- predict(fit, X_new, type = "response")
  l <- predict(fit, X_new, type = "link")
  expect_lt(max(abs(r - plogis(l))), 1e-10)
})

# T-PRED-02: type="terms" rows sum to link prediction (gaussian).
# Uses explicit newdata=X (the NULL-newdata path is blocked separately).
test_that("T-PRED-02: type='terms' with explicit newdata rows sum to link", {
  set.seed(4)
  n <- 300
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n, 0, 0.5)
  fit <- bgam(X, y, family = "gaussian",
              autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  terms_mat  <- predict(fit, X, type = "terms")
  link_pred  <- predict(fit, X, type = "link")
  expect_true(is.matrix(terms_mat))
  expect_true(ncol(terms_mat) == length(fit$predictor_names) ||
              ncol(terms_mat) == length(fit$predictor_names) + 1L)
  row_sums <- rowSums(terms_mat)
  tol <- 1e-8
  diff_no_intercept   <- max(abs(row_sums - (link_pred - fit$intercept_value)))
  diff_with_intercept <- max(abs(row_sums - link_pred))
  expect_true(diff_no_intercept < tol || diff_with_intercept < tol,
    info = paste("max diff (no intercept):", diff_no_intercept,
                 "max diff (with intercept):", diff_with_intercept))
})

test_that("T-PRED-03: se.fit returns list with non-negative se values", {
  set.seed(4)
  n <- 300
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n, 0, 0.5)
  fit <- bgam(X, y, family = "gaussian",
              autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  set.seed(222)
  X_new <- matrix(rnorm(10 * p), 10, p)
  p_se <- predict(fit, X_new, se.fit = TRUE)
  expect_true(is.list(p_se))
  expect_true("fit" %in% names(p_se))
  expect_true("se.fit" %in% names(p_se))
  expect_true(all(p_se$se.fit >= 0))
  expect_equal(length(p_se$fit), 10L)
  expect_equal(length(p_se$se.fit), 10L)
})

# T-PRED-04: se.fit grows with extrapolation distance.
# test-spec: se_far > se_near where X_far = matrix(10, 1, 1)
# BLOCKED: B-splines evaluate to 0 outside training range, so se_far = 0
# (not larger than se_near). The SE implementation does not account for
# extrapolation uncertainty — it produces 0 SE when B-spline basis is zero.
test_that("T-PRED-04: se.fit grows with extrapolation [BLOCKED - zero SE outside range]", {
  set.seed(40)
  X_train <- matrix(rnorm(100), 100, 1)
  y <- sin(X_train) + rnorm(100, 0, 0.3)
  fit <- bgam(X_train, y, autotune = FALSE, mstop = 50, nu = 0.1)
  X_near <- matrix(0, 1, 1)
  X_far  <- matrix(10, 1, 1)
  se_near <- predict(fit, X_near, se.fit = TRUE)$se.fit
  se_far  <- predict(fit, X_far,  se.fit = TRUE)$se.fit
  # test-spec: expect se_far > se_near; actual: se_far == 0 < se_near
  expect_gt(se_far, se_near)
})

test_that("T-PRED-05: predict with newdata has correct length", {
  set.seed(4)
  n <- 300
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n, 0, 0.5)
  fit <- bgam(X, y, family = "gaussian",
              autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  set.seed(333)
  X_new <- matrix(rnorm(7 * p), 7, p)
  preds <- predict(fit, newdata = X_new)
  expect_equal(length(preds), 7L)
})

test_that("T-REG-02: formula and default interfaces give identical results", {
  set.seed(1)
  n <- 50
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- sin(x1) + rnorm(n, 0, 0.3)
  df <- data.frame(y = y, x1 = x1, x2 = x2)
  X <- cbind(x1, x2)
  fit_formula <- bgam(y ~ x1 + x2, data = df, autotune = FALSE, mstop = 20, nu = 0.1)
  fit_default <- bgam(X, y, autotune = FALSE, mstop = 20, nu = 0.1)
  expect_lt(max(abs(fit_formula$fitted.values - fit_default$fitted.values)), 1e-10)
})

test_that("T-REG-03: refit with deterministic path gives identical result", {
  set.seed(80)
  X <- matrix(rnorm(100), 100, 2)
  y <- X[, 1] + rnorm(100)
  fit1 <- bgam(X, y, autotune = FALSE, mstop = 30, nu = 0.1, nthreads = 1)
  fit2 <- bgam(X, y, autotune = FALSE, mstop = 30, nu = 0.1, nthreads = 1)
  expect_identical(fit1$fitted.values, fit2$fitted.values)
  expect_identical(fit1$selection_path, fit2$selection_path)
})
