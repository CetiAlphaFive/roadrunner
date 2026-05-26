# test-bgam-gaussian.R
# T-GAUSS-01 through T-GAUSS-07: gaussian family correctness, selection, loss path
# Derived from test-spec.md Section 3.2

library(testthat)

test_that("T-GAUSS-01: training loss is non-increasing (gaussian)", {
  set.seed(3)
  n <- 200
  x1 <- runif(n)
  y <- sin(2 * pi * x1) + rnorm(n, 0, 0.5)
  fit <- bgam(cbind(x1), y, family = "gaussian",
              autotune = FALSE, mstop = 100, nu = 0.1, nthreads = 1)
  expect_true(all(diff(fit$loss_path) <= 1e-10),
    info = paste("max increase:", max(diff(fit$loss_path))))
})

test_that("T-GAUSS-02: true signal predictors dominate selection_frequency", {
  skip_if_quick()
  set.seed(4)
  n <- 300
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n, 0, 0.5)
  fit <- bgam(X, y, family = "gaussian",
              autotune = FALSE, mstop = 200, nu = 0.1, nthreads = 1)
  expect_gt(fit$selection_frequency[1] + fit$selection_frequency[2], 0.7)
  expect_gt(fit$selection_frequency[1], 0.3)
  expect_gt(fit$selection_frequency[2], 0.15)
})

test_that("T-GAUSS-03: predict(NULL) matches training fitted values", {
  set.seed(4)
  n <- 300
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n, 0, 0.5)
  fit <- bgam(X, y, family = "gaussian",
              autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  expect_lt(max(abs(fit$fitted.values - predict(fit))), 1e-10)
})

test_that("T-GAUSS-04: fitted values equal linear predictors for gaussian", {
  set.seed(4)
  n <- 300
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n, 0, 0.5)
  fit <- bgam(X, y, family = "gaussian",
              autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  expect_lt(max(abs(fit$fitted.values - fit$linear.predictors)), 1e-10)
})

test_that("T-GAUSS-05: predict on new data for gaussian", {
  set.seed(4)
  n <- 300
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n, 0, 0.5)
  fit <- bgam(X, y, family = "gaussian",
              autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  set.seed(99)
  X_new <- matrix(rnorm(50 * p), 50, p)
  preds <- predict(fit, newdata = X_new)
  expect_equal(length(preds), 50L)
  expect_true(all(is.finite(preds)))
  preds_link <- predict(fit, X_new, type = "link")
  preds_resp <- predict(fit, X_new, type = "response")
  expect_lt(max(abs(preds_link - preds_resp)), 1e-10)
})

# T-GAUSS-06: weighted fit changes fitted values.
# test-spec: max(abs(fit_unwt$fitted.values - fit_wt$fitted.values)) > 0.01
# BLOCKED: bgam_boost_cpp has no weights parameter; weights are validated
# and stored but never passed to the C++ engine. Weighted and unweighted
# fits are identical.
test_that("T-GAUSS-06: weighted fit differs from unweighted", {
  set.seed(5)
  n <- 100
  x <- rnorm(n)
  y <- x + rnorm(n, 0, 0.5)
  w <- runif(n, 0.5, 2)
  fit_unwt <- bgam(matrix(x), y, autotune = FALSE, mstop = 30, nu = 0.1)
  fit_wt   <- bgam(matrix(x), y, autotune = FALSE, mstop = 30, nu = 0.1, weights = w)
  expect_gt(max(abs(fit_unwt$fitted.values - fit_wt$fitted.values)), 0.01)
})

# T-GAUSS-07: type="terms" rows sum to link prediction.
# test-spec requires is.matrix(terms_mat) = TRUE
# BLOCKED: predict(fit, type="terms") with newdata=NULL returns fitted.values
# (a vector), not a terms matrix. The type="terms" path only runs when
# newdata is explicitly supplied.
test_that("T-GAUSS-07: type='terms' with newdata=X returns a matrix", {
  set.seed(4)
  n <- 300
  p <- 10
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n, 0, 0.5)
  fit <- bgam(X, y, family = "gaussian",
              autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  # When newdata is NULL, type="terms" is ignored and fitted.values are returned
  terms_mat_null <- predict(fit, type = "terms")
  # When newdata=X is explicit, the matrix is returned correctly
  terms_mat_X <- predict(fit, newdata = X, type = "terms")
  link_preds <- predict(fit, newdata = X, type = "link")
  expect_true(is.matrix(terms_mat_X))
  # With explicit newdata=X, check row sum criterion
  tol <- 1e-8
  row_sums <- rowSums(terms_mat_X)
  diff_no_intercept  <- max(abs(row_sums - (link_preds - fit$intercept_value)))
  diff_with_intercept <- max(abs(row_sums - link_preds))
  expect_true(diff_no_intercept < tol || diff_with_intercept < tol,
    info = paste("max diff (no intercept):", diff_no_intercept,
                 "max diff (with intercept):", diff_with_intercept))
  # BLOCK: test-spec calls predict(fit, type="terms") with NULL newdata
  expect_true(is.matrix(terms_mat_null))
})
