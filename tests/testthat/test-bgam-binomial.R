# test-bgam-binomial.R
# T-BINOM-01 through T-BINOM-06: binomial family correctness
# Derived from test-spec.md Section 3.3

library(testthat)

.make_binomial_fit <- function() {
  set.seed(10)
  n <- 200
  x <- rnorm(n)
  y <- rbinom(n, 1, plogis(1.5 * x))
  fit <- bgam(matrix(x), y, family = "binomial",
              autotune = FALSE, mstop = 50, nu = 0.1, nthreads = 1)
  list(fit = fit, x = x, y = y, n = n)
}

test_that("T-BINOM-01: fitted values in (0,1) for binomial", {
  res <- .make_binomial_fit()
  fit <- res$fit
  expect_true(all(fit$fitted.values > 0))
  expect_true(all(fit$fitted.values < 1))
})

test_that("T-BINOM-02: fitted values = plogis(linear.predictors)", {
  res <- .make_binomial_fit()
  fit <- res$fit
  expect_lt(max(abs(fit$fitted.values - plogis(fit$linear.predictors))), 1e-10)
})

test_that("T-BINOM-03: predict type='link' vs type='response' differ by plogis", {
  res <- .make_binomial_fit()
  fit <- res$fit
  set.seed(111)
  X_new <- matrix(rnorm(20), 20, 1)
  link_preds <- predict(fit, X_new, type = "link")
  resp_preds <- predict(fit, X_new, type = "response")
  expect_lt(max(abs(plogis(link_preds) - resp_preds)), 1e-10)
})

test_that("T-BINOM-04: loss path for binomial is non-increasing", {
  res <- .make_binomial_fit()
  fit <- res$fit
  expect_true(all(diff(fit$loss_path) <= 1e-8),
    info = paste("max increase:", max(diff(fit$loss_path))))
})

test_that("T-BINOM-05: 2-level factor y accepted and matches numeric y", {
  set.seed(10)
  n <- 200
  x <- rnorm(n)
  y_num <- rbinom(n, 1, plogis(1.5 * x))
  fit1 <- bgam(matrix(x), y_num, family = "binomial",
               autotune = FALSE, mstop = 50, nu = 0.1)
  y_factor <- factor(ifelse(y_num == 1, "yes", "no"))
  fit2 <- bgam(matrix(x), y_factor, family = "binomial",
               autotune = FALSE, mstop = 50, nu = 0.1)
  expect_no_error(fit2)
  expect_lt(max(abs(fit1$fitted.values - fit2$fitted.values)), 1e-10)
})

test_that("T-BINOM-06: logical y accepted and matches numeric y", {
  set.seed(10)
  n <- 200
  x <- rnorm(n)
  y_num <- rbinom(n, 1, plogis(1.5 * x))
  fit1 <- bgam(matrix(x), y_num, family = "binomial",
               autotune = FALSE, mstop = 50, nu = 0.1)
  y_logical <- as.logical(y_num)
  fit3 <- bgam(matrix(x), y_logical, family = "binomial",
               autotune = FALSE, mstop = 50, nu = 0.1)
  expect_no_error(fit3)
  expect_lt(max(abs(fit1$fitted.values - fit3$fitted.values)), 1e-10)
})
