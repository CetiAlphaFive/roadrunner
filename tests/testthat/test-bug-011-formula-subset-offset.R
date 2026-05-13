# BUG-011 regression test (v0.0.0.9032).
#
# ares.formula() previously silently dropped `subset = ...` (it fell into
# `...` and went nowhere) and silently absorbed `offset(...)` terms as
# ordinary predictors. Both produced wrong fits with zero indication.
#
# Fix: add explicit `subset` arg to ares.formula and pass to model.frame;
# detect offset() via attr(terms, "offset") before model.frame runs and
# error with an actionable message.

test_that("ares.formula honours subset = ...", {
  set.seed(1)
  df <- data.frame(x1 = rnorm(100), x2 = rnorm(100))
  df$y <- df$x1 - 0.5 * df$x2 + rnorm(100, sd = 0.3)
  fit_full <- ares(y ~ x1 + x2, data = df, nthreads = 2)
  fit_sub  <- ares(y ~ x1 + x2, data = df, subset = 1:50, nthreads = 2)
  # The subset fit was trained on 50 rows; fitted.values should reflect
  # that. The full fit's fitted.values length is 100.
  expect_equal(length(fit_full$fitted.values), 100L)
  expect_equal(length(fit_sub$fitted.values), 50L)
})

test_that("ares.formula errors on offset() terms with a clear message", {
  set.seed(2)
  df <- data.frame(x1 = rnorm(60), x2 = rnorm(60),
                   log_exp = log(runif(60, 0.5, 2)))
  df$y <- rpois(60, exp(0.3 * df$x1 + df$log_exp))
  expect_error(
    ares(y ~ x1 + x2 + offset(log_exp),
         data = df, family = "poisson", nthreads = 2),
    "offset"
  )
})

test_that("plain formula without subset/offset works (no regression)", {
  set.seed(3)
  df <- data.frame(x1 = rnorm(80), x2 = rnorm(80))
  df$y <- df$x1 - 0.5 * df$x2 + rnorm(80, sd = 0.3)
  fit <- ares(y ~ x1 + x2, data = df, nthreads = 2)
  expect_s3_class(fit, "ares")
  expect_equal(length(fit$fitted.values), 80L)
})
