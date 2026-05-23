# BUG-019 regression test (v0.0.0.9055).
#
# meep() silently fit a regression on integer level codes when y was a
# factor (or character) with 3+ levels. The `y <- as.numeric(y)` coerce
# happened before family detection; the resulting 1, 2, 3 codes have
# distinct-count > 2, so the binary check failed and family defaulted to
# "gaussian". This violates meep()'s scope (binary or continuous only)
# and produced silently wrong fits.
#
# Fix: reject factor / character outcomes with > 2 levels upfront with a
# clear message, mirroring the treatment-validation pattern.

test_that("meep() errors on factor outcome with 3+ levels", {
  set.seed(1)
  n <- 60
  X <- matrix(runif(n * 3), n, 3)
  yf <- factor(sample(c("a", "b", "c"), n, replace = TRUE))
  expect_error(meep(X, yf, folds = 4L, tune = "none", seed = 1),
               "Multi-class outcomes are not supported")
})

test_that("meep() errors on character outcome with 3+ levels", {
  set.seed(2)
  n <- 60
  X <- matrix(runif(n * 3), n, 3)
  yc <- sample(c("a", "b", "c"), n, replace = TRUE)
  expect_error(meep(X, yc, folds = 4L, tune = "none", seed = 1),
               "Multi-class outcomes are not supported")
})

test_that("meep() still accepts binary factor outcomes", {
  set.seed(3)
  n <- 60
  X <- matrix(runif(n * 3), n, 3)
  yb <- factor(sample(c("yes", "no"), n, replace = TRUE))
  # `logreg` is the family-compatible linear learner for a binary
  # outcome; `ols` would (correctly) error as family-incompatible.
  expect_no_error(
    meep(X, yb, folds = 4L, learners = "logreg", tune = "none", seed = 1)
  )
})
