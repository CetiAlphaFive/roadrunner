# BUG-007 regression test (v0.0.0.9029).
#
# family = "poisson" used to accept y = rep(0, n) silently, producing a
# degenerate intercept-only fit with deviance = 0 and coefficients tending
# to -Inf under the log link. More broadly, a constant y (all observed
# values equal) is degenerate for ANY family: only the intercept can be
# fit. v0.0.0.9029 rejects both up front so the user gets a clear signal
# instead of a silently degenerate fit.

test_that("family='poisson' rejects all-zero y with a clear error", {
  set.seed(31)
  n <- 100; p <- 3
  x <- matrix(stats::runif(n * p), n, p)
  y0 <- rep(0L, n)
  expect_error(ares(x, y0, family = "poisson", nthreads = 2),
               "poisson")
  msg <- tryCatch(ares(x, y0, family = "poisson", nthreads = 2),
                  error = function(e) conditionMessage(e))
  expect_match(msg, "at least one y > 0", fixed = TRUE)
})

test_that("family='poisson' accepts y with some positive counts", {
  set.seed(32)
  n <- 100; p <- 3
  x <- matrix(stats::runif(n * p), n, p)
  mu <- exp(0.2 + 0.8 * x[, 1])
  y <- stats::rpois(n, mu)
  # Sanity: at least one non-zero count is virtually certain at mu ~ 1.5.
  stopifnot(any(y > 0))
  fit <- ares(x, y, family = "poisson", nthreads = 2)
  expect_true(is.finite(fit$rss))
  expect_true(length(fit$selected.terms) >= 1)
})

test_that("constant y is rejected for gaussian", {
  set.seed(33)
  n <- 60; p <- 2
  x <- matrix(stats::runif(n * p), n, p)
  yc <- rep(1.5, n)
  expect_error(ares(x, yc, nthreads = 2), "constant")
})

test_that("constant y is rejected for binomial (all zeros)", {
  set.seed(34)
  n <- 60; p <- 2
  x <- matrix(stats::runif(n * p), n, p)
  expect_error(ares(x, rep(0L, n), family = "binomial", nthreads = 2),
               "constant")
})

test_that("constant y is rejected for binomial (all ones)", {
  set.seed(35)
  n <- 60; p <- 2
  x <- matrix(stats::runif(n * p), n, p)
  expect_error(ares(x, rep(1L, n), family = "binomial", nthreads = 2),
               "constant")
})

test_that("constant y is rejected for gamma", {
  set.seed(36)
  n <- 60; p <- 2
  x <- matrix(stats::runif(n * p), n, p)
  expect_error(ares(x, rep(2.5, n), family = "gamma", nthreads = 2),
               "constant")
})

test_that("error message identifies the constant value", {
  set.seed(37)
  n <- 50; p <- 2
  x <- matrix(stats::runif(n * p), n, p)
  msg <- tryCatch(ares(x, rep(3.14, n), nthreads = 2),
                  error = function(e) conditionMessage(e))
  expect_match(msg, "3.14", fixed = TRUE)
})
