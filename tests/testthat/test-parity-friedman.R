# Parity vs earth on RSS. v0.4 switched the backward pass from per-trial
# Householder rebuild (with the rank-clamped β path inside `ols_qr`) to a
# Givens downdate of an upfront QR factor. The downdate computes the *true*
# OLS minimum RSS rather than the rank-clamped pseudo-minimum, so on
# rank-deficient designs the two solvers can diverge by O(1%) of RSS even
# though the model selected is the same. We loosen the threshold from
# 5e-3 to 2e-2 to reflect this — full earth-parity in the strict sense
# requires reproducing earth's pivoted-Cholesky rank handling, which is
# v0.5+ work.
test_that("ares matches earth on Friedman-1 (rel-err < 2e-2 on RSS)", {
  skip_if_not_installed("earth")
  set.seed(20260509)
  n <- 300; p <- 10
  x <- matrix(stats::runif(n * p), n, p)
  y <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
       10 * x[, 4] + 5 * x[, 5] + stats::rnorm(n)
  fa <- ares(x, y, degree = 2, nthreads = 2)
  fe <- earth::earth(x, y, degree = 2)
  rel_err <- abs(fa$rss - fe$rss) / fe$rss
  expect_lt(rel_err, 2e-2)
  # Predicted-y RMSE small relative to sd(y)
  pa <- predict(fa, x); pe <- predict(fe, x)
  rmse <- sqrt(mean((pa - pe)^2))
  expect_lt(rmse / stats::sd(y), 0.1)
})

test_that("ares matches earth on Friedman-1 degree=1 (rel-err < 2e-2)", {
  skip_if_not_installed("earth")
  set.seed(20260509)
  n <- 200; p <- 5
  x <- matrix(stats::runif(n * p), n, p)
  y <- 10 * x[, 1] + 5 * x[, 2] + 3 * x[, 3] + stats::rnorm(n)
  fa <- ares(x, y, degree = 1, nthreads = 2)
  fe <- earth::earth(x, y, degree = 1)
  rel_err <- abs(fa$rss - fe$rss) / fe$rss
  expect_lt(rel_err, 2e-2)
})
