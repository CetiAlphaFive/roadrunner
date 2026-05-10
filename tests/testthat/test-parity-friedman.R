test_that("ares matches earth on Friedman-1 (rel-err < 5e-3 on RSS)", {
  skip_if_not_installed("earth")
  set.seed(20260509)
  n <- 300; p <- 10
  x <- matrix(stats::runif(n * p), n, p)
  y <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
       10 * x[, 4] + 5 * x[, 5] + stats::rnorm(n)
  fa <- ares(x, y, degree = 2, nthreads = 2)
  fe <- earth::earth(x, y, degree = 2)
  rel_err <- abs(fa$rss - fe$rss) / fe$rss
  expect_lt(rel_err, 5e-3)
  # Predicted-y RMSE small relative to sd(y)
  pa <- predict(fa, x); pe <- predict(fe, x)
  rmse <- sqrt(mean((pa - pe)^2))
  expect_lt(rmse / stats::sd(y), 0.1)
})

test_that("ares matches earth on Friedman-1 degree=1 (rel-err < 5e-3)", {
  skip_if_not_installed("earth")
  set.seed(20260509)
  n <- 200; p <- 5
  x <- matrix(stats::runif(n * p), n, p)
  y <- 10 * x[, 1] + 5 * x[, 2] + 3 * x[, 3] + stats::rnorm(n)
  fa <- ares(x, y, degree = 1, nthreads = 2)
  fe <- earth::earth(x, y, degree = 1)
  rel_err <- abs(fa$rss - fe$rss) / fe$rss
  expect_lt(rel_err, 5e-3)
})
