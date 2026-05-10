test_that("constant column produces a warning and is dropped", {
  x <- cbind(stats::runif(100), rep(1, 100), stats::runif(100))
  y <- 2 * x[, 1] + stats::rnorm(100)
  expect_warning(fit <- ares(x, y, nthreads = 2), "constant column")
  expect_equal(ncol(fit$dirs), 2L)  # one dropped
})

test_that("single-predictor case fits without error", {
  x <- matrix(stats::runif(100), ncol = 1)
  y <- pmax(0, x - 0.5) * 5 + stats::rnorm(100)
  fit <- ares(x, y, nthreads = 2)
  expect_s3_class(fit, "ares")
  expect_equal(ncol(fit$dirs), 1L)
})

test_that("degree=2 produces RSS no worse than degree=1", {
  set.seed(42)
  n <- 200; p <- 5
  x <- matrix(stats::runif(n * p), n, p)
  y <- 5 * pmax(0, x[, 1] - 0.5) * pmax(0, x[, 2] - 0.3) + stats::rnorm(n)
  f1 <- ares(x, y, degree = 1, nthreads = 2)
  f2 <- ares(x, y, degree = 2, nthreads = 2)
  expect_lte(f2$rss, f1$rss + 1e-9)
})

test_that("input validation: bad x", {
  expect_error(ares("hello", 1:5), "x must be")
  expect_error(ares(matrix(NA_real_, 5, 2), 1:5), "NA")
})

test_that("input validation: bad y", {
  expect_error(ares(matrix(1, 5, 2), c(1, 2, NA, 4, 5)), "NA")
  expect_error(ares(matrix(1, 5, 2), 1:4), "length")
})

test_that("input validation: degree", {
  expect_error(ares(matrix(stats::runif(20), 10, 2), 1:10, degree = 0), "degree")
})
