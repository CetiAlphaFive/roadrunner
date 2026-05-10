test_that("nthreads=1 and nthreads=2 produce identical fits (Friedman-1 deg=1)", {
  set.seed(20260509)
  n <- 200; p <- 5
  x <- matrix(stats::runif(n * p), n, p)
  y <- 10 * x[, 1] + 5 * x[, 2] + stats::rnorm(n)
  f1 <- ares(x, y, degree = 1, nthreads = 1)
  f2 <- ares(x, y, degree = 1, nthreads = 2)
  expect_equal(f1$selected.terms, f2$selected.terms)
  expect_equal(f1$rss, f2$rss, tolerance = 1e-10)
  expect_equal(unname(f1$coefficients), unname(f2$coefficients), tolerance = 1e-10)
})

test_that("nthreads=1 and nthreads=2 produce identical fits (Friedman-1 deg=2)", {
  set.seed(20260509)
  n <- 150; p <- 5
  x <- matrix(stats::runif(n * p), n, p)
  y <- 10 * sin(pi * x[, 1] * x[, 2]) + 5 * x[, 3] + stats::rnorm(n)
  f1 <- ares(x, y, degree = 2, nthreads = 1)
  f2 <- ares(x, y, degree = 2, nthreads = 2)
  expect_equal(f1$selected.terms, f2$selected.terms)
  expect_equal(f1$rss, f2$rss, tolerance = 1e-10)
})
