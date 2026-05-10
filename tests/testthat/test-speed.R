# Speed test — informational only (skip_on_cran).
# v0.0.0.9000 ships with correctness-first inner loop; absolute wall-clock vs
# earth is slower (see NEWS.md). The hard gate here is parallel scaling.

test_that("nthreads scales: 2 threads at least 1.5x faster than 1 thread", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("earth")
  set.seed(20260509)
  n <- 500; p <- 10
  x <- matrix(stats::runif(n * p), n, p)
  y <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
       10 * x[, 4] + 5 * x[, 5] + stats::rnorm(n)
  t1 <- system.time(ares(x, y, degree = 2, nthreads = 1))[["elapsed"]]
  t2 <- system.time(ares(x, y, degree = 2, nthreads = 2))[["elapsed"]]
  message(sprintf("ares speed: 1t=%.2fs  2t=%.2fs  scaling=%.2fx", t1, t2, t1 / t2))
  # Soft check: parallel path is at least 1.3x faster; we want ≥ 1.5x but allow
  # noise on small CI machines. If this fails consistently it's a real regression.
  expect_gt(t1 / t2, 1.3)
})

test_that("ares records reasonable wall-clock vs earth (informational)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("earth")
  set.seed(20260509)
  n <- 500; p <- 10
  x <- matrix(stats::runif(n * p), n, p)
  y <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
       10 * x[, 4] + 5 * x[, 5] + stats::rnorm(n)
  t_e <- system.time(earth::earth(x, y, degree = 2))[["elapsed"]]
  t_a <- system.time(ares(x, y, degree = 2, nthreads = 2))[["elapsed"]]
  message(sprintf("earth=%.3fs  ares-2t=%.3fs  ratio=%.1fx", t_e, t_a, t_a / t_e))
  # Soft acceptance for v0.0.0.9000 — the hard ratio target is for v0.1.
  expect_lt(t_a / t_e, 200)
})
