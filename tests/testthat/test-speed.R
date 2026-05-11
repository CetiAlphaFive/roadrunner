# Speed test — informational only (skip_on_cran).
# Thread-scaling is checked at a problem size large enough that per-thread work
# dominates TBB scheduler overhead. With v0.1 fast-LS the per-pair inner loop
# is much cheaper, so small n can be scheduler-bound.

test_that("nthreads scales: 2 threads >= 1 thread (does not regress)", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("earth")
  set.seed(20260509)
  n <- 800; p <- 10
  x <- matrix(stats::runif(n * p), n, p)
  y <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
       10 * x[, 4] + 5 * x[, 5] + stats::rnorm(n)
  t1 <- system.time(ares(x, y, degree = 2, nthreads = 1))[["elapsed"]]
  t2 <- system.time(ares(x, y, degree = 2, nthreads = 2))[["elapsed"]]
  message(sprintf("ares speed: 1t=%.2fs  2t=%.2fs  scaling=%.2fx", t1, t2, t1 / t2))
  # Post v0.1 (fast-LS): per-pair inner loop dropped from O(K*n*Mq) to
  # O((n+K)*Mq), so the parallelizable region is much smaller and serial parts
  # (build_Q, ols_qr) now constrain scaling via Amdahl. We assert only that
  # 2t does not regress vs 1t; tighter scaling is informational, see inst/sims/.
  expect_gte(t1 / t2, 1.0)
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
