# Tests for autotune (Phase 2 — v0.15+).

test_that("autotune = TRUE returns $autotune slot with grid + winner", {
  set.seed(20260510)
  n <- 200; p <- 5
  x <- matrix(stats::runif(n * p), n, p)
  y <- 10 * sin(pi * x[, 1]) + 3 * x[, 2] + stats::rnorm(n, sd = 1.5)
  fit <- ares(x, y, autotune = TRUE, seed.cv = 1, nthreads = 2)
  expect_true("autotune" %in% names(fit))
  expect_true(is.data.frame(fit$autotune$grid))
  expect_true(all(c("degree", "penalty", "cv_mse") %in% names(fit$autotune$grid)))
  # grid must include the chosen (degree, penalty)
  hit <- with(fit$autotune$grid,
              degree == fit$autotune$degree &
              abs(penalty - fit$autotune$penalty) < 1e-12)
  expect_true(any(hit))
  expect_true(fit$autotune$cv_mse <= max(fit$autotune$grid$cv_mse))
})

test_that("autotune chooses degree 2 on a clearly-interaction DGP", {
  set.seed(20260510)
  n <- 250; p <- 5
  x <- matrix(stats::runif(n * p), n, p)
  # Strong x1*x2 interaction; deg=1 cannot capture it.
  y <- 10 * sin(pi * x[, 1] * x[, 2]) +
       20 * (x[, 3] - 0.5)^2 + stats::rnorm(n)
  fit <- ares(x, y, autotune = TRUE, seed.cv = 7, nthreads = 2)
  expect_gte(fit$autotune$degree, 2L)
})

test_that("autotune is deterministic at fixed seed.cv across nthreads", {
  set.seed(20260510)
  x <- matrix(stats::runif(200 * 5), 200, 5)
  y <- 5 * x[, 1] + stats::rnorm(200)
  s1 <- ares(x, y, autotune = TRUE, seed.cv = 42, nthreads = 1)
  s2 <- ares(x, y, autotune = TRUE, seed.cv = 42, nthreads = 4)
  expect_equal(s1$rss, s2$rss, tolerance = 1e-10)
  expect_equal(s1$autotune$degree, s2$autotune$degree)
  expect_equal(s1$autotune$penalty, s2$autotune$penalty)
  expect_equal(unname(s1$coefficients), unname(s2$coefficients),
               tolerance = 1e-10)
})

test_that("autotune is non-blocking on small n / few features", {
  set.seed(20260510)
  x <- matrix(stats::runif(80 * 3), 80, 3)
  y <- 5 * x[, 1] + stats::rnorm(80)
  fit <- ares(x, y, autotune = TRUE, seed.cv = 1, nthreads = 2)
  expect_s3_class(fit, "ares")
  # Should not include degree=3 (nk_eff is small for p=3).
  expect_true(all(fit$autotune$grid$degree <= 2L))
})

test_that("autotune predict round-trips and pmethod is reported", {
  set.seed(20260510)
  x <- matrix(stats::runif(150 * 4), 150, 4)
  y <- 5 * x[, 1] + 3 * x[, 2] + stats::rnorm(150)
  fit <- ares(x, y, autotune = TRUE, seed.cv = 11, nthreads = 2)
  expect_equal(fit$pmethod, "backward")
  yhat <- predict(fit, x[1:5, ])
  expect_length(yhat, 5L)
})

# ---- v0.16: nk grid + successive halving -----------------------------------

test_that("autotune grid includes nk multipliers when nk_eff is large enough", {
  set.seed(20260510)
  x <- matrix(stats::runif(200 * 6), 200, 6)
  y <- 5 * x[, 1] + 3 * x[, 2] + stats::rnorm(200)
  fit <- ares(x, y, autotune = TRUE, seed.cv = 1, nthreads = 2)
  expect_true("nk" %in% names(fit$autotune$grid))
  # Should sweep at least 2 nk values when default nk is small enough that
  # 2x and 4x do not collapse into the cap of 200.
  expect_gte(length(unique(fit$autotune$grid$nk)), 2L)
  # Winner must report nk that appears in the grid.
  expect_true(fit$autotune$nk %in% fit$autotune$grid$nk)
})

test_that("successive halving eliminates clearly bad cells after fold 1", {
  set.seed(20260510)
  # Strong interaction DGP: deg=1 cells should die fast.
  n <- 250; p <- 6
  x <- matrix(stats::runif(n * p), n, p)
  y <- 10 * sin(pi * x[, 1] * x[, 2]) +
       20 * (x[, 3] - 0.5)^2 + stats::rnorm(n)
  fit <- ares(x, y, autotune = TRUE, seed.cv = 1, nthreads = 2)
  # All degree-1 cells should be eliminated on this DGP.
  d1_elim <- with(fit$autotune$grid, eliminated[degree == 1L])
  expect_true(any(d1_elim))
  # Winner is degree >= 2.
  expect_gte(fit$autotune$degree, 2L)
})

test_that("autotune determinism across nthreads with v0.16 grid", {
  set.seed(20260510)
  x <- matrix(stats::runif(200 * 5), 200, 5)
  y <- 5 * x[, 1] + 3 * x[, 2] + stats::rnorm(200)
  s1 <- ares(x, y, autotune = TRUE, seed.cv = 42, nthreads = 1)
  s2 <- ares(x, y, autotune = TRUE, seed.cv = 42, nthreads = 4)
  expect_identical(s1$autotune$nk, s2$autotune$nk)
  expect_identical(s1$autotune$grid$cv_mse, s2$autotune$grid$cv_mse)
  expect_identical(s1$autotune$grid$eliminated, s2$autotune$grid$eliminated)
  expect_equal(s1$rss, s2$rss, tolerance = 1e-10)
})
