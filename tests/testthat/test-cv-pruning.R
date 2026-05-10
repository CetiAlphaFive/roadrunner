# Tests for pmethod="cv" — CV-based subset selection.

test_that("pmethod='cv' fits without error and returns $cv slot", {
  set.seed(20260510)
  n <- 200; p <- 6
  x <- matrix(stats::runif(n * p), n, p)
  y <- 10 * x[, 1] + 5 * x[, 2] + 3 * sin(pi * x[, 3]) + stats::rnorm(n)
  fit <- ares(x, y, degree = 2, pmethod = "cv", nfold = 5,
              seed.cv = 1, nthreads = 2)
  expect_s3_class(fit, "ares")
  expect_true("cv" %in% names(fit))
  expect_equal(fit$cv$nfold, 5L)
  expect_equal(fit$cv$ncross, 1L)
  expect_true(is.numeric(fit$cv$cv.mse))
  expect_true(fit$cv$size.star >= 1L)
  expect_equal(length(fit$selected.terms), fit$cv$size.star)
  # predict still works
  yhat <- predict(fit, x[1:5, ])
  expect_length(yhat, 5L)
})

test_that("nfold > 0 promotes pmethod to 'cv' automatically", {
  set.seed(20260510)
  x <- matrix(stats::runif(150 * 4), 150, 4)
  y <- 5 * x[, 1] + stats::rnorm(150)
  fit <- ares(x, y, nfold = 5, seed.cv = 7, nthreads = 2)
  expect_equal(fit$pmethod, "cv")
  expect_true("cv" %in% names(fit))
})

test_that("pmethod='cv' is deterministic at fixed seed.cv across nthreads", {
  set.seed(20260510)
  x <- matrix(stats::runif(200 * 5), 200, 5)
  y <- 10 * x[, 1] + 5 * x[, 2] + stats::rnorm(200)
  f1 <- ares(x, y, degree = 1, pmethod = "cv", nfold = 5,
             seed.cv = 42, nthreads = 1)
  f2 <- ares(x, y, degree = 1, pmethod = "cv", nfold = 5,
             seed.cv = 42, nthreads = 2)
  expect_equal(f1$rss, f2$rss, tolerance = 1e-10)
  expect_equal(f1$cv$size.star, f2$cv$size.star)
  expect_equal(unname(f1$coefficients), unname(f2$coefficients),
               tolerance = 1e-10)
})

test_that("ncross > 1 averages over repetitions", {
  set.seed(20260510)
  x <- matrix(stats::runif(200 * 4), 200, 4)
  y <- 5 * x[, 1] + 3 * x[, 2] + stats::rnorm(200)
  f1 <- ares(x, y, pmethod = "cv", nfold = 5, ncross = 1,
             seed.cv = 11, nthreads = 2)
  f3 <- ares(x, y, pmethod = "cv", nfold = 5, ncross = 3,
             seed.cv = 11, nthreads = 2)
  expect_equal(f1$cv$ncross, 1L)
  expect_equal(f3$cv$ncross, 3L)
  expect_s3_class(f3, "ares")
})

test_that("seed.cv preserves the caller's RNG stream", {
  # Build inputs first, OUTSIDE the seed scope, so the only RNG consumer in
  # the seed window is ares() itself.
  x <- matrix(stats::runif(100 * 3), 100, 3)
  y <- stats::runif(100)

  set.seed(123)
  before <- stats::runif(1)

  set.seed(123)
  fit <- ares(x, y, nfold = 3, seed.cv = 777, nthreads = 2)
  after <- stats::runif(1)

  # With seed.cv set, ares saves/restores .Random.seed around the CV path;
  # the next consumer in the user's session must see the same draw as if
  # ares had not been called.
  expect_equal(before, after, tolerance = 1e-12)
})

test_that("seed.cv = NULL does NOT save/restore (uses current RNG)", {
  x <- matrix(stats::runif(80 * 3), 80, 3)
  y <- stats::runif(80)

  set.seed(123)
  fit <- ares(x, y, nfold = 3, nthreads = 2)
  draw_after_with_cv <- stats::runif(1)

  set.seed(123)
  draw_after_no_call <- stats::runif(1)

  # Without seed.cv, the CV partition consumes user RNG, so the next draw
  # must differ from what it would be without calling ares.
  expect_false(isTRUE(all.equal(draw_after_with_cv, draw_after_no_call)))
})

# ---- v0.14: 1-SE rule + per-fold matrix ------------------------------------

test_that("cv.1se = TRUE picks the smallest size within 1 SE of best", {
  set.seed(20260510)
  n <- 250; p <- 5
  x <- matrix(stats::runif(n * p), n, p)
  y <- 10 * sin(pi * x[, 1]) + 5 * x[, 2] + stats::rnorm(n, sd = 1.5)

  fa <- ares(x, y, degree = 2, pmethod = "cv", nfold = 5,
             seed.cv = 1, nthreads = 2, cv.1se = FALSE)
  fb <- ares(x, y, degree = 2, pmethod = "cv", nfold = 5,
             seed.cv = 1, nthreads = 2, cv.1se = TRUE)

  expect_equal(fb$cv$size.argmin, fa$cv$size.argmin)
  expect_lte(fb$cv$size.star, fa$cv$size.argmin)  # 1-SE rule never picks larger
  expect_true(fb$cv$cv.1se)
  expect_false(fa$cv$cv.1se)
})

test_that("cv$cv.mse.mat has the right shape and finite entries", {
  set.seed(20260510)
  x <- matrix(stats::runif(200 * 4), 200, 4)
  y <- 5 * x[, 1] + 3 * x[, 2] + stats::rnorm(200)

  fit <- ares(x, y, pmethod = "cv", nfold = 4, ncross = 2,
              seed.cv = 99, nthreads = 2)
  expect_equal(nrow(fit$cv$cv.mse.mat), 4L * 2L)  # nfold * ncross rows
  expect_equal(ncol(fit$cv$cv.mse.mat), nrow(fit$dirs))
  # At least the smallest sizes should be evaluated by every fold.
  expect_true(any(colSums(!is.na(fit$cv$cv.mse.mat)) == 4L * 2L))
})

test_that("cv$cv.se reports a non-NA SE wherever cv.n > 1", {
  set.seed(20260510)
  x <- matrix(stats::runif(180 * 4), 180, 4)
  y <- 5 * x[, 1] + 3 * x[, 2] + stats::rnorm(180)
  fit <- ares(x, y, pmethod = "cv", nfold = 5, ncross = 2,
              seed.cv = 13, nthreads = 2)
  has_se <- is.finite(fit$cv$cv.se) & fit$cv$cv.n > 1
  has_mean <- is.finite(fit$cv$cv.mse) & fit$cv$cv.n > 1
  expect_true(all(has_se == has_mean))
})

test_that("size.star with 1-SE rule is reproducible at fixed seed.cv", {
  set.seed(20260510)
  x <- matrix(stats::runif(200 * 5), 200, 5)
  y <- 5 * x[, 1] + 2 * x[, 2] + stats::rnorm(200)
  s1 <- ares(x, y, pmethod = "cv", nfold = 5, cv.1se = TRUE,
             seed.cv = 7, nthreads = 1)$cv$size.star
  s2 <- ares(x, y, pmethod = "cv", nfold = 5, cv.1se = TRUE,
             seed.cv = 7, nthreads = 4)$cv$size.star
  expect_equal(s1, s2)
})
