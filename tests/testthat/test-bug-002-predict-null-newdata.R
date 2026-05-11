# BUG-002 regression test (v0.0.0.9029).
#
# `predict(bagged_fit, newdata = NULL)` and `predict(bagged_fit, x_train)`
# must return identical results (both should be the bag mean). Prior to
# v0.0.0.9029 the NULL branch silently short-circuited to the central
# fit's $fitted.values, diverging from the explicit-newdata path.

test_that("predict(bagged_fit) returns the bag mean for gaussian fits", {
  set.seed(70)
  n <- 200; p <- 4
  x <- matrix(stats::runif(n * p), n, p)
  y <- 3 * x[, 1] + 2 * x[, 2]^2 + stats::rnorm(n, sd = 0.4)
  fb <- ares(x, y, n.boot = 5, seed.cv = 5L, nthreads = 2)
  pn   <- predict(fb)
  pxin <- predict(fb, x)
  expect_equal(as.numeric(pn), as.numeric(pxin), tolerance = 1e-12)
})

test_that("predict(bagged_fit) on binomial returns bag mean (response + link)", {
  set.seed(71)
  n <- 200; p <- 4
  x <- matrix(stats::runif(n * p), n, p)
  yb <- as.integer(stats::plogis(2 * x[, 1] - 0.5) > 0.5)
  fb <- ares(x, yb, family = "binomial", n.boot = 4, seed.cv = 7L,
             nthreads = 2)
  expect_equal(as.numeric(predict(fb)),
               as.numeric(predict(fb, x)),
               tolerance = 1e-12)
  expect_equal(as.numeric(predict(fb, type = "link")),
               as.numeric(predict(fb, x, type = "link")),
               tolerance = 1e-12)
})

test_that("predict(bagged_fit, se.fit=TRUE) returns NA-free bag SD on NULL newdata", {
  set.seed(72)
  n <- 200; p <- 3
  x <- matrix(stats::runif(n * p), n, p)
  y <- 2 * x[, 1] + x[, 2]^2 + stats::rnorm(n, sd = 0.3)
  fb <- ares(x, y, n.boot = 5, seed.cv = 11L, nthreads = 2)
  p_null <- predict(fb, se.fit = TRUE)
  p_x    <- predict(fb, x, se.fit = TRUE)
  expect_equal(attr(p_null, "sd"), attr(p_x, "sd"), tolerance = 1e-12)
  expect_true(all(!is.na(attr(p_null, "sd"))))
})

test_that("non-bagged fit: predict(fit) is unchanged (central fitted values)", {
  set.seed(73)
  n <- 80; p <- 3
  x <- matrix(stats::runif(n * p), n, p)
  y <- x[, 1] + x[, 2]^2 + stats::rnorm(n, sd = 0.2)
  f <- ares(x, y, nthreads = 2)
  expect_equal(as.numeric(predict(f)),
               as.numeric(f$fitted.values), tolerance = 0)
})
