# Tests for varmod / prediction intervals.
#
# Quick:
#   - varmod = "none" (default): no $varmod slot stored.
#   - varmod = "const" / "lm" populate $varmod with the right shape.
#   - predict(interval = "pint") returns a (fit, lwr, upr) matrix.
#   - PI matrix shape and ordering (lwr <= fit <= upr).
#   - varmod silently ignored for non-gaussian families.
#   - predict(interval = "pint") errors when no varmod was stored.
# Heavy (skip_if_quick):
#   - empirical coverage of 95% PI on a clean gaussian DGP at n=1000
#     lies inside [0.90, 0.99].

.gauss_dgp <- function(n = 500, p = 5, seed = 1, sigma = 1) {
  set.seed(seed)
  x <- matrix(stats::runif(n * p), n, p)
  f <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
       10 * x[, 4] + 5 * x[, 5]
  y <- f + stats::rnorm(n, sd = sigma)
  list(x = x, y = y, signal = f, sigma = sigma)
}

test_that("varmod = 'none' (default) stores no variance model", {
  d <- .gauss_dgp()
  f <- ares(d$x, d$y, nthreads = 2)
  expect_null(f$varmod)
})

test_that("varmod = 'const' stores sigma_hat and df", {
  d <- .gauss_dgp(seed = 2)
  f <- ares(d$x, d$y, varmod = "const", nthreads = 2)
  expect_identical(f$varmod$type, "const")
  expect_true(is.finite(f$varmod$sigma_hat) && f$varmod$sigma_hat > 0)
  expect_true(f$varmod$df > 0)
})

test_that("varmod = 'lm' stores intercept + slope", {
  d <- .gauss_dgp(seed = 3)
  f <- ares(d$x, d$y, varmod = "lm", nthreads = 2)
  expect_identical(f$varmod$type, "lm")
  expect_true(is.finite(f$varmod$intercept))
  expect_true(is.finite(f$varmod$slope))
  expect_equal(f$varmod$scale, sqrt(pi / 2), tolerance = 1e-12)
})

test_that("predict(interval = 'pint') returns (fit, lwr, upr) matrix", {
  d <- .gauss_dgp(seed = 4)
  f <- ares(d$x, d$y, varmod = "const", nthreads = 2)
  pm <- predict(f, d$x, interval = "pint")
  expect_true(is.matrix(pm))
  expect_identical(colnames(pm), c("fit", "lwr", "upr"))
  expect_equal(nrow(pm), nrow(d$x))
  expect_true(all(pm[, "lwr"] <= pm[, "fit"]))
  expect_true(all(pm[, "fit"] <= pm[, "upr"]))
})

test_that("predict(interval = 'pint') with newdata=NULL also works", {
  d <- .gauss_dgp(seed = 5)
  f <- ares(d$x, d$y, varmod = "const", nthreads = 2)
  pm <- predict(f, interval = "pint")
  expect_true(is.matrix(pm))
  expect_equal(nrow(pm), length(d$y))
})

test_that("predict(interval = 'pint') errors when no varmod", {
  d <- .gauss_dgp(seed = 6)
  f <- ares(d$x, d$y, nthreads = 2)
  expect_error(predict(f, d$x, interval = "pint"),
               "varmod")
})

test_that("varmod silently skipped for non-gaussian families", {
  set.seed(7)
  n <- 200; p <- 4
  x <- matrix(runif(n * p), n, p)
  y <- as.integer(plogis(x[, 1] - x[, 2]) > 0.5)
  fb <- ares(x, y, family = "binomial", varmod = "const", nthreads = 2)
  expect_null(fb$varmod)
  expect_error(predict(fb, x, interval = "pint"),
               "gaussian")
})

test_that("varmod = 'lm' captures heteroscedasticity directionally", {
  # Construct a DGP where sd grows with the signal so the lm slope is
  # positive. The lm-varmod's stored slope should be > 0.
  set.seed(8)
  n <- 800; p <- 3
  x <- matrix(runif(n * p), n, p)
  f <- 5 * x[, 1]
  y <- f + rnorm(n, sd = 0.5 + 2 * x[, 1])  # sd grows with f
  fit <- ares(x, y, varmod = "lm", nthreads = 2)
  # Slope is for |resid| ~ yhat; positive when residual scale grows with yhat.
  expect_gt(fit$varmod$slope, 0)
})

test_that("PI empirical coverage on clean DGP lies in [0.90, 0.99] (heavy)", {
  skip_if_quick()
  set.seed(9)
  n_train <- 1000
  n_test  <- 1000
  d_tr <- .gauss_dgp(n = n_train, seed = 9, sigma = 1)
  d_te <- .gauss_dgp(n = n_test,  seed = 109, sigma = 1)
  f <- ares(d_tr$x, d_tr$y, varmod = "const", nthreads = 2)
  pm <- predict(f, d_te$x, interval = "pint", level = 0.95)
  hit <- d_te$y >= pm[, "lwr"] & d_te$y <= pm[, "upr"]
  cov <- mean(hit)
  expect_gt(cov, 0.90)
  expect_lt(cov, 0.99)
})

test_that("PI composes with weights (sigma estimated under WLS) (heavy)", {
  skip_if_quick()
  set.seed(10)
  n <- 600
  x <- matrix(runif(n * 4), n, 4)
  f <- 3 * x[, 1] + 2 * x[, 2]
  w <- 1 / (0.5 + x[, 1])
  y <- f + rnorm(n)
  fit <- ares(x, y, weights = w, varmod = "const", nthreads = 2)
  expect_true(is.finite(fit$varmod$sigma_hat))
  pm <- predict(fit, x, interval = "pint")
  expect_true(all(pm[, "lwr"] <= pm[, "upr"]))
})
