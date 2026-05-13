# BUG-012 regression test (v0.0.0.9032).
#
# Formula-path fits with derived terms (I(x^2), poly(x, 2), log(x + 10),
# scale(x), splines::bs(x), ...) used to fit successfully but predict()
# would error with "newdata is missing columns: I(x^2)" -- predict.ares
# was looking for the *expanded* column name as a literal column of
# newdata, not re-evaluating the original terms object on newdata.
#
# Fix: when object$terms is non-null (formula path), use
# model.matrix(delete.response(terms), newdata, xlev = object$xlevels)
# to rebuild the design with derived terms re-evaluated. Falls back to
# the prior column-lookup path when no terms object is stored (matrix
# interface).

test_that("predict works for I(x^2) derived term", {
  set.seed(1)
  df <- data.frame(x = stats::rnorm(80))
  df$y <- df$x + df$x^2 + stats::rnorm(80, sd = 0.3)
  fit <- ares(y ~ x + I(x^2), data = df, nthreads = 2)
  df_new <- data.frame(x = seq(-2, 2, length.out = 20))
  p <- predict(fit, df_new)
  expect_equal(length(p), nrow(df_new))
  expect_true(all(is.finite(p)))
})

test_that("predict works for log(x + 10) derived term", {
  set.seed(2)
  df <- data.frame(x = stats::runif(60, 0, 5))
  df$y <- log(df$x + 10) + stats::rnorm(60, sd = 0.2)
  fit <- ares(y ~ log(x + 10), data = df, nthreads = 2)
  df_new <- data.frame(x = stats::runif(15, 0, 5))
  p <- predict(fit, df_new)
  expect_equal(length(p), nrow(df_new))
  expect_true(all(is.finite(p)))
})

test_that("predict works for poly(x, 2)", {
  set.seed(3)
  df <- data.frame(x = stats::rnorm(80))
  df$y <- df$x + df$x^2 + stats::rnorm(80, sd = 0.3)
  fit <- ares(y ~ poly(x, 2), data = df, nthreads = 2)
  df_new <- data.frame(x = seq(-2, 2, length.out = 12))
  p <- predict(fit, df_new)
  expect_equal(length(p), nrow(df_new))
  expect_true(all(is.finite(p)))
})

test_that("predict works for scale(x)", {
  set.seed(4)
  df <- data.frame(x = stats::rnorm(70))
  df$y <- df$x + stats::rnorm(70, sd = 0.4)
  fit <- ares(y ~ scale(x), data = df, nthreads = 2)
  df_new <- data.frame(x = stats::rnorm(20))
  p <- predict(fit, df_new)
  expect_equal(length(p), nrow(df_new))
  expect_true(all(is.finite(p)))
})

test_that("predict reproduces training fitted values on same data", {
  set.seed(5)
  df <- data.frame(x = stats::rnorm(60))
  df$y <- df$x + df$x^2 + stats::rnorm(60, sd = 0.3)
  fit <- ares(y ~ x + I(x^2), data = df, nthreads = 2)
  p_new  <- as.numeric(predict(fit, df))
  expect_equal(p_new, as.numeric(fit$fitted.values), tolerance = 1e-10)
})
