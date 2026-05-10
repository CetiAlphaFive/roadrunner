test_that("ares fits mtcars and predict-RMSE vs earth is small", {
  skip_if_not_installed("earth")
  x <- as.matrix(mtcars[, -1]); y <- mtcars$mpg
  fa <- ares(x, y, nthreads = 2)
  fe <- earth::earth(x, y)
  pa <- predict(fa, x); pe <- predict(fe, x)
  # Predict-RMSE relative to sd(y) is the secondary-tier acceptance
  rmse <- sqrt(mean((pa - pe)^2))
  expect_lt(rmse / stats::sd(y), 0.5)  # generous: small dataset, knot selection diverges
  # GCV is competitive (ares can be slightly worse on tiny n=32)
  expect_lt(fa$gcv, 2 * fe$gcv)
})
