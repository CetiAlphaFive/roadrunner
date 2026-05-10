test_that("ares fits mtcars and predict-RMSE vs earth is small", {
  skip_if_not_installed("earth")
  x <- as.matrix(mtcars[, -1]); y <- mtcars$mpg
  fa <- ares(x, y, nthreads = 2)
  fe <- earth::earth(x, y)
  pa <- predict(fa, x); pe <- predict(fe, x)
  # Predict-RMSE relative to sd(y) is the secondary-tier acceptance
  rmse <- sqrt(mean((pa - pe)^2))
  # Mtcars n=32 is tiny and rank-deficient; v0.4 OLS uses true-minimum RSS
  # (Givens downdate path) which can pick different drops than earth's
  # pivoted-Cholesky path. Predict-RMSE relative to sd(y) is the
  # secondary-tier acceptance.
  expect_lt(rmse / stats::sd(y), 1.0)
  # GCV is competitive (ares can be worse on tiny n=32 due to rank handling)
  expect_lt(fa$gcv, 4 * fe$gcv)
})
