# ols() determinism: byte-identical across thread counts.

test_that("ols fit is byte-identical at nthreads 1 vs N", {
  old <- RcppParallel::defaultNumThreads()
  on.exit(RcppParallel::setThreadOptions(numThreads = old), add = TRUE)

  RcppParallel::setThreadOptions(numThreads = 1L)
  f1 <- ols(mpg ~ wt + hp + disp + qsec, data = mtcars)

  RcppParallel::setThreadOptions(numThreads = 4L)
  f4 <- ols(mpg ~ wt + hp + disp + qsec, data = mtcars)

  expect_identical(f1$coefficients, f4$coefficients)
  expect_identical(f1$vcov, f4$vcov)
  expect_identical(f1$residuals, f4$residuals)
  expect_identical(f1$hatvalues, f4$hatvalues)
})

test_that("bagged ols is byte-identical across thread counts", {
  old <- RcppParallel::defaultNumThreads()
  on.exit(RcppParallel::setThreadOptions(numThreads = old), add = TRUE)

  RcppParallel::setThreadOptions(numThreads = 1L)
  f1 <- ols(mpg ~ wt + hp, data = mtcars, n.boot = 30L, seed = 17L)

  RcppParallel::setThreadOptions(numThreads = 4L)
  f4 <- ols(mpg ~ wt + hp, data = mtcars, n.boot = 30L, seed = 17L)

  expect_identical(f1$boot$coefficients, f4$boot$coefficients)
})

test_that("ols robust vcov is byte-identical across thread counts", {
  old <- RcppParallel::defaultNumThreads()
  on.exit(RcppParallel::setThreadOptions(numThreads = old), add = TRUE)

  fit <- ols(mpg ~ wt + hp + disp, data = mtcars)

  RcppParallel::setThreadOptions(numThreads = 1L)
  v1 <- roadrunner:::.ols_vcov(fit, "HC3")

  RcppParallel::setThreadOptions(numThreads = 4L)
  v4 <- roadrunner:::.ols_vcov(fit, "HC3")

  expect_identical(v1, v4)
})
