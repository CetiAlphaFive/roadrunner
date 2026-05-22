# logreg() determinism: byte-identical across thread counts.

test_that("logreg fit is byte-identical at nthreads 1 vs N", {
  old <- RcppParallel::defaultNumThreads()
  on.exit(RcppParallel::setThreadOptions(numThreads = old), add = TRUE)
  df <- data.frame(y = as.integer(mtcars$vs),
                   mtcars[c("wt", "drat", "gear")])

  RcppParallel::setThreadOptions(numThreads = 1L)
  f1 <- logreg(y ~ wt + drat + gear, data = df)

  RcppParallel::setThreadOptions(numThreads = 4L)
  f4 <- logreg(y ~ wt + drat + gear, data = df)

  expect_identical(f1$coefficients, f4$coefficients)
  expect_identical(f1$vcov, f4$vcov)
  expect_identical(f1$fitted.values, f4$fitted.values)
  expect_identical(f1$linear.predictors, f4$linear.predictors)
  expect_identical(f1$hatvalues, f4$hatvalues)
  expect_identical(f1$deviance, f4$deviance)
})

test_that("bagged logreg is byte-identical across thread counts", {
  old <- RcppParallel::defaultNumThreads()
  on.exit(RcppParallel::setThreadOptions(numThreads = old), add = TRUE)
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])

  RcppParallel::setThreadOptions(numThreads = 1L)
  f1 <- logreg(y ~ wt + hp, data = df, n.boot = 30L, seed = 17L)

  RcppParallel::setThreadOptions(numThreads = 4L)
  f4 <- logreg(y ~ wt + hp, data = df, n.boot = 30L, seed = 17L)

  expect_identical(f1$boot$coefficients, f4$boot$coefficients)
})

test_that("logreg robust vcov is byte-identical across thread counts", {
  old <- RcppParallel::defaultNumThreads()
  on.exit(RcppParallel::setThreadOptions(numThreads = old), add = TRUE)
  df <- data.frame(y = as.integer(mtcars$vs), mtcars[c("wt", "drat", "gear")])
  fit <- logreg(y ~ wt + drat + gear, data = df)

  RcppParallel::setThreadOptions(numThreads = 1L)
  v1 <- roadrunner:::.logreg_vcov(fit, "HC3")

  RcppParallel::setThreadOptions(numThreads = 4L)
  v4 <- roadrunner:::.logreg_vcov(fit, "HC3")

  expect_identical(v1, v4)
})
