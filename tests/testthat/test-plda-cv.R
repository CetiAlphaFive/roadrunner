test_that(".plda_cv selects a lambda and reports a grid", {
  set.seed(5)
  x <- as.matrix(iris[, 1:4]); yint <- as.integer(iris$Species)
  cv <- roadrunner:::.plda_cv(x, yint, G = 3L, K = 2L, penalty = "L1",
                              pen_code = 0L, lam2 = 0, nfold = 5L,
                              lambda_grid = NULL, maxit = 100L, tol = 1e-6)
  expect_true(is.numeric(cv$lambda) && cv$lambda > 0)
  expect_true(cv$K >= 1L && cv$K <= 2L)
  expect_true(length(cv$grid) >= 5L)
  expect_equal(length(cv$errors), length(cv$grid))
})

test_that("plda autotune fit predicts iris well", {
  set.seed(6)
  fit <- plda(Species ~ ., data = iris)   # autotune = TRUE default
  acc <- mean(predict(fit, iris) == iris$Species)
  expect_gt(acc, 0.9)
})

test_that(".plda_cv does not modify the caller RNG stream", {
  set.seed(42); ref <- runif(1)
  set.seed(42)
  roadrunner:::.plda_cv(as.matrix(iris[, 1:4]), as.integer(iris$Species),
                        G = 3L, K = 2L, penalty = "L1", pen_code = 0L, lam2 = 0,
                        nfold = 3L, lambda_grid = NULL, maxit = 50L, tol = 1e-4)
  expect_equal(runif(1), ref)
})

test_that("plda autotune honors an explicit K", {
  set.seed(7)
  fit <- plda(Species ~ ., data = iris, K = 1)   # autotune lambda, K fixed at 1
  expect_equal(fit$K, 1L)
})

test_that(".plda_lambda_grid spans a real data-driven range", {
  g <- roadrunner:::.plda_lambda_grid(as.matrix(iris[, 1:4]),
                                      as.integer(iris$Species), G = 3L)
  expect_length(g, 12L)
  expect_gt(max(g), 0.1)          # not the old collapsed 1e-8 grid
  expect_true(all(diff(g) > 0))   # ascending
})
