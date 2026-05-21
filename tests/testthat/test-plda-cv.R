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
