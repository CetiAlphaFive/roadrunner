test_that("plda.default returns a well-formed plda object", {
  set.seed(3)
  x <- as.matrix(iris[, 1:4]); y <- iris$Species
  fit <- plda(x, y, K = 2, lambda = 0.1, autotune = FALSE)
  expect_s3_class(fit, "plda")
  expect_equal(dim(fit$discrim), c(4L, 2L))
  expect_equal(fit$classes, levels(y))
  expect_equal(fit$penalty, "L1")
})

test_that("plda rejects K greater than G-1", {
  x <- as.matrix(iris[, 1:4]); y <- iris$Species
  expect_error(plda(x, y, K = 3, lambda = 0.1, autotune = FALSE), "K")
})

test_that("plda.formula equals plda.default", {
  f <- plda(Species ~ ., data = iris, K = 2, lambda = 0.1, autotune = FALSE)
  d <- plda(as.matrix(iris[, 1:4]), iris$Species, K = 2, lambda = 0.1, autotune = FALSE)
  expect_equal(unname(f$discrim), unname(d$discrim), tolerance = 1e-10)
  expect_s3_class(f, "plda")
})
