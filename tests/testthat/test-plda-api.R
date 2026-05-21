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

# ---- Regression tests for code-review fixes ----------------------------------

test_that("plda errors when autotune=FALSE and lambda=NULL", {
  x <- as.matrix(iris[, 1:4]); y <- iris$Species
  expect_error(plda(x, y, autotune = FALSE), "lambda")
})

test_that("plda errors on negative lambda", {
  x <- as.matrix(iris[, 1:4]); y <- iris$Species
  expect_error(plda(x, y, lambda = -0.1, autotune = FALSE),
               "lambda.*positive|positive.*lambda")
})

test_that("plda errors on zero lambda", {
  x <- as.matrix(iris[, 1:4]); y <- iris$Species
  expect_error(plda(x, y, lambda = 0, autotune = FALSE),
               "lambda.*positive|positive.*lambda")
})

test_that("plda errors when x contains NA", {
  x <- as.matrix(iris[, 1:4]); y <- iris$Species
  x[1, 1] <- NA
  expect_error(plda(x, y, lambda = 0.1, autotune = FALSE), "NA")
})

test_that("plda.formula works when data is omitted and vars are in calling env", {
  x_local <- as.matrix(iris[, 1:4])
  y_local <- iris$Species
  f <- plda(y_local ~ x_local, lambda = 0.1, autotune = FALSE)
  expect_s3_class(f, "plda")
})

test_that("plda print and summary work", {
  fit <- plda(Species ~ ., data = iris, K = 2, lambda = 0.1, autotune = FALSE)
  expect_output(print(fit), "Penalized LDA")
  s <- summary(fit)
  expect_s3_class(s, "summary.plda")
  expect_output(print(s), "nonzero")
})

test_that("predict.plda errors on matrix newdata with wrong column count", {
  fit <- plda(as.matrix(iris[, 1:4]), iris$Species,
              K = 2, lambda = 0.1, autotune = FALSE)
  bad <- as.matrix(iris[, 1:3])
  expect_snapshot(error = TRUE, predict(fit, bad))
})

test_that("predict.plda returns class / posterior / projection", {
  set.seed(11)
  x <- as.matrix(iris[, 1:4]); y <- iris$Species
  fit <- plda(x, y, K = 2, lambda = 0.05, autotune = FALSE)
  cl <- predict(fit, x)
  expect_s3_class(cl, "factor")
  expect_equal(levels(cl), levels(y))
  expect_gt(mean(cl == y), 0.9)
  po <- predict(fit, x, type = "posterior")
  expect_equal(dim(po), c(150L, 3L))
  expect_equal(unname(rowSums(po)), rep(1, 150), tolerance = 1e-8)
  pr <- predict(fit, x, type = "projection")
  expect_equal(dim(pr), c(150L, 2L))
})

test_that("plda handles edge cases", {
  set.seed(12)
  # p >> n
  x <- matrix(rnorm(20 * 200), 20, 200); y <- factor(rep(1:2, each = 10))
  expect_s3_class(plda(x, y, K = 1, lambda = 0.2, autotune = FALSE), "plda")
  # huge lambda -> all-zero discriminant
  f0 <- plda(x, y, K = 1, lambda = 1e6, autotune = FALSE)
  expect_true(all(f0$discrim[, 1] == 0))
  # constant feature does not crash
  xc <- x; xc[, 1] <- 3
  expect_s3_class(plda(xc, y, K = 1, lambda = 0.2, autotune = FALSE), "plda")
  # OOV / missing predictor column in predict errors clearly
  fit <- plda(Species ~ ., data = iris, K = 2, lambda = 0.1, autotune = FALSE)
  expect_error(predict(fit, iris[1:5, 1:3]), "not found")
})
