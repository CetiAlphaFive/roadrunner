# ols() bagging.

test_that("bagged ols stores a coefficient matrix", {
  fit <- ols(mpg ~ wt + hp, data = mtcars, n.boot = 20L, seed = 42L)
  expect_false(is.null(fit$boot))
  expect_equal(dim(fit$boot$coefficients), c(3L, 20L))
  expect_equal(fit$boot$n.boot, 20L)
})

test_that("bagged ols is reproducible at a fixed seed", {
  f1 <- ols(mpg ~ wt + hp, data = mtcars, n.boot = 25L, seed = 7L)
  f2 <- ols(mpg ~ wt + hp, data = mtcars, n.boot = 25L, seed = 7L)
  expect_identical(f1$boot$coefficients, f2$boot$coefficients)
})

test_that("bagged ols seed does not perturb the user's RNG", {
  set.seed(123L)
  before <- .Random.seed
  ols(mpg ~ wt + hp, data = mtcars, n.boot = 10L, seed = 99L)
  expect_identical(.Random.seed, before)
})

test_that("predict on a bagged ols returns the bag mean", {
  fit <- ols(mpg ~ wt + hp, data = mtcars, n.boot = 50L, seed = 5L)
  nd <- mtcars[1:6, ]
  bc <- fit$boot$coefficients
  xnew <- cbind(1, as.matrix(nd[, c("wt", "hp")]))
  preds <- cbind(drop(xnew %*% fit$coefficients), xnew %*% bc)
  expect_equal(unname(predict(fit, nd)), unname(rowMeans(preds)))
})

test_that("bagged ols predict se.fit returns the bag SD", {
  fit <- ols(mpg ~ wt + hp, data = mtcars, n.boot = 40L, seed = 3L)
  nd <- mtcars[1:6, ]
  p <- predict(fit, nd, se.fit = TRUE)
  expect_length(attr(p, "se.fit"), 6L)
  expect_true(all(attr(p, "se.fit") > 0))
})

test_that("non-bagged ols has no boot slot", {
  fit <- ols(mpg ~ wt + hp, data = mtcars)
  expect_null(fit$boot)
})

test_that("bagged ols composes with weights", {
  set.seed(31)
  w <- runif(nrow(mtcars), 0.5, 2)
  fit <- ols(mpg ~ wt + hp, data = mtcars, weights = w,
             n.boot = 15L, seed = 8L)
  expect_equal(dim(fit$boot$coefficients), c(3L, 15L))
})
