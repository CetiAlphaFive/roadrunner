# logreg() bagging.

test_that("bagged logreg stores a coefficient matrix", {
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  fit <- logreg(y ~ wt + hp, data = df, n.boot = 20L, seed = 42L)
  expect_false(is.null(fit$boot))
  expect_equal(dim(fit$boot$coefficients), c(3L, 20L))
  expect_equal(fit$boot$n.boot, 20L)
})

test_that("bagged logreg is reproducible at a fixed seed", {
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  f1 <- logreg(y ~ wt + hp, data = df, n.boot = 25L, seed = 7L)
  f2 <- logreg(y ~ wt + hp, data = df, n.boot = 25L, seed = 7L)
  expect_identical(f1$boot$coefficients, f2$boot$coefficients)
})

test_that("bagged logreg seed does not perturb the user's RNG", {
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  set.seed(123L)
  before <- .Random.seed
  logreg(y ~ wt + hp, data = df, n.boot = 10L, seed = 99L)
  expect_identical(.Random.seed, before)
})

test_that("predict on a bagged logreg returns the bag-mean probability", {
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  fit <- logreg(y ~ wt + hp, data = df, n.boot = 50L, seed = 5L)
  nd <- df[1:6, ]
  bc <- fit$boot$coefficients
  ok <- which(!apply(is.na(bc), 2L, any))
  xnew <- cbind(1, as.matrix(nd[, c("wt", "hp")]))
  eta_mat <- cbind(drop(xnew %*% fit$coefficients),
                   xnew %*% bc[, ok, drop = FALSE])
  expect_equal(unname(predict(fit, nd, type = "response")),
               unname(rowMeans(plogis(eta_mat))), tolerance = 1e-10)
})

test_that("bagged logreg predict se.fit returns the bag SD", {
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  fit <- logreg(y ~ wt + hp, data = df, n.boot = 40L, seed = 3L)
  nd <- df[1:6, ]
  p <- predict(fit, nd, type = "response", se.fit = TRUE)
  expect_length(attr(p, "se.fit"), 6L)
  expect_true(all(attr(p, "se.fit") >= 0))
})

test_that("non-bagged logreg has no boot slot", {
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  fit <- logreg(y ~ wt + hp, data = df)
  expect_null(fit$boot)
})

test_that("bagged logreg composes with weights", {
  set.seed(31)
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  w <- runif(nrow(df), 0.5, 2)
  fit <- logreg(y ~ wt + hp, data = df, weights = w,
                n.boot = 15L, seed = 8L)
  expect_equal(dim(fit$boot$coefficients), c(3L, 15L))
})
