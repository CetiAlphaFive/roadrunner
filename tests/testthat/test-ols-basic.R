# ols() basic parity vs stats::lm / lm.wfit.

test_that("ols matches lm on coefficients, fitted, residuals (formula)", {
  fit <- ols(mpg ~ wt + hp + disp, data = mtcars)
  lm0 <- lm(mpg ~ wt + hp + disp, data = mtcars)
  expect_equal(unname(coef(fit)), unname(coef(lm0)))
  expect_equal(unname(fit$fitted.values), unname(fitted(lm0)))
  expect_equal(unname(fit$residuals), unname(residuals(lm0)))
  expect_equal(fit$df.residual, lm0$df.residual)
})

test_that("ols matches lm on sigma and classical vcov", {
  fit <- ols(mpg ~ wt + hp, data = mtcars)
  lm0 <- lm(mpg ~ wt + hp, data = mtcars)
  expect_equal(fit$sigma, summary(lm0)$sigma)
  expect_equal(unname(fit$vcov), unname(vcov(lm0)))
  expect_equal(fit$rank, lm0$rank)
})

test_that("ols default (matrix) interface matches lm", {
  x <- as.matrix(mtcars[, c("wt", "hp", "disp")])
  y <- mtcars$mpg
  fit <- ols(x, y)
  lm0 <- lm(y ~ x)
  expect_equal(unname(coef(fit)), unname(coef(lm0)))
  expect_equal(names(coef(fit)), c("(Intercept)", "wt", "hp", "disp"))
})

test_that("ols default interface honours intercept = FALSE", {
  x <- as.matrix(mtcars[, c("wt", "hp")])
  y <- mtcars$mpg
  fit <- ols(x, y, intercept = FALSE)
  lm0 <- lm(y ~ x - 1)
  expect_equal(unname(coef(fit)), unname(coef(lm0)))
  expect_false("(Intercept)" %in% names(coef(fit)))
})

test_that("ols weighted matches lm.wfit / weighted lm", {
  set.seed(1)
  w <- runif(nrow(mtcars), 0.5, 3)
  fit <- ols(mpg ~ wt + hp, data = mtcars, weights = w)
  lm0 <- lm(mpg ~ wt + hp, data = mtcars, weights = w)
  expect_equal(unname(coef(fit)), unname(coef(lm0)))
  expect_equal(unname(fit$fitted.values), unname(fitted(lm0)))
  expect_equal(fit$sigma, summary(lm0)$sigma)
  expect_equal(unname(fit$vcov), unname(vcov(lm0)))
})

test_that("ols weighted (default interface) matches lm.wfit", {
  set.seed(2)
  x <- as.matrix(mtcars[, c("wt", "hp")])
  y <- mtcars$mpg
  w <- runif(length(y), 0.5, 3)
  fit <- ols(x, y, weights = w)
  wf <- lm.wfit(cbind(1, x), y, w)
  expect_equal(unname(coef(fit)), unname(wf$coefficients))
})

test_that("ols handles factor predictors via data frame", {
  df <- mtcars
  df$cyl <- factor(df$cyl)
  fit <- ols(mpg ~ wt + cyl, data = df)
  lm0 <- lm(mpg ~ wt + cyl, data = df)
  expect_equal(unname(coef(fit)), unname(coef(lm0)))
})

test_that("ols default interface expands a data frame with a factor", {
  df <- mtcars[, c("wt", "hp")]
  df$cyl <- factor(mtcars$cyl)
  fit <- ols(df, mtcars$mpg)
  lm0 <- lm(mtcars$mpg ~ wt + hp + cyl, data = df)
  expect_equal(unname(coef(fit)), unname(coef(lm0)))
})

test_that("ols errors clearly on a rank-deficient design", {
  x <- as.matrix(mtcars[, c("wt", "hp")])
  x <- cbind(x, dup = x[, "wt"])          # exact collinearity.
  expect_error(ols(x, mtcars$mpg), "rank-deficient")
})

test_that("ols rejects mismatched lengths and bad weights", {
  x <- as.matrix(mtcars[, c("wt", "hp")])
  expect_error(ols(x, mtcars$mpg[-1]), "length")
  expect_error(ols(x, mtcars$mpg, weights = c(1, 2)), "length")
  expect_error(ols(x, mtcars$mpg, weights = rep(-1, nrow(x))),
               "strictly positive")
})

test_that("ols default interface drops rows with NA", {
  x <- as.matrix(mtcars[, c("wt", "hp")])
  y <- mtcars$mpg
  x[3, 1] <- NA
  expect_warning(fit <- ols(x, y), "dropped")
  lm0 <- lm(y[-3] ~ x[-3, ])
  expect_equal(unname(coef(fit)), unname(coef(lm0)))
})

test_that("ols formula path handles NA via na.action", {
  df <- mtcars
  df$wt[5] <- NA
  fit <- ols(mpg ~ wt + hp, data = df)
  lm0 <- lm(mpg ~ wt + hp, data = df)
  expect_equal(unname(coef(fit)), unname(coef(lm0)))
})
