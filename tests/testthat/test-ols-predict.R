# ols() predict: intervals vs predict.lm(interval = ).

test_that("predict.ols point predictions match predict.lm", {
  fit <- ols(mpg ~ wt + hp, data = mtcars)
  lm0 <- lm(mpg ~ wt + hp, data = mtcars)
  nd <- mtcars[1:8, ]
  expect_equal(unname(predict(fit, nd)),
               unname(predict(lm0, nd)))
})

test_that("predict.ols with newdata = NULL returns fitted values", {
  fit <- ols(mpg ~ wt + hp, data = mtcars)
  expect_equal(predict(fit), fit$fitted.values)
})

test_that("predict.ols confidence intervals match predict.lm", {
  fit <- ols(mpg ~ wt + hp, data = mtcars)
  lm0 <- lm(mpg ~ wt + hp, data = mtcars)
  nd <- mtcars[1:10, ]
  ci_rr <- predict(fit, nd, interval = "confidence")
  ci_lm <- predict(lm0, nd, interval = "confidence")
  expect_equal(unname(ci_rr), unname(ci_lm))
})

test_that("predict.ols prediction intervals match predict.lm", {
  fit <- ols(mpg ~ wt + hp, data = mtcars)
  lm0 <- lm(mpg ~ wt + hp, data = mtcars)
  nd <- mtcars[1:10, ]
  pi_rr <- predict(fit, nd, interval = "prediction")
  pi_lm <- predict(lm0, nd, interval = "prediction")
  expect_equal(unname(pi_rr), unname(pi_lm))
})

test_that("predict.ols se.fit matches predict.lm se.fit", {
  fit <- ols(mpg ~ wt + hp, data = mtcars)
  lm0 <- lm(mpg ~ wt + hp, data = mtcars)
  nd <- mtcars[1:10, ]
  p_rr <- predict(fit, nd, se.fit = TRUE)
  p_lm <- predict(lm0, nd, se.fit = TRUE)
  expect_equal(unname(attr(p_rr, "se.fit")), unname(p_lm$se.fit))
})

test_that("predict.ols robust intervals widen with HC3", {
  skip_if_not_installed("sandwich")
  fit <- ols(mpg ~ wt + hp, data = mtcars)
  nd <- mtcars[1:10, ]
  ci_cl <- predict(fit, nd, interval = "confidence")
  ci_rb <- predict(fit, nd, interval = "confidence", robust = "HC3")
  # The HC3 SEs differ; widths must change (not identical).
  expect_false(isTRUE(all.equal(ci_cl[, "upr"] - ci_cl[, "lwr"],
                                ci_rb[, "upr"] - ci_rb[, "lwr"])))
})

test_that("predict.ols matrix newdata works", {
  x <- as.matrix(mtcars[, c("wt", "hp")])
  fit <- ols(x, mtcars$mpg)
  p <- predict(fit, x[1:5, ])
  expect_length(p, 5L)
  expect_equal(unname(p), unname(fit$fitted.values[1:5]))
})

test_that("predict.ols rejects newdata with wrong column count", {
  x <- as.matrix(mtcars[, c("wt", "hp")])
  fit <- ols(x, mtcars$mpg)
  expect_error(predict(fit, x[, 1, drop = FALSE]), "columns")
})

test_that("predict.ols prediction interval uses varmod = 'const'", {
  fit_none <- ols(mpg ~ wt + hp, data = mtcars)
  fit_vm <- ols(mpg ~ wt + hp, data = mtcars, varmod = "const")
  nd <- mtcars[1:6, ]
  pi_none <- predict(fit_none, nd, interval = "prediction")
  pi_vm <- predict(fit_vm, nd, interval = "prediction")
  # varmod = "const" stores a (df-adjusted) residual SD; the PI is a
  # valid finite interval bracketing the point estimate.
  expect_true(all(pi_vm[, "lwr"] < pi_vm[, "fit"]))
  expect_true(all(pi_vm[, "upr"] > pi_vm[, "fit"]))
  expect_equal(dim(pi_vm), dim(pi_none))
})
