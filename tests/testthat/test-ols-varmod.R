# ols() residual-variance model (shared ares varmod helper).

test_that("varmod = 'none' stores no variance model", {
  fit <- ols(mpg ~ wt + hp, data = mtcars)
  expect_null(fit$varmod)
})

test_that("varmod = 'const' stores a single residual SD", {
  fit <- ols(mpg ~ wt + hp, data = mtcars, varmod = "const")
  expect_false(is.null(fit$varmod))
  expect_identical(fit$varmod$type, "const")
  expect_true(is.finite(fit$varmod$sigma_hat))
  expect_true(fit$varmod$sigma_hat > 0)
})

test_that("varmod = 'lm' stores intercept + slope", {
  fit <- ols(mpg ~ wt + hp, data = mtcars, varmod = "lm")
  expect_false(is.null(fit$varmod))
  expect_identical(fit$varmod$type, "lm")
  expect_true(is.finite(fit$varmod$intercept))
  expect_true(is.finite(fit$varmod$slope))
})

test_that("varmod = 'const' sigma_hat matches the residual SE", {
  fit <- ols(mpg ~ wt + hp, data = mtcars, varmod = "const")
  # The varmod helper uses df = max(n - p, 1), the same df as the fit;
  # so sigma_hat coincides with the fitted residual standard error.
  expect_equal(fit$varmod$sigma_hat, fit$sigma)
})

test_that("prediction interval with varmod = 'lm' is finite and ordered", {
  fit <- ols(mpg ~ wt + hp, data = mtcars, varmod = "lm")
  nd <- mtcars[1:8, ]
  pi <- predict(fit, nd, interval = "prediction")
  expect_true(all(is.finite(pi)))
  expect_true(all(pi[, "lwr"] < pi[, "fit"]))
  expect_true(all(pi[, "upr"] > pi[, "fit"]))
})

test_that("varmod composes with weighted ols", {
  set.seed(21)
  w <- runif(nrow(mtcars), 0.5, 2)
  fit <- ols(mpg ~ wt + hp, data = mtcars, weights = w, varmod = "const")
  expect_identical(fit$varmod$type, "const")
  expect_true(fit$varmod$sigma_hat > 0)
})
