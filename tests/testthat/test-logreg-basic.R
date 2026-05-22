# logreg() basic parity vs stats::glm(family = binomial).

test_that("logreg matches glm on coefficients, fitted, deviance (formula)", {
  # vs ~ wt + drat + gear is a comfortably convergent binary model;
  # several am ~ * models in mtcars are perfectly separated.
  df <- data.frame(y = as.integer(mtcars$vs), mtcars[c("wt", "drat", "gear")])
  fit <- logreg(y ~ wt + drat + gear, data = df)
  g <- glm(y ~ wt + drat + gear, data = df, family = binomial())
  expect_equal(unname(coef(fit)), unname(coef(g)), tolerance = 1e-7)
  expect_equal(unname(fit$fitted.values), unname(fitted(g)),
               tolerance = 1e-8)
  expect_equal(fit$deviance, g$deviance, tolerance = 1e-8)
  expect_equal(fit$null.deviance, g$null.deviance, tolerance = 1e-8)
  expect_equal(fit$df.residual, g$df.residual)
  expect_equal(fit$df.null, g$df.null)
})

test_that("logreg matches glm on AIC", {
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  fit <- logreg(y ~ wt + hp, data = df)
  g <- glm(y ~ wt + hp, data = df, family = binomial())
  expect_equal(fit$aic, AIC(g), tolerance = 1e-7)
})

test_that("logreg classical vcov matches the Fisher information", {
  # glm()'s own vcov is computed from the QR of the second-to-last IRLS
  # iterate, so it is itself slightly off the exact (X'WX)^-1 at the
  # converged mu; logreg uses the converged-mu Fisher information, which
  # is the textbook MLE covariance. Compare to that directly.
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  fit <- logreg(y ~ wt + hp, data = df)
  g <- glm(y ~ wt + hp, data = df, family = binomial())
  mu <- g$fitted.values
  x <- model.matrix(g)
  fisher_inv <- solve(t(x) %*% (mu * (1 - mu) * x))
  expect_equal(unname(fit$vcov), unname(fisher_inv), tolerance = 1e-7)
  # And it agrees with glm's vcov to glm's own internal precision.
  expect_equal(unname(fit$vcov), unname(vcov(g)), tolerance = 1e-3)
})

test_that("logreg default (matrix) interface matches glm", {
  x <- as.matrix(mtcars[, c("wt", "hp")])
  y <- as.integer(mtcars$am)
  fit <- logreg(x, y)
  g <- glm(y ~ x, family = binomial())
  expect_equal(unname(coef(fit)), unname(coef(g)), tolerance = 1e-7)
  expect_equal(names(coef(fit)), c("(Intercept)", "wt", "hp"))
})

test_that("logreg default interface honours intercept = FALSE", {
  x <- as.matrix(mtcars[, c("wt", "hp")])
  y <- as.integer(mtcars$am)
  fit <- logreg(x, y, intercept = FALSE)
  g <- glm(y ~ x - 1, family = binomial())
  expect_equal(unname(coef(fit)), unname(coef(g)), tolerance = 1e-7)
  expect_false("(Intercept)" %in% names(coef(fit)))
})

test_that("logreg weighted matches glm with prior weights", {
  set.seed(1)
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  w <- runif(nrow(df), 0.5, 3)
  fit <- logreg(y ~ wt + hp, data = df, weights = w)
  # glm warns "non-integer #successes" for non-integer prior weights;
  # the weighted IRLS fit is unaffected.
  g <- suppressWarnings(
    glm(y ~ wt + hp, data = df, family = binomial(), weights = w))
  expect_equal(unname(coef(fit)), unname(coef(g)), tolerance = 1e-7)
  expect_equal(unname(fit$fitted.values), unname(fitted(g)),
               tolerance = 1e-8)
  expect_equal(fit$deviance, g$deviance, tolerance = 1e-7)
})

test_that("logreg accepts a logical response", {
  df <- data.frame(y = mtcars$am > 0, mtcars[c("wt", "hp")])
  fit <- logreg(y ~ wt + hp, data = df)
  g <- glm(y ~ wt + hp, data = df, family = binomial())
  expect_equal(unname(coef(fit)), unname(coef(g)), tolerance = 1e-7)
})

test_that("logreg accepts a two-level factor response", {
  df <- data.frame(y = factor(ifelse(mtcars$am > 0, "auto", "manual")),
                   mtcars[c("wt", "hp")])
  fit <- logreg(y ~ wt + hp, data = df)
  g <- glm(y ~ wt + hp, data = df, family = binomial())
  expect_equal(unname(coef(fit)), unname(coef(g)), tolerance = 1e-7)
  expect_equal(fit$y.levels, levels(df$y))
})

test_that("logreg handles factor predictors via data frame", {
  df <- data.frame(y = as.integer(mtcars$am), wt = mtcars$wt,
                   cyl = factor(mtcars$cyl))
  fit <- logreg(y ~ wt + cyl, data = df)
  g <- glm(y ~ wt + cyl, data = df, family = binomial())
  expect_equal(unname(coef(fit)), unname(coef(g)), tolerance = 1e-7)
})

test_that("logreg rejects non-binary numeric y with a clear message", {
  x <- as.matrix(mtcars[, c("wt", "hp")])
  expect_error(logreg(x, mtcars$mpg), "0/1")
  expect_error(logreg(x, rep(c(0, 1, 2), length.out = nrow(x))), "0/1")
})

test_that("logreg rejects a factor response with the wrong level count", {
  x <- as.matrix(mtcars[, c("wt", "hp")])
  y3 <- factor(rep(c("a", "b", "c"), length.out = nrow(x)))
  expect_error(logreg(x, y3), "two levels")
})

test_that("logreg rejects mismatched lengths and bad weights", {
  x <- as.matrix(mtcars[, c("wt", "hp")])
  y <- as.integer(mtcars$am)
  expect_error(logreg(x, y[-1]), "length")
  expect_error(logreg(x, y, weights = c(1, 2)), "length")
  expect_error(logreg(x, y, weights = rep(-1, nrow(x))),
               "strictly positive")
})

test_that("logreg errors when a class is absent after NA removal", {
  x <- as.matrix(mtcars[, c("wt", "hp")])
  y <- rep(0L, nrow(x))
  expect_error(logreg(x, y), "0/1|one class")
})

test_that("logreg warns on perfect separation (non-convergence)", {
  # A predictor that perfectly separates the classes drives the MLE to
  # +/-Inf; IRLS cannot converge.
  set.seed(7)
  n <- 40L
  xsep <- c(rnorm(n / 2, -3), rnorm(n / 2, 3))
  ysep <- c(rep(0L, n / 2), rep(1L, n / 2))
  expect_warning(fit <- logreg(matrix(xsep, ncol = 1L), ysep),
                 "converge|separation")
  expect_false(fit$converged)
})

test_that("logreg default interface drops rows with NA", {
  x <- as.matrix(mtcars[, c("wt", "hp")])
  y <- as.integer(mtcars$am)
  x[3, 1] <- NA
  expect_warning(fit <- logreg(x, y), "dropped")
  g <- glm(y[-3] ~ x[-3, ], family = binomial())
  expect_equal(unname(coef(fit)), unname(coef(g)), tolerance = 1e-7)
})

test_that("logreg formula path handles NA via na.action", {
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  df$wt[5] <- NA
  fit <- logreg(y ~ wt + hp, data = df)
  g <- glm(y ~ wt + hp, data = df, family = binomial())
  expect_equal(unname(coef(fit)), unname(coef(g)), tolerance = 1e-7)
})

test_that("logreg rejects offset() terms", {
  df <- data.frame(y = as.integer(mtcars$am), mtcars[c("wt", "hp")])
  df$off <- 1
  expect_error(logreg(y ~ wt + offset(off), data = df), "offset")
})
