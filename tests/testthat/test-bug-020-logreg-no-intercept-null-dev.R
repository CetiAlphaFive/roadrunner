# BUG-020 regression test (v0.0.0.9055).
#
# logreg() reported a `null.deviance` that disagreed with stats::glm()
# for no-intercept formulas (y ~ 0 + x). The C++ engine always used
# mu0 = weighted mean of y for the null-deviance baseline; the canonical
# null model under the link is eta = 0 -> mu = 0.5 when there is no
# intercept. The R-side df.null was already adjusted, so the result was
# internally inconsistent.
#
# Fix: pass has_intercept (detected R-side as "(Intercept)" %in%
# colnames(X) or intercept = TRUE) into logreg_fit_cpp. When false, set
# mu0 = 0.5 in the null-deviance accumulator.

test_that("logreg null.deviance matches glm() for no-intercept formula", {
  df <- data.frame(y = as.integer(mtcars$vs),
                   wt = mtcars$wt,
                   drat = mtcars$drat)
  fit <- logreg(y ~ 0 + wt + drat, data = df)
  g   <- stats::glm(y ~ 0 + wt + drat, data = df, family = stats::binomial())
  expect_equal(fit$null.deviance, g$null.deviance, tolerance = 1e-6)
})

test_that("logreg null.deviance still matches glm() with an intercept", {
  df <- data.frame(y = as.integer(mtcars$vs),
                   wt = mtcars$wt,
                   drat = mtcars$drat)
  fit <- logreg(y ~ wt + drat, data = df)
  g   <- stats::glm(y ~ wt + drat, data = df, family = stats::binomial())
  expect_equal(fit$null.deviance, g$null.deviance, tolerance = 1e-6)
})
