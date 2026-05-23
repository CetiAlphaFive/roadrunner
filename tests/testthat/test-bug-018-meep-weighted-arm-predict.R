# BUG-018 regression test (v0.0.0.9055).
#
# predict.meep(nuisance = "mu0" | "mu1") returned all-NA for any
# weighted meep() fit. The arm-restricted refit inside predict.meep
# passed object$weights (length n) without slicing to the treatment-arm
# subset rows. The base learner errored length(weights) != length(y),
# the error was swallowed by tryCatch, the OOF column stayed NA, and
# .meep_apply_weights returned NA_real_ per row.
#
# Fix: subset weights by `rows` before passing to the learner refit.

test_that("predict.meep mu0 / mu1 are non-NA for weighted fits", {
  set.seed(1)
  n <- 200
  p <- 3
  X <- matrix(runif(n * p), n, p)
  D <- rbinom(n, 1, plogis(X[, 1] - 0.5))
  y <- 1 + 0.5 * X[, 1] + 0.6 * D + rnorm(n, sd = 0.3)
  w <- runif(n, 0.5, 1.5)
  fit <- meep(X, y, treatment = D, weights = w, folds = 4L,
              learners = c("ols", "logreg"), tune = "none", seed = 1)
  newx <- matrix(runif(5 * p), 5, p)
  p_mu0 <- predict(fit, newx, nuisance = "mu0")
  p_mu1 <- predict(fit, newx, nuisance = "mu1")
  expect_false(all(is.na(p_mu0)))
  expect_false(all(is.na(p_mu1)))
  expect_equal(length(p_mu0), 5L)
  expect_equal(length(p_mu1), 5L)
})
