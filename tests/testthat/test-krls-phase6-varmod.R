# Phase 6 tests: varmod = 'const' + interval = 'pint'.

# --- 1. varmod = 'none' (default) preserves back-compat -------------
test_that("Phase 6: varmod='none' (default) is byte-identical", {
  set.seed(61L)
  X <- matrix(rnorm(50 * 2), 50, 2); y <- X[, 1] + rnorm(50, sd = 0.1)
  f0 <- krls(X, y, sigma = 2, lambda = 0.3,
             derivative = FALSE, vcov = TRUE)
  f1 <- krls(X, y, sigma = 2, lambda = 0.3,
             derivative = FALSE, vcov = TRUE,
             varmod = "none")
  expect_equal(as.numeric(f0$coeffs), as.numeric(f1$coeffs),
               tolerance = 1e-14)
  expect_null(f0$varmod)
  expect_null(f1$varmod)
})

# --- 2. varmod = 'const' stashes sigma_hat + effdf ----------------
test_that("Phase 6: varmod='const' stores sigma_hat + effdf", {
  set.seed(62L)
  X <- matrix(rnorm(80 * 2), 80, 2); y <- X[, 1] + rnorm(80, sd = 0.3)
  fit <- krls(X, y, sigma = 2, lambda = 0.3,
              derivative = FALSE, vcov = TRUE,
              varmod = "const")
  expect_true(!is.null(fit$varmod))
  expect_equal(fit$varmod$type, "const")
  expect_true(fit$varmod$sigma_hat > 0)
  expect_true(fit$varmod$effdf > 0)
  expect_true(fit$varmod$effdf < nrow(X))
})

# --- 3. interval='pint' errors without varmod ----------------------
test_that("Phase 6: interval='pint' without varmod errors", {
  set.seed(63L)
  X <- matrix(rnorm(40 * 2), 40, 2); y <- X[, 1] + rnorm(40, sd = 0.1)
  fit <- krls(X, y, sigma = 2, lambda = 0.3,
              derivative = FALSE, vcov = TRUE)
  expect_error(predict(fit, X, interval = "pint"), "varmod")
})

# --- 4. interval='pint' returns (fit, lwr, upr) -------------------
test_that("Phase 6: interval='pint' returns (fit, lwr, upr)", {
  set.seed(64L)
  X <- matrix(rnorm(50 * 2), 50, 2); y <- X[, 1] + rnorm(50, sd = 0.2)
  fit <- krls(X, y, sigma = 2, lambda = 0.3,
              derivative = FALSE, vcov = TRUE,
              varmod = "const")
  pi <- predict(fit, X, interval = "pint", level = 0.95)
  expect_equal(colnames(pi), c("fit", "lwr", "upr"))
  expect_true(all(pi[, "lwr"] < pi[, "fit"]))
  expect_true(all(pi[, "fit"] < pi[, "upr"]))
})

# --- 5. empirical 95% coverage on clean DGP in [0.90, 0.99] -------
test_that("Phase 6: 95% PI coverage on clean DGP within [0.90, 0.99]", {
  set.seed(65L)
  n_tr <- 300L
  Xtr <- matrix(rnorm(n_tr * 2), n_tr, 2)
  ytr <- sin(Xtr[, 1]) + 0.3 * Xtr[, 2] + rnorm(n_tr, sd = 0.25)
  fit <- krls(Xtr, ytr, sigma = 2, derivative = FALSE, vcov = TRUE,
              varmod = "const")
  n_te <- 1000L
  Xte <- matrix(rnorm(n_te * 2), n_te, 2)
  yte <- sin(Xte[, 1]) + 0.3 * Xte[, 2] + rnorm(n_te, sd = 0.25)
  pi <- predict(fit, Xte, interval = "pint", level = 0.95)
  cov <- mean(yte >= pi[, "lwr"] & yte <= pi[, "upr"])
  expect_gte(cov, 0.90)
  expect_lte(cov, 0.99)
})

# --- 6. level controls width ---------------------------------------
test_that("Phase 6: level=0.99 gives wider intervals than level=0.95", {
  set.seed(66L)
  X <- matrix(rnorm(60 * 2), 60, 2); y <- X[, 1] + rnorm(60, sd = 0.3)
  fit <- krls(X, y, sigma = 2, lambda = 0.3,
              derivative = FALSE, vcov = TRUE,
              varmod = "const")
  p95 <- predict(fit, X, interval = "pint", level = 0.95)
  p99 <- predict(fit, X, interval = "pint", level = 0.99)
  expect_true(mean(p99[, "upr"] - p99[, "lwr"]) >
              mean(p95[, "upr"] - p95[, "lwr"]))
})

# --- 7. interval='pint' on NULL newdata returns in-sample PI -------
test_that("Phase 6: predict(fit, NULL, interval='pint') returns in-sample PI", {
  set.seed(67L)
  X <- matrix(rnorm(40 * 2), 40, 2); y <- X[, 1] + rnorm(40, sd = 0.2)
  fit <- krls(X, y, sigma = 2, lambda = 0.3,
              derivative = FALSE, vcov = TRUE,
              varmod = "const")
  pi <- predict(fit, NULL, interval = "pint", level = 0.95)
  expect_equal(nrow(pi), nrow(X))
  expect_equal(colnames(pi), c("fit", "lwr", "upr"))
})
