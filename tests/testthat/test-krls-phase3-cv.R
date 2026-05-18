# Phase 3 tests: CV-based lambda selection.

# --- 1. lambda.method default = 'loo' is byte-identical to pre-phase --
test_that("Phase 3: lambda.method='loo' (default) preserves back-compat", {
  set.seed(31L)
  X <- matrix(rnorm(50 * 2), 50, 2); y <- X[, 1] + rnorm(50, sd = 0.1)
  f0 <- krls(X, y, sigma = 2, derivative = FALSE, vcov = TRUE)
  f1 <- krls(X, y, sigma = 2, derivative = FALSE, vcov = TRUE,
             lambda.method = "loo")
  expect_equal(f0$lambda, f1$lambda)
  expect_equal(as.numeric(f0$coeffs), as.numeric(f1$coeffs),
               tolerance = 1e-14)
})

# --- 2. lambda.method='cv' returns a usable lambda and diagnostics ----
test_that("Phase 3: lambda.method='cv' returns lambda + diagnostics", {
  set.seed(32L)
  X <- matrix(rnorm(60 * 3), 60, 3); y <- sin(X[, 1]) + rnorm(60, sd = 0.2)
  fit <- krls(X, y, sigma = 3, derivative = FALSE, vcov = TRUE,
              lambda.method = "cv", nfold = 5L, seed.cv = 1L)
  expect_true(is.numeric(fit$lambda) && fit$lambda > 0)
  expect_true(!is.null(fit$cv))
  expect_equal(fit$cv$method, "cv")
  expect_equal(fit$cv$nfold, 5L)
  expect_equal(length(fit$cv$lambda.grid), 30L)
  expect_equal(length(fit$cv$mean_mse), 30L)
  expect_true(is.matrix(fit$cv$mse_per_fold))
})

# --- 3. seed.cv reproducibility ---------------------------------------
test_that("Phase 3: seed.cv produces reproducible lambda", {
  set.seed(33L)
  X <- matrix(rnorm(60 * 3), 60, 3); y <- sin(X[, 1]) + rnorm(60, sd = 0.2)
  f1 <- krls(X, y, sigma = 3, derivative = FALSE, vcov = TRUE,
             lambda.method = "cv", nfold = 5L, seed.cv = 42L)
  f2 <- krls(X, y, sigma = 3, derivative = FALSE, vcov = TRUE,
             lambda.method = "cv", nfold = 5L, seed.cv = 42L)
  expect_equal(f1$lambda, f2$lambda)
  expect_equal(f1$cv$mean_mse, f2$cv$mean_mse)
})

# --- 4. nfold = n approximates LOO ------------------------------------
test_that("Phase 3: nfold = n yields lambda close to LOO closed-form", {
  set.seed(34L)
  n <- 40L
  X <- matrix(rnorm(n * 2), n, 2); y <- X[, 1] + rnorm(n, sd = 0.1)
  # Tighten the grid by passing a custom range around the LOO answer.
  f_loo <- krls(X, y, sigma = 2, derivative = FALSE, vcov = TRUE)
  fine_grid <- exp(seq(log(f_loo$lambda * 0.1), log(f_loo$lambda * 10),
                       length.out = 80L))
  f_loocv <- krls(X, y, sigma = 2, derivative = FALSE, vcov = TRUE,
                  lambda.method = "cv", nfold = n, seed.cv = 1L,
                  lambda.grid = fine_grid)
  # Closed-form LOO and per-fold n-fold LOO are both LOO estimators
  # but compute on different (full vs leave-one) kernel decompositions
  # and use a discrete grid. They are typically within an order of
  # magnitude; allow a generous factor of 30 to keep this test robust.
  expect_true(abs(log(f_loocv$lambda) - log(f_loo$lambda)) < log(30))
})

# --- 5. 1-SE rule picks a larger lambda -------------------------------
test_that("Phase 3: cv.1se picks lambda >= best", {
  set.seed(35L)
  X <- matrix(rnorm(80 * 3), 80, 3); y <- sin(X[, 1]) + rnorm(80, sd = 0.3)
  f_min <- krls(X, y, sigma = 3, derivative = FALSE, vcov = TRUE,
                lambda.method = "cv", nfold = 5L, seed.cv = 7L,
                cv.1se = FALSE)
  f_1se <- krls(X, y, sigma = 3, derivative = FALSE, vcov = TRUE,
                lambda.method = "cv", nfold = 5L, seed.cv = 7L,
                cv.1se = TRUE)
  # Same grid + seed; 1-SE is at least as conservative (i.e., >= argmin).
  expect_true(f_1se$lambda >= f_min$lambda)
})

# --- 6. custom lambda.grid honoured ----------------------------------
test_that("Phase 3: custom lambda.grid is used as-is", {
  set.seed(36L)
  X <- matrix(rnorm(40 * 2), 40, 2); y <- X[, 1] + rnorm(40, sd = 0.1)
  grid <- c(0.01, 0.1, 1, 10)
  fit <- krls(X, y, sigma = 2, derivative = FALSE, vcov = TRUE,
              lambda.method = "cv", nfold = 5L, seed.cv = 9L,
              lambda.grid = grid)
  expect_equal(fit$cv$lambda.grid, grid)
  expect_true(fit$lambda %in% grid)
})

# --- 7. invalid lambda.grid errors ------------------------------------
test_that("Phase 3: negative or non-numeric lambda.grid errors", {
  set.seed(37L)
  X <- matrix(rnorm(20), 10, 2); y <- rnorm(10)
  expect_error(krls(X, y, sigma = 2, lambda.method = "cv",
                    lambda.grid = c(-1, 0.1)))
})

# --- 8. nfold > n errors ----------------------------------------------
test_that("Phase 3: nfold > n is rejected", {
  set.seed(38L)
  X <- matrix(rnorm(20), 10, 2); y <- rnorm(10)
  expect_error(krls(X, y, sigma = 2, lambda.method = "cv", nfold = 100L),
               "exceeds nrow")
})

# --- 9. RNG isolation: CV doesn't change user's RNG stream ------------
test_that("Phase 3: seed.cv does not pollute caller's RNG stream", {
  set.seed(39L)
  X <- matrix(rnorm(30 * 2), 30, 2); y <- X[, 1] + rnorm(30, sd = 0.1)
  set.seed(123L)
  invisible(krls(X, y, sigma = 2, derivative = FALSE, vcov = TRUE,
                  lambda.method = "cv", nfold = 5L, seed.cv = 1L))
  r1 <- runif(1)
  set.seed(123L)
  r2 <- runif(1)
  expect_equal(r1, r2)
})

# --- 10. CV-cv ~ matches lambda quality on a clean DGP ----------------
test_that("Phase 3: CV lambda yields reasonable holdout MSE", {
  set.seed(40L)
  n <- 100L
  Xtr <- matrix(rnorm(n * 3), n, 3); ytr <- sin(Xtr[, 1]) + rnorm(n, sd = 0.2)
  Xte <- matrix(rnorm(50 * 3), 50, 3); yte <- sin(Xte[, 1]) + rnorm(50, sd = 0.2)
  f_loo <- krls(Xtr, ytr, sigma = 3, derivative = FALSE, vcov = TRUE)
  f_cv  <- krls(Xtr, ytr, sigma = 3, derivative = FALSE, vcov = TRUE,
                lambda.method = "cv", nfold = 5L, seed.cv = 1L)
  p_loo <- predict(f_loo, newdata = Xte)$fit
  p_cv  <- predict(f_cv,  newdata = Xte)$fit
  mse_loo <- mean((yte - p_loo)^2)
  mse_cv  <- mean((yte - p_cv)^2)
  # CV shouldn't be dramatically worse than LOO on a clean DGP.
  expect_true(mse_cv < 2 * mse_loo)
})
