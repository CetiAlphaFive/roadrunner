# Phase 5 tests: bagging via n.boot.

# --- 1. n.boot = 0 (default) preserves back-compat -------------------
test_that("Phase 5: n.boot = 0 (default) is byte-identical", {
  set.seed(51L)
  X <- matrix(rnorm(50 * 2), 50, 2); y <- X[, 1] + rnorm(50, sd = 0.1)
  f0 <- krls(X, y, sigma = 2, lambda = 0.3,
             derivative = FALSE, vcov = TRUE)
  f1 <- krls(X, y, sigma = 2, lambda = 0.3,
             derivative = FALSE, vcov = TRUE,
             n.boot = 0L)
  expect_equal(as.numeric(f0$coeffs), as.numeric(f1$coeffs),
               tolerance = 1e-14)
  expect_null(f0$boot)
  expect_null(f1$boot)
})

# --- 2. n.boot > 0 stores replicates ---------------------------------
test_that("Phase 5: n.boot > 0 stores replicates + idx", {
  set.seed(52L)
  X <- matrix(rnorm(40 * 2), 40, 2); y <- X[, 1] + rnorm(40, sd = 0.1)
  fit <- krls(X, y, sigma = 2, lambda = 0.3,
              derivative = FALSE, vcov = TRUE,
              n.boot = 8L, seed.cv = 1L)
  expect_equal(fit$boot$n.boot, 8L)
  expect_equal(length(fit$boot$replicates), 8L)
  expect_equal(length(fit$boot$idx), 8L)
  expect_equal(length(fit$boot$idx[[1]]), 40L)
  # Each replicate carries the inherited sigma + lambda.
  for (b in seq_along(fit$boot$replicates)) {
    expect_equal(fit$boot$replicates[[b]]$sigma, 2)
    expect_equal(fit$boot$replicates[[b]]$lambda, 0.3)
  }
})

# --- 3. predict on bagged fit averages across replicates -------------
test_that("Phase 5: predict on bagged fit returns bag mean", {
  set.seed(53L)
  X <- matrix(rnorm(40 * 2), 40, 2); y <- X[, 1] + rnorm(40, sd = 0.1)
  f0 <- krls(X, y, sigma = 2, lambda = 0.3,
             derivative = FALSE, vcov = TRUE)
  fb <- krls(X, y, sigma = 2, lambda = 0.3,
             derivative = FALSE, vcov = TRUE,
             n.boot = 12L, seed.cv = 1L)
  p0 <- as.numeric(predict(f0, X[1:5, ])$fit)
  pb <- as.numeric(predict(fb, X[1:5, ])$fit)
  # Central pred and bag pred should differ but be in the same ballpark.
  expect_false(isTRUE(all.equal(p0, pb)))
  expect_true(all(abs(p0 - pb) < 1))
})

# --- 4. bag_sd is returned ------------------------------------------
test_that("Phase 5: predict returns bag_sd attribute", {
  set.seed(54L)
  X <- matrix(rnorm(40 * 2), 40, 2); y <- X[, 1] + rnorm(40, sd = 0.1)
  fb <- krls(X, y, sigma = 2, lambda = 0.3,
             derivative = FALSE, vcov = TRUE,
             n.boot = 8L, seed.cv = 1L)
  pr <- predict(fb, X[1:5, ])
  expect_true(!is.null(pr$bag_sd))
  expect_equal(length(pr$bag_sd), 5L)
  expect_true(all(pr$bag_sd >= 0))
})

# --- 5. reproducibility via seed.cv ---------------------------------
test_that("Phase 5: seed.cv reproducibility for bagging", {
  set.seed(55L)
  X <- matrix(rnorm(40 * 2), 40, 2); y <- X[, 1] + rnorm(40, sd = 0.1)
  f1 <- krls(X, y, sigma = 2, lambda = 0.3,
             derivative = FALSE, vcov = TRUE,
             n.boot = 6L, seed.cv = 11L)
  f2 <- krls(X, y, sigma = 2, lambda = 0.3,
             derivative = FALSE, vcov = TRUE,
             n.boot = 6L, seed.cv = 11L)
  for (b in seq_along(f1$boot$replicates)) {
    expect_equal(f1$boot$replicates[[b]]$coeffs,
                 f2$boot$replicates[[b]]$coeffs)
    expect_equal(f1$boot$idx[[b]], f2$boot$idx[[b]])
  }
})

# --- 6. predict shape stays correct on bagged fit -------------------
test_that("Phase 5: predict shape matches non-bagged", {
  set.seed(56L)
  X <- matrix(rnorm(50 * 2), 50, 2); y <- X[, 1] + rnorm(50, sd = 0.1)
  Xn <- matrix(rnorm(10 * 2), 10, 2)
  fb <- krls(X, y, sigma = 2, lambda = 0.3,
             derivative = FALSE, vcov = TRUE,
             n.boot = 5L, seed.cv = 1L)
  pr <- predict(fb, Xn)
  expect_equal(nrow(pr$fit), 10L)
})
