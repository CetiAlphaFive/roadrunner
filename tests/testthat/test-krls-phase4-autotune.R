# Phase 4 tests: autotune on sigma.

# --- 1. autotune = FALSE (default) is byte-identical -----------------
test_that("Phase 4: autotune=FALSE (default) preserves back-compat", {
  set.seed(41L)
  X <- matrix(rnorm(50 * 2), 50, 2); y <- X[, 1] + rnorm(50, sd = 0.1)
  f0 <- krls(X, y, sigma = 3, lambda = 0.3,
             derivative = FALSE, vcov = TRUE)
  f1 <- krls(X, y, sigma = 3, lambda = 0.3,
             derivative = FALSE, vcov = TRUE,
             autotune = FALSE)
  expect_equal(as.numeric(f0$coeffs), as.numeric(f1$coeffs),
               tolerance = 1e-14)
  expect_null(f1$autotune)
})

# --- 2. autotune returns a valid sigma + diagnostics ----------------
test_that("Phase 4: autotune returns a sigma and diagnostics", {
  set.seed(42L)
  X <- matrix(rnorm(80 * 3), 80, 3); y <- sin(X[, 1]) + rnorm(80, sd = 0.2)
  fit <- krls(X, y, autotune = TRUE, nfold = 5L, seed.cv = 1L,
              derivative = FALSE, vcov = TRUE)
  expect_true(!is.null(fit$autotune))
  expect_true(fit$sigma %in% fit$autotune$grid)
  expect_equal(length(fit$autotune$grid), length(fit$autotune$mse))
  expect_equal(fit$autotune$winner, fit$sigma)
})

# --- 3. autotune picks a sigma close to a known-good one -----------
test_that("Phase 4: autotune prefers sigma that matches DGP scale", {
  # On a DGP where the optimal kernel is wide, autotune should pick a
  # larger-than-default sigma. We seed the DGP so dim(X) = 3 but the
  # signal lives on a slow scale -> wide kernel is best.
  set.seed(43L)
  n <- 100L
  X <- matrix(rnorm(n * 3), n, 3)
  y <- 0.5 * X[, 1] + 0.3 * X[, 2] + rnorm(n, sd = 0.2)  # nearly linear
  fit <- krls(X, y, autotune = TRUE, nfold = 5L, seed.cv = 1L,
              derivative = FALSE, vcov = TRUE)
  # autotune sigma should be at least as large as default (3) since
  # linear DGPs prefer wider kernels (smoother fits).
  expect_true(fit$autotune$winner >= 3)
})

# --- 4. autotune reproducibility via seed.cv ------------------------
test_that("Phase 4: seed.cv reproducibility for autotune", {
  set.seed(44L)
  X <- matrix(rnorm(60 * 3), 60, 3); y <- sin(X[, 1]) + rnorm(60, sd = 0.2)
  f1 <- krls(X, y, autotune = TRUE, nfold = 5L, seed.cv = 9L,
             derivative = FALSE, vcov = TRUE)
  f2 <- krls(X, y, autotune = TRUE, nfold = 5L, seed.cv = 9L,
             derivative = FALSE, vcov = TRUE)
  expect_equal(f1$sigma, f2$sigma)
  expect_equal(f1$autotune$mse, f2$autotune$mse)
})

# --- 5. custom autotune.grid honoured -------------------------------
test_that("Phase 4: custom autotune.grid is used as-is", {
  set.seed(45L)
  X <- matrix(rnorm(60 * 3), 60, 3); y <- sin(X[, 1]) + rnorm(60, sd = 0.2)
  grid <- c(0.5, 1, 2, 5, 10)
  fit <- krls(X, y, autotune = TRUE, autotune.grid = grid,
              nfold = 5L, seed.cv = 1L,
              derivative = FALSE, vcov = TRUE)
  expect_equal(fit$autotune$grid, grid)
  expect_true(fit$sigma %in% grid)
})

# --- 6. invalid autotune.grid errors --------------------------------
test_that("Phase 4: negative autotune.grid errors", {
  set.seed(46L)
  X <- matrix(rnorm(40), 20, 2); y <- rnorm(20)
  expect_error(krls(X, y, autotune = TRUE,
                    autotune.grid = c(-1, 0.5)))
})

# --- 7. RNG isolation -----------------------------------------------
test_that("Phase 4: autotune does not pollute caller's RNG stream", {
  set.seed(47L)
  X <- matrix(rnorm(30 * 2), 30, 2); y <- X[, 1] + rnorm(30, sd = 0.1)
  set.seed(123L)
  invisible(krls(X, y, autotune = TRUE, nfold = 5L, seed.cv = 1L,
                 derivative = FALSE, vcov = TRUE))
  r1 <- runif(1)
  set.seed(123L)
  r2 <- runif(1)
  expect_equal(r1, r2)
})

# --- 8. autotune + lambda.method = 'cv' compose ---------------------
test_that("Phase 4: autotune works with lambda.method='cv'", {
  set.seed(48L)
  X <- matrix(rnorm(60 * 3), 60, 3); y <- sin(X[, 1]) + rnorm(60, sd = 0.3)
  fit <- krls(X, y, autotune = TRUE, lambda.method = "cv",
              nfold = 5L, seed.cv = 1L,
              derivative = FALSE, vcov = TRUE)
  expect_true(!is.null(fit$autotune))
  expect_true(!is.null(fit$cv))
})
