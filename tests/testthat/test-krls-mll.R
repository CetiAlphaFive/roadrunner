# Phase Q2 (v0.0.0.9049): Type-II marginal log-likelihood lambda selection.

test_that(".krls_mll_loss matches hand-computed formula", {
  set.seed(21)
  n <- 30
  dvals <- sort(runif(n, 0.01, 5), decreasing = TRUE)
  Vty <- rnorm(n)
  lam <- 0.5

  # roadrunner's helper.
  loss_r <- roadrunner:::.krls_mll_loss(dvals, Vty, lam, n)
  # Hand-computed reference (omit constants).
  ref <- sum(Vty^2 / (dvals + lam)) + sum(log(dvals + lam))
  expect_lt(abs(loss_r - ref), 1e-12)
})

test_that("MLL lambda agrees with hand-computed 100-pt grid", {
  set.seed(22)
  n <- 50; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(scale(sin(X[, 1]) + 0.5 * X[, 2] + 0.2 * rnorm(n)))
  fit <- krls(X, y, lambda.method = "mll",
              derivative = FALSE, vcov = FALSE)

  # Build the same eigen-basis hand-side.
  Xs <- matrix(scale(X), n, p)
  ys <- as.numeric(scale(y))
  sigma_vec <- fit$sigma_vec
  K <- krls_kernel_cpp(Xs, sigma_vec)
  eo <- krls_eig_cpp(K)
  dvals <- eo$values
  V <- eo$vectors
  Vty <- as.numeric(crossprod(V, ys))

  # 100-pt log-spaced grid spanning a wide bracket around fit$lambda.
  lam_lo <- max(.Machine$double.eps, fit$lambda * 1e-4)
  lam_hi <- fit$lambda * 1e4
  grid <- exp(seq(log(lam_lo), log(lam_hi), length.out = 100L))
  losses <- vapply(grid,
                   function(l) roadrunner:::.krls_mll_loss(dvals, Vty, l, n),
                   numeric(1L))
  lam_grid <- grid[which.min(losses)]
  # The fit lambda found via golden-section should be within a factor
  # of the nearest grid point (the grid is only 100 cells wide on log
  # scale; allow factor 2 wiggle).
  expect_lt(abs(log10(fit$lambda) - log10(lam_grid)), 0.5)
})

test_that("MLL and LOO usually pick different lambdas on real data", {
  set.seed(23)
  n <- 80; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(scale(sin(X[, 1]) + 0.5 * X[, 2]^2 + 0.3 * rnorm(n)))
  fit_loo <- krls(X, y, lambda.method = "loo",
                  derivative = FALSE, vcov = FALSE)
  fit_mll <- krls(X, y, lambda.method = "mll",
                  derivative = FALSE, vcov = FALSE)
  # Allow them to be different (they should be, in general).
  # Don't enforce a sign — just that they aren't bit-identical.
  expect_false(identical(fit_loo$lambda, fit_mll$lambda))
})

test_that("MLL + Nystrom errors with clear message", {
  set.seed(24)
  n <- 50; p <- 2
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  expect_error(
    krls(X, y, lambda.method = "mll", approx = "nystrom",
         derivative = FALSE, vcov = FALSE),
    regexp = "mll.*not supported.*nystrom"
  )
})

test_that("MLL composes with autotune (Gaussian, scalar sigma)", {
  skip_if_quick()
  set.seed(25)
  n <- 70; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(scale(sin(X[, 1]) + 0.3 * rnorm(n)))
  fit <- krls(X, y, lambda.method = "mll",
              autotune = TRUE,
              autotune.grid = c(2, 5, 10, 20, 50),
              autotune.warmstart = FALSE,
              derivative = FALSE, vcov = FALSE,
              seed.cv = 1L)
  expect_true(!is.null(fit$autotune))
  expect_true(fit$sigma %in% c(2, 5, 10, 20, 50))
})

test_that("MLL composes with whichkernel = 'poly2'", {
  set.seed(26)
  n <- 50; p <- 2
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(scale(0.3 * X[, 1]^2 + X[, 2] + 0.2 * rnorm(n)))
  fit <- krls(X, y, whichkernel = "poly2",
              lambda.method = "mll",
              derivative = FALSE, vcov = FALSE)
  expect_equal(fit$whichkernel, "poly2")
  expect_true(fit$lambda > 0)
})
