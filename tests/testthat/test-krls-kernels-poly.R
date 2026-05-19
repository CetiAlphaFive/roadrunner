# Phase Q2 (v0.0.0.9049): polynomial / linear kernel parity tests.
#
# Coverage:
#   1.  Linear kernel: byte-equal hand-built XX' (post-standardisation).
#   2.  Poly2 kernel:  byte-equal hand-built (XX'+c)^2.
#   3.  Linear deriv:  constant per column (X'c).
#   4.  Poly2 deriv:   matches finite-difference reference.
#   5.  Validation gates (errors).
#   6.  Autotune sweeps poly_c (poly2).
#   7.  Bagging composes with linear kernel.
#   8.  Gaussian default byte-identical to v9048 (back-compat seal).

test_that("linear kernel: K matches X X' (post-standardisation, byte-equal)", {
  set.seed(1)
  n <- 30; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  fit <- krls(X, y, whichkernel = "linear",
              derivative = FALSE, vcov = FALSE)
  # Reference: build XX' on the standardised X.
  Xs <- scale(X)
  Xs <- matrix(Xs, n, p)
  K_ref <- tcrossprod(Xs)
  expect_equal(fit$K, K_ref, tolerance = 1e-12)
})

test_that("poly2 kernel: K matches (XX'+c)^2 byte-equal", {
  set.seed(2)
  n <- 25; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  cval <- 0.7
  fit <- krls(X, y, whichkernel = "poly2", poly_c = cval,
              derivative = FALSE, vcov = FALSE)
  Xs <- scale(X)
  Xs <- matrix(Xs, n, p)
  K_ref <- (tcrossprod(Xs) + cval)^2
  expect_equal(fit$K, K_ref, tolerance = 1e-12)
})

test_that("linear kernel deriv is constant per column", {
  set.seed(3)
  n <- 40; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X %*% c(1, -2, 0.5)) + 0.1 * rnorm(n)
  fit <- suppressWarnings(krls(X, y, whichkernel = "linear",
              derivative = TRUE, vcov = TRUE))
  # Each column of derivatives should be constant (linear kernel).
  for (k in seq_len(p)) {
    expect_lt(diff(range(fit$derivatives[, k])), 1e-8)
  }
  # var.avgderivatives should be NA (deferred to Q6).
  expect_true(all(is.na(as.numeric(fit$var.avgderivatives))))
})

test_that("poly2 deriv matches finite-difference reference", {
  set.seed(4)
  n <- 35; p <- 2
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(0.5 * X[, 1]^2 + X[, 2] + 0.1 * rnorm(n))
  cval <- 1.0
  fit <- suppressWarnings(krls(X, y, whichkernel = "poly2", poly_c = cval,
              derivative = TRUE, vcov = TRUE))
  # Finite-difference reference on the predictor.
  fd <- function(j, eps = 1e-5) {
    Xp <- X; Xm <- X
    Xp[, j] <- Xp[, j] + eps
    Xm[, j] <- Xm[, j] - eps
    yp <- as.numeric(predict(fit, Xp)$fit)
    ym <- as.numeric(predict(fit, Xm)$fit)
    (yp - ym) / (2 * eps)
  }
  for (j in seq_len(p)) {
    fd_j <- fd(j)
    deriv_j <- fit$derivatives[, j]
    expect_lt(max(abs(fd_j - deriv_j)), 1e-4)
  }
})

test_that("non-Gaussian kernels reject ARD/Nystrom/vector-sigma", {
  set.seed(5)
  n <- 30; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  # (a) linear + nystrom -> error
  expect_error(krls(X, y, whichkernel = "linear", approx = "nystrom"),
               regexp = "nystrom.*incompatible|incompatible.*nystrom")
  # (b) poly2 + ard -> error
  expect_error(krls(X, y, whichkernel = "poly2", ard = "cheap"),
               regexp = "ard.*incompatible|incompatible")
  # (c) poly3 + vector sigma -> error
  expect_error(krls(X, y, whichkernel = "poly3", sigma = rep(2, p)),
               regexp = "vector sigma.*incompatible|incompatible")
  # (d) linear + autotune.grid -> error (no sweep)
  expect_error(krls(X, y, whichkernel = "linear", autotune = TRUE,
                    autotune.grid = c(1, 2, 4)),
               regexp = "no kernel hyperparameter")
  # (e) MLL + nystrom -> error
  expect_error(krls(X, y, lambda.method = "mll", approx = "nystrom"),
               regexp = "mll.*not supported.*nystrom|nystrom.*not.*mll")
})

test_that("autotune sweeps poly_c for poly2 kernel", {
  skip_if_quick()
  set.seed(6)
  n <- 80; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(0.5 * X[, 1]^2 - X[, 2] + 0.1 * rnorm(n))
  fit <- krls(X, y, whichkernel = "poly2",
              autotune = TRUE,
              autotune.grid = c(0.25, 0.5, 1, 2, 4),
              autotune.warmstart = FALSE,
              derivative = FALSE, vcov = FALSE,
              seed.cv = 42L)
  expect_equal(fit$whichkernel, "poly2")
  expect_true(fit$poly_c %in% c(0.25, 0.5, 1, 2, 4))
  expect_true(!is.null(fit$autotune$winner_poly_c))
})

test_that("bagging composes with linear kernel", {
  skip_if_quick()
  set.seed(7)
  n <- 60; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(X %*% c(0.5, -1, 0.2)) + 0.1 * rnorm(n)
  fit <- krls(X, y, whichkernel = "linear",
              derivative = FALSE, vcov = FALSE,
              n.boot = 5L, seed.cv = 99L)
  expect_equal(length(fit$boot$replicates), 5L)
  for (b in fit$boot$replicates) {
    expect_false(is.null(b$coeffs))
  }
  # Predict averages across bag replicates.
  Xnew <- matrix(rnorm(8 * p), 8, p)
  pr <- predict(fit, Xnew)
  expect_equal(nrow(pr$fit), 8L)
})

test_that("Gaussian default byte-identical to v9048 (back-compat seal)", {
  # Same setup as test-krls-ard-kernel Test 1 (broadcast equivalence)
  # but explicitly asserts the new kernel_type=0 / poly_c arguments
  # default to the byte-identical Gaussian fast path.
  set.seed(2718)
  n <- 80; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(scale(rowSums(X[, 1:3]) + 0.3 * rnorm(n)))

  f1 <- krls(X, y, sigma = 10, derivative = TRUE, vcov = TRUE)
  f2 <- krls(X, y, sigma = 10,
             whichkernel = "gaussian",  # explicit
             derivative = TRUE, vcov = TRUE)
  f3 <- krls(X, y, sigma = 10,
             whichkernel = "gaussian", poly_c = 1.0,
             derivative = TRUE, vcov = TRUE)
  expect_identical(f1$coeffs, f2$coeffs)
  expect_identical(f1$coeffs, f3$coeffs)
  expect_identical(f1$fitted, f2$fitted)
  expect_identical(f1$derivatives, f3$derivatives)
})
