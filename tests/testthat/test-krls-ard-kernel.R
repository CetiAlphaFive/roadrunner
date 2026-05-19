# Phase 2a tests: per-feature Gaussian bandwidth (ARD) plumbing.
# Spec: inst/plans/0002-phase2a-ard-kernel.md
#
# Tests:
#   1. Broadcast equivalence (scalar == length-p constant vector).
#   2. Back-compat byte-equal vs v0.0.0.9044 snapshot.
#   3. Selective lengthscale recovery (kill irrelevant dims with huge sigma).
#   4. Per-feature marginal-effect scale (closed-form reference).
#   5. Determinism across threads under ARD.
#   6. Predict under ARD matches an outer-loop reference.
#   7. Validation errors.

test_that("ARD: scalar and length-p-constant sigma are byte-identical", {
  set.seed(2718)
  n <- 80; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(scale(rowSums(X[, 1:3]) + 0.3 * rnorm(n)))

  f1 <- krls(X, y, sigma = 10, derivative = TRUE, vcov = TRUE)
  f2 <- krls(X, y, sigma = rep(10, p),
             derivative = TRUE, vcov = TRUE)

  expect_identical(f1$coeffs, f2$coeffs)
  expect_identical(f1$fitted, f2$fitted)
  expect_identical(f1$derivatives, f2$derivatives)
  expect_identical(f1$avgderivatives, f2$avgderivatives)
  expect_identical(f1$var.avgderivatives, f2$var.avgderivatives)
})

test_that("ARD: scalar sigma path is byte-identical to v0.0.0.9044 snapshot", {
  snap_file <- testthat::test_path("data", "baseline_v9044.rds")
  skip_if_not(file.exists(snap_file),
              "Baseline snapshot baseline_v9044.rds missing.")

  baseline <- readRDS(snap_file)
  fit <- krls(baseline$X, baseline$y, derivative = TRUE, vcov = TRUE)

  # Back-compat: BLAS-derived fields use a tight numeric tolerance to
  # survive ULP-level differences across BLAS implementations (the
  # baseline_v9044.rds snapshot was captured on one machine; CI runners
  # use different BLAS builds). Within-machine byte-identity is still
  # verified by Tests 1 and 5. Scalars sigma+lambda derive from FP-
  # ordered code with no BLAS gemm path, so identical() remains.
  expect_equal(fit$coeffs,             baseline$coeffs,             tolerance = 1e-10)
  expect_equal(fit$fitted,             baseline$fitted,             tolerance = 1e-10)
  expect_equal(fit$R2,                 baseline$R2,                 tolerance = 1e-10)
  expect_equal(fit$Looe,               baseline$Looe,               tolerance = 1e-10)
  expect_identical(fit$sigma,              baseline$sigma)
  expect_identical(fit$lambda,             baseline$lambda)
  expect_equal(fit$avgderivatives,     baseline$avgderivatives,     tolerance = 1e-10)
  expect_equal(fit$var.avgderivatives, baseline$var.avgderivatives, tolerance = 1e-10)
})

test_that("ARD: huge sigmas on irrelevant dims recover single-feature fit", {
  set.seed(13)
  n <- 120; p <- 6
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(X[, 1])

  # Standardise X to compare with the projected fit on its raw scale.
  f_ard <- krls(X, y, sigma = c(1, rep(1e10, p - 1L)),
                derivative = FALSE, vcov = FALSE)
  f_x1  <- krls(X[, 1, drop = FALSE], y, sigma = 1,
                derivative = FALSE, vcov = FALSE)

  # When the kernel ignores x_{2..p} (huge sigma -> distance contribution
  # -> 0), the ARD fit on the full design must reduce to a univariate
  # kernel regression on x_1.
  # Equivalence is approximate: the huge sigma leaves O(d^2/sigma)
  # residual kernel contribution that propagates through lambda. Across
  # platforms (BLAS / LAPACK eigen) this floor is O(1e-2). Assert near-
  # perfect correlation plus a loose abs-diff floor.
  expect_gt(cor(f_ard$fitted, f_x1$fitted), 0.999)
  expect_lt(max(abs(f_ard$fitted - f_x1$fitted)), 0.1)
})

test_that("ARD: per-feature marginal-effect scale matches closed-form", {
  # Use a hand-coded n=4, p=2 design with known X, c, sigma_vec.
  # Reference: derivmat[i, k] = (-2/sigma_k) * (X_ik * (Kc)_i -
  #                                              (K diag(c) X)_ik).
  X <- rbind(c(0.5, -0.5),
             c(1.0,  0.5),
             c(-0.5, 0.0),
             c(0.0,  1.0))
  cvec <- c(1, 0, -1, 0)
  sigma_vec <- c(2, 8)

  # Build K in R (Gaussian, ARD).
  K <- matrix(0, 4, 4)
  for (i in 1:4) for (j in 1:4) {
    s <- sum((X[i, ] - X[j, ])^2 / sigma_vec)
    K[i, j] <- exp(-s)
  }

  # R reference: (-2/sigma_k) * (X_ik * (Kc)_i - (K diag(c) X)_ik).
  Kc <- as.numeric(K %*% cvec)
  KdiagcX <- K %*% diag(cvec) %*% X    # n x p
  Xs1 <- X * Kc                        # broadcast n-vec * each col
  D_ref <- Xs1 - KdiagcX               # n x p, unscaled
  scl <- -2.0 / sigma_vec
  derv_ref <- sweep(D_ref, 2L, scl, "*")

  # C++ call.
  derv_cpp <- krls_deriv_cpp(X, K, cvec, sigma_vec)

  expect_lt(max(abs(derv_cpp - derv_ref)), 1e-12)
})

test_that("ARD: nthreads=1 and nthreads=4 produce byte-equal coeffs", {
  set.seed(5)
  n <- 200; p <- 8
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(scale(rowSums(X[, 1:3]^2) + 0.3 * rnorm(n)))
  sv <- exp(seq(log(2), log(20), length.out = p))   # heterogeneous ARD

  f1 <- krls(X, y, sigma = sv, nthreads = 1L,
             derivative = FALSE, vcov = FALSE)
  f4 <- krls(X, y, sigma = sv, nthreads = 4L,
             derivative = FALSE, vcov = FALSE)

  expect_identical(f1$coeffs, f4$coeffs)
  expect_identical(f1$fitted, f4$fitted)
})

test_that("ARD: predict() under ARD matches outer-loop reference", {
  set.seed(7)
  n <- 60; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(scale(X[, 1]^2 - 0.5 * X[, 2] + 0.3 * rnorm(n)))
  sv <- c(1.5, 3.0, 7.0, 12.0)

  fit <- krls(X, y, sigma = sv, derivative = FALSE, vcov = FALSE)

  Xnew <- matrix(rnorm(10 * p), 10, p)

  # Reference: standardise Xnew against training centres/scales.
  Xmeans <- colMeans(X)
  Xsd    <- apply(X, 2L, sd)
  ymean  <- mean(y); ysd <- sd(y)
  Xs    <- scale(X, center = Xmeans, scale = Xsd)
  Xs    <- matrix(Xs, n, p)
  Xn_s  <- scale(Xnew, center = Xmeans, scale = Xsd)
  Xn_s  <- matrix(Xn_s, nrow(Xnew), p)

  K_ref <- matrix(0, nrow(Xnew), n)
  for (i in seq_len(nrow(Xnew))) {
    for (j in seq_len(n)) {
      s <- sum((Xn_s[i, ] - Xs[j, ])^2 / sv)
      K_ref[i, j] <- exp(-s)
    }
  }
  # Unstandardise the prediction in raw y scale.
  yhat_ref <- as.numeric(K_ref %*% fit$coeffs) * ysd + ymean

  pr <- predict(fit, Xnew)$fit
  expect_lt(max(abs(as.numeric(pr) - yhat_ref)), 1e-10)
})

test_that("ARD: sigma validation rejects bad input and bad compositions", {
  set.seed(1)
  n <- 30; p <- 5
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  # (a) Wrong length.
  expect_error(krls(X, y, sigma = c(1, 2, 3)),
               regexp = "length 1 or ncol\\(X\\)")
  # (b) Non-positive element.
  expect_error(krls(X, y, sigma = c(1, 0, 2, 1, 1)),
               regexp = "strictly positive")
  # (c) NA element.
  expect_error(krls(X, y, sigma = c(1, NA, 2, 1, 1)),
               regexp = "finite|NA")
  # (d) ARD + autotune.
  expect_error(krls(X, y, sigma = c(1, 1, 1, 1, 2), autotune = TRUE),
               regexp = "autotune over per-feature sigma")
  # (e) ARD + Nystrom.
  expect_error(krls(X, y, sigma = c(1, 1, 1, 1, 2), approx = "nystrom"),
               regexp = "approx = 'nystrom' does not yet support per-feature sigma")
})
