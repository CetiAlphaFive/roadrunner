# Phase Q2 (v0.0.0.9049): GP posterior variance via predict(type="variance").

test_that("posterior var matches brute-force K** - K*(K+lI)^{-1}K*'", {
  set.seed(11)
  n <- 40; p <- 2
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(X[, 1]) + 0.5 * X[, 2] + 0.2 * rnorm(n)
  fit <- krls(X, y, derivative = FALSE, vcov = FALSE)
  Xnew <- matrix(rnorm(15 * p), 15, p)

  v_cpp <- predict(fit, Xnew, type = "variance")

  # Brute-force reference on standardised scale.
  Xmeans <- colMeans(X)
  Xsd    <- apply(X, 2L, sd)
  Xs <- matrix(scale(X, Xmeans, Xsd), n, p)
  Xn <- matrix(scale(Xnew, Xmeans, Xsd), nrow(Xnew), p)
  sigma_vec <- fit$sigma_vec
  K <- krls_kernel_cpp(Xs, sigma_vec)
  Kstar <- krls_kernel_pred_cpp(Xn, Xs, sigma_vec)
  Kss <- rep(1.0, nrow(Xnew))
  v_brute <- Kss - diag(Kstar %*% solve(K + fit$lambda * diag(n), t(Kstar)))
  expect_lt(max(abs(v_cpp - v_brute)), 1e-8)
})

test_that("posterior var at training points = 1 - diag(H)", {
  set.seed(12)
  n <- 30; p <- 2
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(X[, 1]) + 0.5 * X[, 2] + 0.2 * rnorm(n)
  fit <- krls(X, y, derivative = FALSE, vcov = FALSE)

  v <- predict(fit, X, type = "variance")

  # Build hat-matrix diagonal: H = V diag(d/(d+lambda)) V', diag = sum_k V_ik^2 * d_k/(d_k+lam)
  K <- fit$K
  eo <- krls_eig_cpp(K)
  d  <- eo$values
  V  <- eo$vectors
  shrink <- d / (d + fit$lambda)
  h <- as.numeric((V * V) %*% shrink)
  # var at training is 1 - diag(K(K+lI)^{-1}K') in eigen-basis.
  # That equals 1 - sum_k (V_ik^2 * d_k * d_k/(d_k+lam))/(? ) ... rederive:
  # K = V D V'; K*(K+lI)^{-1}K' on training i = sum_k V_ik^2 * d_k^2/(d_k+lam)
  inner <- as.numeric((V * V) %*% (d * d / (d + fit$lambda)))
  var_ref <- 1 - inner
  var_ref[var_ref < 0] <- 0
  expect_lt(max(abs(v - var_ref)), 1e-10)
})

test_that("posterior var is non-negative and finite (FP floor at 0)", {
  set.seed(13)
  n <- 25; p <- 2
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  fit <- krls(X, y, derivative = FALSE, vcov = FALSE)
  Xnew <- matrix(rnorm(50 * p), 50, p)
  v <- predict(fit, Xnew, type = "variance")
  expect_true(all(is.finite(v)))
  expect_true(all(v >= 0))
})

test_that("posterior var increases away from training (1D toy)", {
  set.seed(14)
  n <- 20
  X <- matrix(seq(-2, 2, length.out = n), n, 1)
  y <- as.numeric(sin(X)) + 0.05 * rnorm(n)
  fit <- krls(X, y, derivative = FALSE, vcov = FALSE)
  # Test points: at training centre vs at far extrapolation.
  Xnew <- matrix(c(0, 10), 2, 1)
  v <- predict(fit, Xnew, type = "variance")
  # Variance at extrapolation should be larger than at center.
  expect_gt(v[2], v[1])
})

test_that("unscale = TRUE multiplies by var(y_train)", {
  set.seed(15)
  n <- 30; p <- 2
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(X[, 1]) + rnorm(n)
  fit <- krls(X, y, derivative = FALSE, vcov = FALSE)
  Xnew <- matrix(rnorm(10 * p), 10, p)
  v_std <- predict(fit, Xnew, type = "variance", unscale = FALSE)
  v_raw <- predict(fit, Xnew, type = "variance", unscale = TRUE)
  expect_lt(max(abs(v_raw - v_std * var(y))), 1e-12)
})

test_that("Nystrom + type='variance' errors with clear message", {
  set.seed(16)
  n <- 80; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  fit <- krls(X, y, approx = "nystrom",
              derivative = FALSE, vcov = FALSE,
              landmark_seed = 1L)
  expect_error(predict(fit, X, type = "variance"),
               regexp = "not supported for Nystrom")
})

test_that("bagging averages posterior variance across replicates", {
  skip_if_quick()
  set.seed(17)
  n <- 50; p <- 2
  X <- matrix(rnorm(n * p), n, p)
  y <- sin(X[, 1]) + 0.2 * rnorm(n)
  fit <- krls(X, y, derivative = FALSE, vcov = FALSE,
              n.boot = 4L, seed.cv = 7L)
  Xnew <- matrix(rnorm(8 * p), 8, p)
  v <- predict(fit, Xnew, type = "variance")
  expect_equal(length(v), 8L)
  expect_true(all(is.finite(v) & v >= 0))
})

test_that("poly2 posterior var diagonal = (||x*||^2 + c)^2 on standardised scale", {
  set.seed(18)
  n <- 25; p <- 2
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  cval <- 1.0
  fit <- krls(X, y, whichkernel = "poly2", poly_c = cval,
              derivative = FALSE, vcov = FALSE)
  Xnew <- matrix(rnorm(5 * p), 5, p)
  v <- predict(fit, Xnew, type = "variance")
  # Reference: kss = (||x*||^2 + c)^2 - K*(K+lI)^{-1}K*'.
  Xmeans <- colMeans(X)
  Xsd    <- apply(X, 2L, sd)
  Xs <- matrix(scale(X, Xmeans, Xsd), n, p)
  Xn <- matrix(scale(Xnew, Xmeans, Xsd), nrow(Xnew), p)
  K <- krls_kernel_cpp(Xs, fit$sigma_vec, 3L, cval)
  Kstar <- krls_kernel_pred_cpp(Xn, Xs, fit$sigma_vec, 3L, cval)
  kss <- (rowSums(Xn^2) + cval)^2
  v_ref <- kss - diag(Kstar %*% solve(K + fit$lambda * diag(n), t(Kstar)))
  v_ref[v_ref < 0] <- 0
  expect_lt(max(abs(v - v_ref)), 1e-8)
})
