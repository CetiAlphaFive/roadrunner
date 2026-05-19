## Phase 2 (v0.0.0.9043) tests: Nystrom approximation.
## Mirrors inst/specs/2026-05-19-krls-nystrom-design.md.

test_that("LAND-RAND-1: NULL landmarks with random method returns valid indices", {
  set.seed(1L)
  n <- 100L; d <- 3L
  X <- matrix(rnorm(n * d), n, d)
  out <- roadrunner:::.resolve_landmarks(
    landmarks = NULL, landmark_method = "random", nystrom_m = 25L,
    X_std = X, X_centers = rep(0, d), X_scales = rep(1, d),
    landmark_seed = NULL
  )
  expect_equal(nrow(out$matrix), 25L)
  expect_equal(ncol(out$matrix), d)
  expect_length(out$indices, 25L)
  expect_true(all(out$indices >= 1L & out$indices <= n))
  expect_false(anyDuplicated(out$indices) > 0)
  expect_identical(out$method_used, "random")
})

test_that("LAND-RAND-SEED-1: same landmark_seed produces identical indices", {
  set.seed(99L)
  n <- 100L; d <- 3L
  X <- matrix(rnorm(n * d), n, d)
  out1 <- roadrunner:::.resolve_landmarks(NULL, "random", 25L, X,
                                          rep(0, d), rep(1, d),
                                          landmark_seed = 42L)
  out2 <- roadrunner:::.resolve_landmarks(NULL, "random", 25L, X,
                                          rep(0, d), rep(1, d),
                                          landmark_seed = 42L)
  expect_identical(out1$indices, out2$indices)
  expect_identical(out1$matrix,  out2$matrix)
})

test_that("LAND-IDX-1: integer vector input returns matching X rows", {
  set.seed(2L)
  n <- 50L; d <- 4L
  X <- matrix(rnorm(n * d), n, d)
  idx <- c(3L, 7L, 12L, 25L, 40L)
  out <- roadrunner:::.resolve_landmarks(idx, "random", 5L, X,
                                         rep(0, d), rep(1, d), NULL)
  expect_identical(out$indices, idx)
  expect_equal(out$matrix, X[idx, , drop = FALSE])
  expect_identical(out$method_used, "user_indices")

  # Out-of-range indices error
  expect_error(
    roadrunner:::.resolve_landmarks(c(0L, 5L), "random", 2L, X,
                                    rep(0, d), rep(1, d), NULL),
    "unique integers in 1"
  )
})

test_that("LAND-MAT-1: numeric matrix input is standardized", {
  set.seed(3L)
  n <- 50L; d <- 3L
  X <- matrix(rnorm(n * d), n, d)
  X_centers <- colMeans(X)
  X_scales  <- apply(X, 2, sd)
  X_std     <- scale(X, center = X_centers, scale = X_scales)
  attr(X_std, "scaled:center") <- NULL
  attr(X_std, "scaled:scale")  <- NULL

  # Pass landmarks in original scale; helper should standardize internally
  Z_orig <- X[c(1, 5, 10), , drop = FALSE]
  out <- roadrunner:::.resolve_landmarks(Z_orig, "random", 3L, X_std,
                                         X_centers, X_scales, NULL)
  expect_null(out$indices)
  expect_equal(out$matrix, X_std[c(1, 5, 10), , drop = FALSE],
               tolerance = 1e-12)
  expect_identical(out$method_used, "user_matrix")

  # ncol mismatch error
  expect_error(
    roadrunner:::.resolve_landmarks(matrix(0, 3, 2), "random", 3L, X_std,
                                    X_centers, X_scales, NULL),
    "ncol\\(landmarks\\) must equal"
  )
})

test_that("LAND-VALIDATE-1: nystrom_m and nystrom_eps validators reject bad input", {
  expect_error(roadrunner:::.validate_nystrom_m(NA, 100), "finite integer")
  expect_error(roadrunner:::.validate_nystrom_m(0L, 100), "1 <= nystrom_m")
  expect_error(roadrunner:::.validate_nystrom_m(101L, 100), "1 <= nystrom_m")
  expect_error(roadrunner:::.validate_nystrom_eps(0), "positive scalar")
  expect_error(roadrunner:::.validate_nystrom_eps(-1), "positive scalar")
  expect_identical(roadrunner:::.validate_nystrom_eps(1e-9), 1e-9)
})

test_that("NYS-XK-1: krls_nystrom_predict_cpp matches naive reference", {
  set.seed(10L)
  n_new <- 20L; m <- 8L; d <- 3L
  Xn <- matrix(rnorm(n_new * d), n_new, d)
  Z  <- matrix(rnorm(m * d), m, d)
  alpha <- rnorm(m)
  sigma <- 5.0

  yhat_cpp <- roadrunner:::krls_nystrom_predict_cpp(Xn, Z, alpha, sigma)

  # Naive reference
  K_new <- matrix(NA_real_, n_new, m)
  for (i in seq_len(n_new)) for (j in seq_len(m)) {
    K_new[i, j] <- exp(-sum((Xn[i, ] - Z[j, ])^2) / sigma)
  }
  yhat_ref <- as.numeric(K_new %*% alpha)
  expect_equal(yhat_cpp, yhat_ref, tolerance = 1e-10)
})

test_that("NYS-FIT-1: krls_nystrom_fit_cpp matches reference R impl", {
  # Reference impl lives at /tmp/krls-reference/R/nystrom.R
  ref_path <- "/tmp/krls-reference/R/nystrom.R"
  skip_if_not(file.exists(ref_path), "reference Nystrom impl missing")
  ref_env <- new.env()
  sys.source(ref_path, envir = ref_env)

  set.seed(12L)
  n <- 80L; d <- 4L; m <- 16L
  X <- matrix(rnorm(n * d), n, d)
  y <- as.numeric(rowSums(X[, 1:2]) + 0.3 * rnorm(n))
  Z <- X[seq(1L, n, length.out = m), , drop = FALSE]
  sigma <- 6.0; eps <- 1e-9

  lambda_args <- list(tol = 1e-6,
                      L0  = .Machine$double.eps,
                      L_step = 10.0,
                      U_start_from_n = TRUE)

  out_cpp <- roadrunner:::krls_nystrom_fit_cpp(
    X, Z, y, sigma, lambda_args, eps, compute_vcov = FALSE
  )

  # Reference
  ref_resolved <- list(matrix = Z, indices = NULL, method_used = "user")
  out_ref <- ref_env$.fit_krls_nystrom(
    X, y, sigma, ref_resolved,
    lambda = NULL, nystrom_eps = eps,
    L = NULL, U = NULL, tol = 1e-6, noisy = FALSE,
    compute_vcov = FALSE, lambda_method = "loo"
  )

  expect_equal(as.numeric(out_cpp$coeffs), as.numeric(out_ref$coeffs),
               tolerance = 1e-7)
  expect_equal(out_cpp$lambda, out_ref$lambda, tolerance = 1e-6)
})

test_that("NYS-FIT-2: byte-identical fits at fixed sigma+landmarks", {
  set.seed(13L)
  n <- 60L; d <- 3L; m <- 12L
  X <- matrix(rnorm(n * d), n, d)
  y <- rnorm(n)
  Z <- X[1:m, , drop = FALSE]
  sigma <- 4.0

  lambda_args <- list(tol = 1e-6, L0 = .Machine$double.eps,
                      L_step = 10.0, U_start_from_n = TRUE)

  o1 <- roadrunner:::krls_nystrom_fit_cpp(X, Z, y, sigma, lambda_args,
                                          1e-9, compute_vcov = FALSE)
  o2 <- roadrunner:::krls_nystrom_fit_cpp(X, Z, y, sigma, lambda_args,
                                          1e-9, compute_vcov = FALSE)
  expect_identical(o1$coeffs, o2$coeffs)
  expect_identical(o1$lambda, o2$lambda)
})

test_that("NYS-GETLM-1: get_landmarks returns landmarks in requested scale", {
  # Construct a minimal fake Nystrom fit
  set.seed(11L)
  n <- 50L; d <- 3L
  X <- matrix(rnorm(n * d), n, d)
  X_centers <- colMeans(X)
  X_scales  <- apply(X, 2, sd)
  Z_std <- X[c(1, 5, 10), , drop = FALSE]
  Z_std <- scale(Z_std, center = X_centers, scale = X_scales)
  attr(Z_std, "scaled:center") <- NULL
  attr(Z_std, "scaled:scale")  <- NULL

  fit <- list(approx = "nystrom", landmarks = Z_std,
              X_means = X_centers, X_sds = X_scales)
  class(fit) <- "krls_rr"

  z_std <- get_landmarks(fit, scale = "standardized")
  z_orig <- get_landmarks(fit, scale = "original")
  expect_equal(z_std, Z_std)
  expect_equal(z_orig,
               sweep(sweep(Z_std, 2, X_scales, `*`), 2, X_centers, `+`),
               tolerance = 1e-12)

  # Errors on exact fit
  fit2 <- list(approx = "exact"); class(fit2) <- "krls_rr"
  expect_error(get_landmarks(fit2), "not built with approx")
})

test_that("NYS-DET-1: same landmark_seed yields byte-identical fits", {
  set.seed(20L)
  n <- 120L; p <- 4L
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(rowSums(X[, 1:2]) + 0.3 * rnorm(n))
  f1 <- roadrunner::krls(X, y, approx = "nystrom",
                         nystrom_m = 30L, landmark_seed = 42L,
                         derivative = FALSE, vcov = FALSE)
  f2 <- roadrunner::krls(X, y, approx = "nystrom",
                         nystrom_m = 30L, landmark_seed = 42L,
                         derivative = FALSE, vcov = FALSE)
  expect_identical(f1$sigma,  f2$sigma)
  expect_identical(f1$lambda, f2$lambda)
  expect_identical(f1$coeffs, f2$coeffs)
})

test_that("NYS-DET-2: landmark_seed does not disturb caller's RNG state", {
  set.seed(21L)
  n <- 80L; p <- 3L
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)

  set.seed(999L)
  v_before <- rnorm(5L)
  set.seed(999L)
  fit <- roadrunner::krls(X, y, approx = "nystrom", nystrom_m = 20L,
                          landmark_seed = 42L, derivative = FALSE,
                          vcov = FALSE)
  v_after <- rnorm(5L)

  expect_identical(v_before, v_after)
})

test_that("NYS-PRED-1: predict on a Nystrom fit returns length-n_new numeric", {
  set.seed(22L)
  n <- 100L; p <- 3L
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(rowSums(X) + 0.3 * rnorm(n))
  fit <- roadrunner::krls(X, y, approx = "nystrom", nystrom_m = 25L,
                          landmark_seed = 42L, derivative = FALSE,
                          vcov = FALSE)
  Xn <- matrix(rnorm(15 * p), 15, p)
  yhat <- predict(fit, newdata = Xn)$fit
  expect_length(yhat, 15L)
  expect_true(all(is.finite(yhat)))
})

test_that("NYS-AT-1: krls(approx='nystrom', autotune=TRUE) produces finite fit", {
  set.seed(30L)
  n <- 200L; p <- 5L
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(rowSums(X[, 1:2]) + 0.3 * rnorm(n))
  fit <- roadrunner::krls(X, y, approx = "nystrom", autotune = TRUE,
                          landmark_seed = 42L, seed.cv = 42L,
                          derivative = FALSE, vcov = FALSE,
                          autotune.nthreads = 1L)
  expect_true(is.finite(fit$sigma))
  expect_true(is.finite(fit$lambda))
  expect_identical(fit$approx, "nystrom")
  expect_true(fit$nystrom_m >= 1L)
})

test_that("NYS-AT-EQUIV: autotune Nystrom nthreads=1 vs 4 byte-identical", {
  set.seed(31L)
  n <- 180L; p <- 4L
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(rowSums(X[, 1:2]) + 0.3 * rnorm(n))
  f1 <- roadrunner::krls(X, y, approx = "nystrom", autotune = TRUE,
                         landmark_seed = 42L, seed.cv = 42L,
                         derivative = FALSE, vcov = FALSE,
                         autotune.nthreads = 1L)
  f4 <- roadrunner::krls(X, y, approx = "nystrom", autotune = TRUE,
                         landmark_seed = 42L, seed.cv = 42L,
                         derivative = FALSE, vcov = FALSE,
                         autotune.nthreads = 4L)
  expect_identical(f1$sigma,  f4$sigma)
  expect_identical(f1$lambda, f4$lambda)
  expect_identical(f1$coeffs, f4$coeffs)
})

test_that("NYS-AT-FALL: autotune Nystrom with weights errors clearly", {
  set.seed(32L)
  n <- 100L; p <- 3L
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  w <- runif(n, 0.5, 2)
  expect_error(
    roadrunner::krls(X, y, weights = w, approx = "nystrom",
                     autotune = TRUE, landmark_seed = 42L),
    "weights.*not supported.*nystrom"
  )
})
