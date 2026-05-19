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
