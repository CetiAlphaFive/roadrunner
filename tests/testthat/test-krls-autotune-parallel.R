## Phase 1 (v0.0.0.9042) tests: shared distance + parallel autotune.
## Mirrors inst/specs/2026-05-18-krls-speedup-phase1-design.md.

test_that("DIST-1: krls_pairwise_sqdist_cpp matches as.matrix(dist(X))^2", {
  set.seed(1L)
  X <- matrix(rnorm(30), 10, 3)
  D_cpp <- roadrunner:::krls_pairwise_sqdist_cpp(X, X)
  D_ref <- as.matrix(stats::dist(X))^2
  expect_equal(unname(D_cpp), unname(D_ref), tolerance = 1e-10)
})

test_that("DIST-2: cross-matrix (n_a x n_b) matches naive double loop", {
  set.seed(2L)
  A <- matrix(rnorm(20), 5, 4)
  B <- matrix(rnorm(32), 8, 4)
  D_cpp <- roadrunner:::krls_pairwise_sqdist_cpp(A, B)
  D_ref <- matrix(NA_real_, 5, 8)
  for (i in seq_len(5)) for (j in seq_len(8)) {
    D_ref[i, j] <- sum((A[i, ] - B[j, ])^2)
  }
  expect_equal(D_cpp, D_ref, tolerance = 1e-10)
})

test_that("DIST-3: no negative leakage on near-zero rows (FP rounding)", {
  X <- matrix(0, 4, 3) + matrix(rnorm(12) * 1e-12, 4, 3)
  D <- roadrunner:::krls_pairwise_sqdist_cpp(X, X)
  expect_true(all(D >= 0))
  expect_true(all(diag(D) == 0 | diag(D) < 1e-20))
})
