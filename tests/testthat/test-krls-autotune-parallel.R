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

test_that("INNER-1: krls_autotune_inner_cpp matches sequential R reference", {
  set.seed(3L)
  n_tr <- 60L; n_te <- 20L; p <- 4L
  X_tr <- matrix(rnorm(n_tr * p), n_tr, p)
  X_te <- matrix(rnorm(n_te * p), n_te, p)
  y_tr <- as.numeric(rowSums(X_tr) + 0.3 * rnorm(n_tr))
  y_te <- as.numeric(rowSums(X_te) + 0.3 * rnorm(n_te))

  D_tr <- roadrunner:::krls_pairwise_sqdist_cpp(X_tr, X_tr)
  D_te <- roadrunner:::krls_pairwise_sqdist_cpp(X_te, X_tr)

  sigma_grid <- c(2, 4, 8, 16)

  ref_mse <- numeric(length(sigma_grid))
  ref_lam <- numeric(length(sigma_grid))
  for (i in seq_along(sigma_grid)) {
    s <- sigma_grid[i]
    K_tr <- exp(-D_tr / s)
    K_te <- exp(-D_te / s)
    e    <- roadrunner:::krls_eig_cpp(K_tr)
    d    <- e$values
    V    <- e$vectors
    Vty  <- as.numeric(crossprod(V, y_tr))
    Vsq  <- roadrunner:::krls_vsq_cpp(V)
    L    <- .Machine$double.eps
    q    <- which.min(abs(d - max(d) / 1000))
    while (sum(d / (d + L)) > q) L <- L * 10
    U    <- length(y_tr)
    while (sum(d / (d + U)) < 1) U <- U - 1
    gr   <- (sqrt(5) - 1) / 2
    X1 <- L + (1 - gr) * (U - L); X2 <- L + gr * (U - L)
    S1 <- roadrunner:::krls_loo_loss_cpp(d, V, Vsq, Vty, X1)
    S2 <- roadrunner:::krls_loo_loss_cpp(d, V, Vsq, Vty, X2)
    while (abs(S1 - S2) > 1e-6) {
      if (S1 < S2) { U <- X2; X2 <- X1; X1 <- L + (1 - gr) * (U - L)
                     S2 <- S1; S1 <- roadrunner:::krls_loo_loss_cpp(d, V, Vsq, Vty, X1)
      } else      { L <- X1; X1 <- X2; X2 <- L + gr * (U - L)
                     S1 <- S2; S2 <- roadrunner:::krls_loo_loss_cpp(d, V, Vsq, Vty, X2) }
    }
    lam <- (X1 + X2) / 2
    sol <- roadrunner:::krls_solve_cpp(d, V, Vsq, Vty, lam)
    alpha <- sol$coeffs
    yhat  <- as.numeric(K_te %*% alpha)
    ref_mse[i] <- mean((y_te - yhat)^2)
    ref_lam[i] <- lam
  }

  out <- roadrunner:::krls_autotune_inner_cpp(
    D_tr, D_te, y_tr, y_te, sigma_grid,
    list(tol = 1e-6, L0 = .Machine$double.eps, L_step = 10,
         U_start_from_n = TRUE),
    nthreads = 1L
  )
  expect_equal(out$mse_per_sigma, ref_mse, tolerance = 1e-9)
  expect_equal(out$lambda_per_sigma, ref_lam, tolerance = 1e-9)
})

test_that("INNER-2: sigma_grid length 1 degenerate path works", {
  set.seed(4L)
  n_tr <- 40L; n_te <- 10L; p <- 3L
  X_tr <- matrix(rnorm(n_tr * p), n_tr, p)
  X_te <- matrix(rnorm(n_te * p), n_te, p)
  y_tr <- rnorm(n_tr); y_te <- rnorm(n_te)
  D_tr <- roadrunner:::krls_pairwise_sqdist_cpp(X_tr, X_tr)
  D_te <- roadrunner:::krls_pairwise_sqdist_cpp(X_te, X_tr)
  out <- roadrunner:::krls_autotune_inner_cpp(
    D_tr, D_te, y_tr, y_te, c(5.0),
    list(tol = 1e-6, L0 = .Machine$double.eps, L_step = 10,
         U_start_from_n = TRUE),
    nthreads = 1L
  )
  expect_length(out$mse_per_sigma, 1L)
  expect_length(out$lambda_per_sigma, 1L)
  expect_true(is.finite(out$mse_per_sigma))
  expect_true(out$lambda_per_sigma > 0)
})

test_that("INNER-3: nthreads > nsigma clamps without spawning extra workers", {
  set.seed(5L)
  n_tr <- 30L; n_te <- 10L; p <- 2L
  X_tr <- matrix(rnorm(n_tr * p), n_tr, p)
  X_te <- matrix(rnorm(n_te * p), n_te, p)
  y_tr <- rnorm(n_tr); y_te <- rnorm(n_te)
  D_tr <- roadrunner:::krls_pairwise_sqdist_cpp(X_tr, X_tr)
  D_te <- roadrunner:::krls_pairwise_sqdist_cpp(X_te, X_tr)
  sg   <- c(2, 4, 8)
  out1 <- roadrunner:::krls_autotune_inner_cpp(
    D_tr, D_te, y_tr, y_te, sg,
    list(tol = 1e-6, L0 = .Machine$double.eps, L_step = 10,
         U_start_from_n = TRUE),
    nthreads = 32L
  )
  out2 <- roadrunner:::krls_autotune_inner_cpp(
    D_tr, D_te, y_tr, y_te, sg,
    list(tol = 1e-6, L0 = .Machine$double.eps, L_step = 10,
         U_start_from_n = TRUE),
    nthreads = 1L
  )
  expect_equal(out1$mse_per_sigma, out2$mse_per_sigma, tolerance = 0)
  expect_equal(out1$lambda_per_sigma, out2$lambda_per_sigma, tolerance = 0)
})

test_that("INNER-EQUIV: nthreads=1 vs nthreads=4 byte-identical on inner_cpp", {
  set.seed(6L)
  n_tr <- 80L; n_te <- 30L; p <- 4L
  X_tr <- matrix(rnorm(n_tr * p), n_tr, p)
  X_te <- matrix(rnorm(n_te * p), n_te, p)
  y_tr <- rnorm(n_tr); y_te <- rnorm(n_te)
  D_tr <- roadrunner:::krls_pairwise_sqdist_cpp(X_tr, X_tr)
  D_te <- roadrunner:::krls_pairwise_sqdist_cpp(X_te, X_tr)
  sg   <- c(1, 2, 4, 8, 16, 32)
  o1 <- roadrunner:::krls_autotune_inner_cpp(
    D_tr, D_te, y_tr, y_te, sg,
    list(tol = 1e-6, L0 = .Machine$double.eps, L_step = 10,
         U_start_from_n = TRUE),
    nthreads = 1L
  )
  o4 <- roadrunner:::krls_autotune_inner_cpp(
    D_tr, D_te, y_tr, y_te, sg,
    list(tol = 1e-6, L0 = .Machine$double.eps, L_step = 10,
         U_start_from_n = TRUE),
    nthreads = 4L
  )
  expect_identical(o1$mse_per_sigma, o4$mse_per_sigma)
  expect_identical(o1$lambda_per_sigma, o4$lambda_per_sigma)
})

test_that("EQUIV-1: byte-identical at nthreads=1 vs v0.0.0.9041 snapshot", {
  baseline_path <- testthat::test_path("..", "..", "inst", "testdata",
                                       "autotune-baseline-9041.rds")
  skip_if_not(file.exists(baseline_path),
              "autotune-baseline-9041.rds missing")
  baseline <- readRDS(baseline_path)

  set.seed(baseline$seed)
  X <- matrix(rnorm(baseline$n * baseline$p), baseline$n, baseline$p)
  y <- as.numeric(rowSums(X[, 1:2]) + 0.5 * rnorm(baseline$n))

  fit <- roadrunner::krls(
    X = X, y = y, autotune = TRUE,
    derivative = FALSE, vcov = FALSE,
    autotune.nthreads = 1L
  )

  expect_identical(fit$sigma,  baseline$sigma)
  expect_identical(fit$lambda, baseline$lambda)
  expect_identical(fit$coeffs, baseline$coeffs)
})

test_that("EQUIV-2: nthreads=1 vs nthreads=4 byte-identical (determinism contract)", {
  set.seed(7L)
  n <- 150L; p <- 4L
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(rowSums(X[, 1:2]) + 0.5 * rnorm(n))
  # Pass seed.cv so fold construction is deterministic across calls.
  fit1 <- roadrunner::krls(X, y, autotune = TRUE, derivative = FALSE,
                           vcov = FALSE, autotune.nthreads = 1L,
                           seed.cv = 7L)
  fit4 <- roadrunner::krls(X, y, autotune = TRUE, derivative = FALSE,
                           vcov = FALSE, autotune.nthreads = 4L,
                           seed.cv = 7L)
  expect_identical(fit1$sigma,  fit4$sigma)
  expect_identical(fit1$lambda, fit4$lambda)
  expect_identical(fit1$coeffs, fit4$coeffs)
})

test_that("EQUIV-3: shuffled autotune.grid input yields identical fit", {
  set.seed(8L)
  n <- 120L; p <- 4L
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(rowSums(X[, 1:2]) + 0.5 * rnorm(n))
  grid_sorted   <- c(0.5, 1, 2, 4, 8, 16)
  grid_shuffled <- sample(grid_sorted)
  fit_s <- roadrunner::krls(X, y, autotune = TRUE, derivative = FALSE,
                            vcov = FALSE, autotune.grid = grid_sorted,
                            seed.cv = 8L)
  fit_h <- roadrunner::krls(X, y, autotune = TRUE, derivative = FALSE,
                            vcov = FALSE, autotune.grid = grid_shuffled,
                            seed.cv = 8L)
  expect_identical(fit_s$sigma,  fit_h$sigma)
  expect_identical(fit_s$lambda, fit_h$lambda)
})

test_that("EQUIV-4: (ncross x nfold) variations all stable across nthreads", {
  set.seed(9L)
  n <- 120L; p <- 3L
  X <- matrix(rnorm(n * p), n, p)
  y <- rnorm(n)
  for (cfg in list(c(2, 2), c(2, 5), c(3, 10))) {
    f1 <- roadrunner::krls(X, y, autotune = TRUE, ncross = cfg[1],
                           nfold = cfg[2], derivative = FALSE, vcov = FALSE,
                           autotune.nthreads = 1L, seed.cv = 9L)
    f4 <- roadrunner::krls(X, y, autotune = TRUE, ncross = cfg[1],
                           nfold = cfg[2], derivative = FALSE, vcov = FALSE,
                           autotune.nthreads = 4L, seed.cv = 9L)
    expect_identical(f1$sigma,  f4$sigma)
    expect_identical(f1$lambda, f4$lambda)
  }
})
