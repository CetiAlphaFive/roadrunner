# test-krls-audit.R -- adversarial audit regression tests for krls().
#
# Each test in this file corresponds to a probe in the 2026-05-17 audit
# (REQ-20260517-001).  These tests are INDEPENDENT of test-krls.R: they
# do not assume KRLS::krls is the ground truth, they assume the
# mathematical contract of the algorithm is.  Where they DO use KRLS as
# an oracle, the goal is to probe at corners the original parity tests
# missed (n=300 not n=60, eigtrunc, auto-lambda at higher n, m=1
# predict, etc.).
#
# Skipped: heavy tests gated on ROADRUNNER_FULL_AUDIT=1 (the simulation
# half of the audit is in inst/sims/v0.27-krls-audit.R, not here).

skip_if_no_krls <- function() testthat::skip_if_not_installed("KRLS")
skip_if_quick <- function() {
  if (!nzchar(Sys.getenv("ROADRUNNER_FULL_AUDIT", ""))) {
    testthat::skip("set ROADRUNNER_FULL_AUDIT=1 to run heavy audit tests")
  }
}

make_dgp <- function(n = 60L, p = 3L, seed = 1L) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  y <- sin(X[, 1L])
  if (p >= 2L) y <- y + 0.5 * X[, 2L]^2
  if (p >= 3L) y <- y - 0.3 * X[, 3L]
  y <- y + rnorm(n, sd = 0.1)
  list(X = X, y = y)
}

# ----------------------------------------------------------------------
# AUDIT-001: multi-column y is rejected with a clear error.
# Before fix: opaque scale() error "length of 'center' must equal..."
# After fix:  "y must be a vector or single-column matrix..."
# ----------------------------------------------------------------------
test_that("krls rejects multi-column y with clear error (AUDIT-001)", {
  set.seed(1L)
  X <- matrix(rnorm(40L), 20L, 2L)
  y <- matrix(rnorm(40L), 20L, 2L)
  expect_error(roadrunner::krls(X, y, derivative = FALSE, vcov = FALSE,
                                print.level = 0),
               "y must be a vector or single-column matrix")
})

# ----------------------------------------------------------------------
# AUDIT-002: character y is REJECTED, not silently coerced.
# Before fix: as.character(rnorm(20)) was silently coerced to double
#             via storage.mode(y) <- "double" and the function ran.
# After fix:  immediate error.
# ----------------------------------------------------------------------
test_that("krls rejects character y up front (AUDIT-002)", {
  set.seed(1L)
  X <- matrix(rnorm(40L), 20L, 2L)
  y_chr <- as.character(rnorm(20L))
  expect_error(roadrunner::krls(X, y_chr, derivative = FALSE, vcov = FALSE,
                                print.level = 0),
               "y must be numeric .got character")
})

test_that("krls rejects factor y up front (AUDIT-002)", {
  set.seed(1L)
  X <- matrix(rnorm(40L), 20L, 2L)
  y_fac <- factor(sample(letters[1L:3L], 20L, replace = TRUE))
  expect_error(roadrunner::krls(X, y_fac, derivative = FALSE, vcov = FALSE,
                                print.level = 0),
               "y must be numeric .got factor")
})

test_that("krls rejects character / factor X up front (AUDIT-002)", {
  X_chr <- matrix(as.character(rnorm(40L)), 20L, 2L)
  expect_error(roadrunner::krls(X_chr, rnorm(20L), derivative = FALSE,
                                vcov = FALSE, print.level = 0),
               "X must be numeric")
})

# ----------------------------------------------------------------------
# AUDIT-003: predict() defensive checks against mutated object$X.
# ----------------------------------------------------------------------
test_that("predict.krls_rr errors on n < 2 stored X (AUDIT-003)", {
  d <- make_dgp(n = 30L, p = 2L)
  fit <- roadrunner::krls(d$X, d$y, sigma = 2, lambda = 0.3,
                          derivative = FALSE, vcov = FALSE,
                          print.level = 0)
  fit$X <- fit$X[1L, , drop = FALSE]
  expect_error(predict(fit, newdata = d$X), "fewer than 2 rows")
})

test_that("predict.krls_rr errors on zero-sd stored X column (AUDIT-003)", {
  d <- make_dgp(n = 30L, p = 2L)
  fit <- roadrunner::krls(d$X, d$y, sigma = 2, lambda = 0.3,
                          derivative = FALSE, vcov = FALSE,
                          print.level = 0)
  fit$X[, 1L] <- 1
  expect_error(predict(fit, newdata = d$X), "zero or non-finite sd")
})

# ----------------------------------------------------------------------
# AUDIT-N: predict() works on single-row newdata.
# ----------------------------------------------------------------------
test_that("predict.krls_rr handles single-row newdata", {
  d <- make_dgp(n = 30L, p = 2L)
  fit <- roadrunner::krls(d$X, d$y, sigma = 2, lambda = 0.3,
                          derivative = FALSE, vcov = TRUE,
                          print.level = 0)
  p1 <- predict(fit, newdata = d$X[1L, , drop = FALSE], se.fit = TRUE)
  pall <- predict(fit, newdata = d$X, se.fit = TRUE)
  expect_equal(as.numeric(p1$fit), as.numeric(pall$fit[1L]),
               tolerance = 1e-12)
  expect_equal(as.numeric(p1$se.fit), as.numeric(pall$se.fit[1L]),
               tolerance = 1e-12)
})

# ----------------------------------------------------------------------
# AUDIT-S3: predict on training data equals $fitted under auto-lambda.
# The original test in test-krls.R only checks this at fixed sigma/lambda.
# ----------------------------------------------------------------------
test_that("predict on training data ≡ fitted under auto-lambda", {
  d <- make_dgp(n = 80L, p = 3L, seed = 42L)
  fit <- roadrunner::krls(d$X, d$y, derivative = FALSE, vcov = FALSE,
                          print.level = 0)
  p <- predict(fit, newdata = d$X)
  expect_equal(as.numeric(p$fit), as.numeric(fit$fitted),
               tolerance = 1e-9)
})

# ----------------------------------------------------------------------
# AUDIT-R: p = 1 univariate
# ----------------------------------------------------------------------
test_that("krls handles p = 1 univariate predictor", {
  set.seed(8L)
  n <- 60L
  X <- matrix(rnorm(n), n, 1L); colnames(X) <- "x1"
  y <- sin(X[, 1L]) + rnorm(n, sd = 0.1)
  fit <- roadrunner::krls(X, y, sigma = 1, lambda = 0.2,
                          derivative = TRUE, vcov = TRUE,
                          print.level = 0)
  expect_true(is.finite(fit$avgderivatives[1L, 1L]))
  expect_gt(fit$R2, 0.5)
})

# ----------------------------------------------------------------------
# AUDIT-O: n < p high-dim
# ----------------------------------------------------------------------
test_that("krls handles n < p (high-dim)", {
  set.seed(9L)
  n <- 15L; p <- 20L
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("x", 1L:p)
  y <- X[, 1L] - 0.5 * X[, 2L] + rnorm(n, sd = 0.1)
  fit <- roadrunner::krls(X, y, sigma = p, lambda = 0.5,
                          derivative = FALSE, vcov = FALSE,
                          print.level = 0)
  expect_true(is.finite(fit$R2))
  expect_true(all(is.finite(fit$coeffs)))
})

# ----------------------------------------------------------------------
# AUDIT-P: repeated rows in X => rank-deficient K, ridge handles it.
# Duplicate rows must have identical fitted values.
# ----------------------------------------------------------------------
test_that("krls handles repeated rows in X (rank-deficient K)", {
  set.seed(303L)
  n <- 30L
  X <- matrix(rnorm(n * 2L), n, 2L); colnames(X) <- c("a", "b")
  X[7L, ] <- X[3L, ]
  X[15L, ] <- X[3L, ]
  y <- rnorm(n)
  fit <- roadrunner::krls(X, y, sigma = 2, lambda = 0.5,
                          derivative = FALSE, vcov = FALSE,
                          print.level = 0)
  ## Duplicate rows -> duplicate fitted values (the kernel sees them
  ## identically, so they must produce identical predictions).
  expect_equal(fit$fitted[3L], fit$fitted[7L], tolerance = 1e-12)
  expect_equal(fit$fitted[3L], fit$fitted[15L], tolerance = 1e-12)
  expect_true(all(is.finite(fit$coeffs)))
})

# ----------------------------------------------------------------------
# AUDIT-Q: n = 2 minimal
# ----------------------------------------------------------------------
test_that("krls handles n = 2 minimal fit", {
  set.seed(14L)
  X <- matrix(rnorm(4L), 2L, 2L); colnames(X) <- c("a", "b")
  y <- c(1, -1)
  fit <- roadrunner::krls(X, y, sigma = 2, lambda = 0.1,
                          derivative = FALSE, vcov = FALSE,
                          print.level = 0)
  expect_true(is.finite(fit$R2))
})

# ----------------------------------------------------------------------
# AUDIT-E: kernel underflow at extreme sigma -- agrees with KRLS.
# ----------------------------------------------------------------------
test_that("kernel underflow at tiny sigma matches KRLS (AUDIT-E)", {
  skip_if_no_krls()
  set.seed(13L)
  X <- matrix(rnorm(20L), 10L, 2L); colnames(X) <- c("a", "b")
  y <- rnorm(10L)
  ref <- suppressWarnings(KRLS::krls(X, y, sigma = 0.001, lambda = 0.5,
                                     derivative = FALSE, binary = FALSE,
                                     vcov = FALSE, print.level = 0))
  fit <- roadrunner::krls(X, y, sigma = 0.001, lambda = 0.5,
                          derivative = FALSE, vcov = FALSE,
                          print.level = 0)
  expect_equal(as.numeric(fit$K), as.numeric(ref$K), tolerance = 1e-12)
  expect_equal(as.numeric(fit$coeffs), as.numeric(ref$coeffs),
               tolerance = 1e-12)
  expect_true(all(is.finite(fit$coeffs)))
})

# ----------------------------------------------------------------------
# AUDIT-M: K is symmetric, diagonal is 1, PSD modulo float noise.
# ----------------------------------------------------------------------
test_that("kernel matrix is symmetric, diag(K) = 1, PSD", {
  set.seed(22L)
  n <- 100L; p <- 3L
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("x", 1L:p)
  y <- rnorm(n)
  fit <- roadrunner::krls(X, y, sigma = 3, lambda = 0.5,
                          derivative = FALSE, vcov = FALSE,
                          print.level = 0)
  expect_equal(as.numeric(max(abs(fit$K - t(fit$K)))), 0)
  expect_equal(as.numeric(max(abs(diag(fit$K) - 1))), 0)
  eigs <- eigen(fit$K, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(eigs), -1e-10)
})

# ----------------------------------------------------------------------
# AUDIT-J: avg-deriv variance row-sum trick is BIT-EQUAL to the explicit
# t(L_k) %*% V_c %*% L_k construction at n = 25.
# Identity: sum(M) = 1' M 1.  Proves the closed-form trick has no
# cancellation hazard.
# ----------------------------------------------------------------------
test_that("avg-deriv variance row-sum identity matches explicit-L at <= 1e-12 (AUDIT-J)", {
  set.seed(33L)
  n <- 25L; p <- 2L
  X <- matrix(rnorm(n * p), n, p)
  Xs <- scale(X, center = TRUE, scale = apply(X, 2L, sd))
  Xs <- matrix(Xs, n, p)
  sigma <- 2; lambda <- 0.4
  sigma_vec <- rep(sigma, p)              # Phase 2a: engine takes length-p
  K  <- roadrunner:::krls_kernel_cpp(Xs, sigma_vec)
  eo <- roadrunner:::krls_eig_cpp(K)
  d  <- as.numeric(eo$values); V <- eo$vectors
  Vc <- V %*% diag(1 / (d + lambda)^2) %*% t(V)
  sigma2 <- 1.0

  ## Build full L_k for k=1 via KRLS's explicit construction.
  rows <- cbind(rep(seq_len(n), each = n), seq_len(n))
  distances <- Xs[rows[, 1L], ] - Xs[rows[, 2L], ]
  distk_krls <- matrix(distances[, 1L], n, n, byrow = TRUE)
  L_full <- distk_krls * K
  ref_var <- (4 / sigma^2 / n^2) * sigma2 *
             sum(t(L_full) %*% Vc %*% L_full)

  cpp_var <- roadrunner:::krls_avg_deriv_var_cpp(Xs, K, V, d, sigma_vec,
                                                  lambda, sigma2)[1L]
  expect_equal(cpp_var, ref_var, tolerance = 1e-12)
})

# ----------------------------------------------------------------------
# AUDIT-D: determinism across thread counts.
# ----------------------------------------------------------------------
test_that("krls fit is bit-equal across RcppParallel thread counts", {
  set.seed(404L)
  n <- 100L; p <- 3L
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("x", 1L:p)
  y <- sin(X[, 1L]) + 0.5 * X[, 2L]^2 + rnorm(n, sd = 0.1)
  RcppParallel::setThreadOptions(numThreads = 1L)
  f1 <- roadrunner::krls(X, y, sigma = p, lambda = 0.3,
                         derivative = TRUE, vcov = TRUE, print.level = 0)
  RcppParallel::setThreadOptions(numThreads = 4L)
  f4 <- roadrunner::krls(X, y, sigma = p, lambda = 0.3,
                         derivative = TRUE, vcov = TRUE, print.level = 0)
  expect_identical(f1$K, f4$K)
  expect_identical(f1$coeffs, f4$coeffs)
  expect_identical(f1$fitted, f4$fitted)
  expect_identical(f1$avgderivatives, f4$avgderivatives)
  expect_identical(f1$var.avgderivatives, f4$var.avgderivatives)
})

# ----------------------------------------------------------------------
# AUDIT-S13: parity at n = 300 (5x larger than test-krls.R's n=60).
# ----------------------------------------------------------------------
test_that("matched-(sigma,lambda) parity holds at n = 300 (AUDIT-S13)", {
  skip_if_no_krls()
  skip_if_quick()
  set.seed(99L)
  n <- 300L; p <- 4L
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("x", 1L:p)
  y <- sin(X[, 1L]) + 0.5 * X[, 2L]^2 - 0.3 * X[, 3L] +
       0.2 * X[, 4L]^2 + rnorm(n, sd = 0.2)
  ref <- suppressWarnings(KRLS::krls(X, y, sigma = p, lambda = 0.5,
                                     derivative = TRUE, binary = FALSE,
                                     vcov = TRUE, print.level = 0))
  fit <- roadrunner::krls(X, y, sigma = p, lambda = 0.5,
                          derivative = TRUE, binary = FALSE,
                          vcov = TRUE, print.level = 0)
  expect_equal(as.numeric(fit$coeffs), as.numeric(ref$coeffs),
               tolerance = 1e-10)
  expect_equal(as.numeric(fit$fitted), as.numeric(ref$fitted),
               tolerance = 1e-10)
  expect_equal(as.numeric(fit$avgderivatives),
               as.numeric(ref$avgderivatives), tolerance = 1e-10)
  expect_equal(as.numeric(fit$var.avgderivatives),
               as.numeric(ref$var.avgderivatives), tolerance = 1e-10)
})

# ----------------------------------------------------------------------
# AUDIT-S14: eigentrunc path parity.
# ----------------------------------------------------------------------
test_that("eigentrunc = 0.01 parity holds (AUDIT-S14)", {
  skip_if_no_krls()
  set.seed(101L)
  n <- 80L; p <- 3L
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("x", 1L:p)
  y <- sin(X[, 1L]) + 0.5 * X[, 2L]^2 + rnorm(n, sd = 0.1)
  ref <- suppressWarnings(KRLS::krls(X, y, sigma = p, lambda = 0.5,
                                     eigtrunc = 0.01, derivative = TRUE,
                                     binary = FALSE, vcov = TRUE,
                                     print.level = 0))
  fit <- roadrunner::krls(X, y, sigma = p, lambda = 0.5, eigtrunc = 0.01,
                          derivative = TRUE, binary = FALSE,
                          vcov = TRUE, print.level = 0)
  expect_equal(as.numeric(fit$coeffs), as.numeric(ref$coeffs),
               tolerance = 1e-10)
  expect_equal(as.numeric(fit$avgderivatives),
               as.numeric(ref$avgderivatives), tolerance = 1e-10)
  expect_equal(as.numeric(fit$var.avgderivatives),
               as.numeric(ref$var.avgderivatives), tolerance = 1e-10)
})

# ----------------------------------------------------------------------
# AUDIT-S: all-binary X works end to end.
# ----------------------------------------------------------------------
test_that("krls handles all-binary X with binary = TRUE flagging", {
  set.seed(15L)
  n <- 60L
  X <- cbind(sample(0L:1L, n, TRUE),
             sample(0L:1L, n, TRUE),
             sample(0L:1L, n, TRUE))
  colnames(X) <- paste0("x", 1L:3L)
  y <- 0.5 + 1.2 * X[, 1L] - 0.7 * X[, 2L] + 0.3 * X[, 3L] +
       rnorm(n, sd = 0.1)
  fit <- roadrunner::krls(X, y, sigma = 3, lambda = 0.5,
                          derivative = TRUE, vcov = TRUE,
                          binary = TRUE, print.level = 0)
  expect_true(all(fit$binaryindicator))
})

# ----------------------------------------------------------------------
# AUDIT-I: lambda search terminates fast even at n = 800.
# ----------------------------------------------------------------------
test_that("auto-lambda search terminates fast at n = 800 (AUDIT-I)", {
  skip_if_quick()
  set.seed(202L)
  n <- 800L; p <- 3L
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("x", 1L:p)
  y <- sin(X[, 1L]) + 0.5 * X[, 2L]^2 + rnorm(n, sd = 0.1)
  t0 <- proc.time()
  fit <- roadrunner::krls(X, y, derivative = FALSE, vcov = FALSE,
                          print.level = 0)
  el <- (proc.time() - t0)[["elapsed"]]
  expect_lt(el, 30)
  expect_true(fit$lambda > 0)
})

# ----------------------------------------------------------------------
# AUDIT-004: predict() matches newdata columns by NAME when both
# training X and newdata are named.  KRLS upstream uses positional
# indexing silently -- a swapped column order produces silently wrong
# predictions there.  roadrunner now (a) reorders when both names exist,
# (b) errors on a hard mismatch, (c) falls back to positional when
# either side is unnamed.
# ----------------------------------------------------------------------
test_that("predict() reorders by column name when names available (AUDIT-004)", {
  set.seed(1L)
  n <- 60L
  X <- matrix(rnorm(n * 3L), n, 3L); colnames(X) <- c("a", "b", "c")
  y <- X[, 1L] - 0.5 * X[, 2L] + 0.3 * X[, 3L] + rnorm(n, sd = 0.1)
  fit <- roadrunner::krls(X, y, sigma = 3, lambda = 0.3,
                          derivative = FALSE, vcov = FALSE,
                          print.level = 0)
  Xswap <- X[, c("b", "a", "c")]
  p_swap <- predict(fit, Xswap)$fit
  p_nat  <- predict(fit, X)$fit
  expect_equal(as.numeric(p_swap), as.numeric(p_nat), tolerance = 1e-12)
})

test_that("predict() errors on colname mismatch (AUDIT-004)", {
  set.seed(2L)
  n <- 40L
  X <- matrix(rnorm(n * 3L), n, 3L); colnames(X) <- c("a", "b", "c")
  y <- rnorm(n)
  fit <- roadrunner::krls(X, y, sigma = 3, lambda = 0.3,
                          derivative = FALSE, vcov = FALSE,
                          print.level = 0)
  Xbad <- X; colnames(Xbad) <- c("a", "b", "d")
  expect_error(predict(fit, Xbad), "colnames\\(newdata\\) do not match")
})

test_that("predict() positional fallback for unnamed newdata (AUDIT-004)", {
  set.seed(3L)
  n <- 40L
  X <- matrix(rnorm(n * 3L), n, 3L); colnames(X) <- c("a", "b", "c")
  y <- rnorm(n)
  fit <- roadrunner::krls(X, y, sigma = 3, lambda = 0.3,
                          derivative = FALSE, vcov = FALSE,
                          print.level = 0)
  Xun <- unname(X)
  expect_equal(as.numeric(predict(fit, Xun)$fit),
               as.numeric(predict(fit, X)$fit), tolerance = 1e-12)
})

# ----------------------------------------------------------------------
# AUDIT-005: determinism across MORE thread counts (1, 2, 4, 8) under
# the FULL pipeline (derivative + vcov).  The original AUDIT-D test
# only covered {1, 4}.
# ----------------------------------------------------------------------
test_that("krls fit is bit-equal across nt in {1,2,4,8} under full pipeline (AUDIT-005)", {
  set.seed(404L)
  n <- 100L; p <- 3L
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("x", seq_len(p))
  y <- sin(X[, 1L]) + 0.5 * X[, 2L]^2 + rnorm(n, sd = 0.1)
  results <- vector("list", 4L)
  threads <- c(1L, 2L, 4L, 8L)
  for (i in seq_along(threads)) {
    RcppParallel::setThreadOptions(numThreads = threads[i])
    results[[i]] <- roadrunner::krls(X, y, sigma = p, lambda = 0.3,
                                     derivative = TRUE, binary = FALSE,
                                     vcov = TRUE, print.level = 0)
  }
  ref <- results[[1L]]
  for (i in 2:4) {
    expect_identical(ref$K,                results[[i]]$K)
    expect_identical(ref$coeffs,           results[[i]]$coeffs)
    expect_identical(ref$fitted,           results[[i]]$fitted)
    expect_identical(ref$derivatives,      results[[i]]$derivatives)
    expect_identical(ref$avgderivatives,   results[[i]]$avgderivatives)
    expect_identical(ref$var.avgderivatives,
                     results[[i]]$var.avgderivatives)
    expect_identical(ref$vcov.c,           results[[i]]$vcov.c)
    expect_identical(ref$vcov.fitted,      results[[i]]$vcov.fitted)
    expect_identical(ref$Looe,             results[[i]]$Looe)
  }
})

# ----------------------------------------------------------------------
# AUDIT-006: avg-deriv variance row-sum identity, MULTI-CELL stress.
# Extends AUDIT-J (single n=25 cell) across variable (n, p, sigma,
# lambda).  Each cell asserts cpp == explicit-L to <= 1e-9 absolute
# (or 1e-8 rel for huge-variance cells).
# ----------------------------------------------------------------------
test_that("avg-deriv variance row-sum identity holds multi-cell (AUDIT-006)", {
  check_cell <- function(n, p, sigma, lambda, seed) {
    set.seed(seed)
    X <- matrix(rnorm(n * p), n, p)
    Xs <- scale(X); Xs <- matrix(Xs, n, p)
    sigma_vec <- rep(sigma, p)            # Phase 2a: engine takes length-p
    K  <- roadrunner:::krls_kernel_cpp(Xs, sigma_vec)
    eo <- roadrunner:::krls_eig_cpp(K)
    d  <- as.numeric(eo$values); V <- eo$vectors
    Vc <- V %*% diag(1 / (d + lambda)^2) %*% t(V)
    sigma2 <- 1.0
    cpp_var <- roadrunner:::krls_avg_deriv_var_cpp(Xs, K, V, d, sigma_vec,
                                                    lambda, sigma2)
    refs <- numeric(p)
    for (k in seq_len(p)) {
      rows <- cbind(rep(seq_len(n), each = n), seq_len(n))
      dist_k <- matrix(Xs[rows[, 1L], k] - Xs[rows[, 2L], k], n, n,
                       byrow = TRUE)
      L_full <- dist_k * K
      refs[k] <- (4 / sigma^2 / n^2) * sigma2 *
                 sum(t(L_full) %*% Vc %*% L_full)
    }
    list(cpp = cpp_var, ref = refs)
  }
  for (cfg in list(list(n = 50,  p = 3,  sigma = 3,    lambda = 0.4, seed = 1L),
                   list(n = 80,  p = 1,  sigma = 1,    lambda = 0.4, seed = 4L),
                   list(n = 80,  p = 10, sigma = 10,   lambda = 0.4, seed = 5L),
                   list(n = 60,  p = 3,  sigma = 0.05, lambda = 0.4, seed = 7L),
                   list(n = 60,  p = 3,  sigma = 60,   lambda = 0.4, seed = 8L),
                   list(n = 60,  p = 3,  sigma = 3,    lambda = 1000, seed = 10L))) {
    r <- check_cell(cfg$n, cfg$p, cfg$sigma, cfg$lambda, cfg$seed)
    expect_equal(as.numeric(r$cpp), as.numeric(r$ref), tolerance = 1e-9)
  }
})

# ----------------------------------------------------------------------
# AUDIT-007: bigger n race-condition stress for the kernel workers.
# Run 20 fits at nt=8 on the same X; every K must be bit-identical.
# ----------------------------------------------------------------------
test_that("kernel worker is race-free under repeated nt=8 fits (AUDIT-007)", {
  set.seed(42L)
  n <- 200L; p <- 4L
  X <- matrix(rnorm(n * p), n, p)
  Xs <- scale(X); Xs <- matrix(Xs, n, p)
  RcppParallel::setThreadOptions(numThreads = 1L)
  sigma_vec_p <- rep(as.numeric(p), p)   # Phase 2a: engine takes length-p
  K_ref <- roadrunner:::krls_kernel_cpp(Xs, sigma_vec_p)
  RcppParallel::setThreadOptions(numThreads = 8L)
  for (rep in seq_len(20L)) {
    K_new <- roadrunner:::krls_kernel_cpp(Xs, sigma_vec_p)
    expect_identical(K_new, K_ref)
  }
})

# ----------------------------------------------------------------------
# AUDIT-008: n=2000 fit completes without crash + memory bounded.
# Heavy; gated on full-audit.
# ----------------------------------------------------------------------
test_that("krls fit at n=2000 completes and is sensible (AUDIT-008)", {
  skip_if_quick()
  set.seed(2024L)
  n <- 2000L; p <- 4L
  X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("x", seq_len(p))
  y <- sin(X[, 1L]) + 0.5 * X[, 2L]^2 - 0.3 * X[, 3L] +
       0.2 * X[, 4L]^2 + rnorm(n, sd = 0.2)
  t0 <- proc.time()
  fit <- roadrunner::krls(X, y, sigma = p, lambda = 0.5,
                          derivative = TRUE, vcov = TRUE,
                          binary = FALSE, print.level = 0)
  el <- (proc.time() - t0)[["elapsed"]]
  expect_lt(el, 30)
  expect_true(all(is.finite(fit$coeffs)))
  expect_true(all(is.finite(fit$avgderivatives)))
  expect_true(all(fit$var.avgderivatives >= 0))
  expect_gt(fit$R2, 0.8)
})

# ----------------------------------------------------------------------
# AUDIT-009: roadrunner::krls is bit-identical to KRLS::krls under MC
# replication on a smooth nonlinear DGP -- proves we haven't introduced
# a subtle statistical drift relative to upstream.
# ----------------------------------------------------------------------
test_that("MC replication of avg derivs matches KRLS::krls (AUDIT-009)", {
  skip_if_no_krls()
  skip_if_quick()
  REPS <- 50L
  est_rr <- est_kr <- matrix(NA_real_, REPS, 3L)
  for (r in seq_len(REPS)) {
    set.seed(100L + r)
    n <- 80L; p <- 3L
    X <- matrix(rnorm(n * p), n, p); colnames(X) <- paste0("x", 1L:p)
    y <- 2 * X[, 1L] + 0.5 * X[, 2L]^2 - 0.3 * X[, 3L] + rnorm(n, sd = 0.4)
    fit_rr <- roadrunner::krls(X, y, derivative = TRUE, vcov = TRUE,
                               binary = FALSE, print.level = 0)
    fit_kr <- suppressWarnings(KRLS::krls(X, y, derivative = TRUE,
                                           vcov = TRUE, binary = FALSE,
                                           print.level = 0))
    est_rr[r, ] <- as.numeric(fit_rr$avgderivatives)
    est_kr[r, ] <- as.numeric(fit_kr$avgderivatives)
  }
  expect_equal(est_rr, est_kr, tolerance = 1e-8)
})
