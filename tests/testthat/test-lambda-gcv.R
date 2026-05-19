test_that("krls_gcv_loss_cpp matches pure-R GCV from eigendecomp", {
  set.seed(1)
  n <- 40; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(scale(rnorm(n)))
  sig <- p
  K <- exp(-as.matrix(dist(X))^2 / sig)
  eo <- eigen(K, symmetric = TRUE)
  d <- eo$values; V <- eo$vectors
  Vty <- as.numeric(crossprod(V, y))
  yty <- sum(y^2)
  for (lam in c(1e-6, 1e-3, 0.1, 1, 10, 100)) {
    trH <- sum(d / (d + lam))
    num <- sum((lam / (d + lam))^2 * Vty^2) +
           max(0, yty - sum(Vty^2))
    denom <- max(1 - trH / n, 1e-8)
    ref <- num / n / denom^2
    expect_equal(krls_gcv_loss_cpp(d, Vty, yty, n, lam), ref,
                 tolerance = 1e-10)
  }
})

test_that("krls(lambda.method='gcv') returns a sane fit", {
  set.seed(2)
  n <- 60; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1]^2 - X[, 2] + 0.3 * rnorm(n)
  fit <- krls(X, y, lambda.method = "gcv", vcov = FALSE,
              derivative = FALSE)
  expect_true(is.finite(fit$lambda) && fit$lambda > 0)
  expect_true(fit$R2 > 0)
})

test_that("gcv loss is deterministic across nthreads", {
  set.seed(3)
  n <- 50; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(scale(rnorm(n)))
  f1 <- krls(X, y, lambda.method = "gcv", vcov = FALSE,
             derivative = FALSE, nthreads = 1L)
  fN <- krls(X, y, lambda.method = "gcv", vcov = FALSE,
             derivative = FALSE, nthreads = 4L)
  expect_identical(f1$lambda, fN$lambda)
  expect_identical(f1$coeffs, fN$coeffs)
})

test_that("loo path is byte-identical to v0.0.0.9043 baseline", {
  set.seed(4)
  n <- 50; p <- 3
  X <- matrix(rnorm(n * p), n, p); y <- as.numeric(scale(rnorm(n)))
  f <- krls(X, y, lambda.method = "loo", vcov = FALSE,
            derivative = FALSE)
  expect_snapshot(round(f$lambda, 12))
  expect_snapshot(round(f$coeffs[1:5], 12))
})

test_that("gcv rejected when combined with approx='nystrom'", {
  set.seed(5)
  n <- 80; p <- 3
  X <- matrix(rnorm(n * p), n, p); y <- as.numeric(scale(rnorm(n)))
  expect_error(
    krls(X, y, lambda.method = "gcv", approx = "nystrom", nystrom_m = 30),
    "not supported"
  )
})
