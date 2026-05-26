# test-bgam-bagging.R
# T-BAG-01 through T-BAG-03 + T-REG-01: bagging behavior
# Derived from test-spec.md Sections 3.6 and 6
#
# NOTE: T-BAG-01 and T-BAG-02 BLOCK due to a bagging bug:
# .bgam_bagged_fit() assigns boot_fits[[b]] <- NULL for failed replications
# (when rpois weights are zero). In R, assigning NULL to a list element
# removes it, corrupting the boot_fits list length. Simultaneously, rpois
# draws produce zeros which are rejected by the weights validator, causing
# nearly all replications to fail, so the bag may contain 0-2 valid fits.

library(testthat)

# T-BAG-01: BLOCKED -- boot$fits length is wrong due to NULL-assignment bug
test_that("T-BAG-01: n.boot > 0 returns bag mean from predict(fit, newdata=NULL)", {
  set.seed(30)
  n <- 60
  X <- matrix(rnorm(n * 3), n, 3)
  y <- X[, 1] + rnorm(n)
  fit_bag <- bgam(X, y, autotune = FALSE, mstop = 20, nu = 0.1,
                  n.boot = 5, seed = 123, nthreads = 1)
  expect_false(is.null(fit_bag$boot))
  expect_equal(fit_bag$boot$n.boot, 5L)
  # BLOCK: length(boot$fits) should be 5 but is typically 1-2 due to NULL bug
  expect_equal(length(fit_bag$boot$fits), 5L)
})

# T-BAG-02: BLOCKED -- bag mean cannot be computed (no valid boot fits)
test_that("T-BAG-02: bagged fit differs from non-bagged fit", {
  set.seed(30)
  n <- 60
  X <- matrix(rnorm(n * 3), n, 3)
  y <- X[, 1] + rnorm(n)
  fit_bag <- bgam(X, y, autotune = FALSE, mstop = 20, nu = 0.1,
                  n.boot = 5, seed = 123, nthreads = 1)
  fit_nob <- bgam(X, y, autotune = FALSE, mstop = 20, nu = 0.1,
                  n.boot = 0, nthreads = 1)
  # BLOCK: when all boot fits are NULL, the bag mean equals the central fit
  expect_gt(max(abs(fit_bag$fitted.values - fit_nob$fitted.values)), 1e-6)
})

test_that("T-BAG-03: bagging is reproducible with same seed", {
  set.seed(30)
  n <- 60
  X <- matrix(rnorm(n * 3), n, 3)
  y <- X[, 1] + rnorm(n)
  fit_a <- bgam(X, y, autotune = FALSE, mstop = 20, nu = 0.1, n.boot = 5, seed = 42)
  fit_b <- bgam(X, y, autotune = FALSE, mstop = 20, nu = 0.1, n.boot = 5, seed = 42)
  expect_lt(max(abs(predict(fit_a) - predict(fit_b))), 1e-14)
})

# T-REG-01: BLOCKED -- NULL newdata path returns NaN when no boot fits survive
test_that("T-REG-01: NULL newdata bag-mean matches newdata=X bag-mean", {
  set.seed(30)
  n <- 60
  X <- matrix(rnorm(n * 3), n, 3)
  y <- X[, 1] + rnorm(n)
  fit_bag <- bgam(X, y, autotune = FALSE, mstop = 20, nu = 0.1,
                  n.boot = 5, seed = 123, nthreads = 1)
  # BLOCK: predict(fit_bag, newdata=NULL) returns NaN when boot fits are NULL
  expect_identical(predict(fit_bag, newdata = NULL),
                   predict(fit_bag, newdata = X))
})
