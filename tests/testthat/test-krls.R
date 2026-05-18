# test-krls.R -- parity tests for roadrunner::krls vs KRLS::krls.
#
# Validates that fits agree to numerical precision when sigma and lambda
# are matched, then validates the auto-lambda golden search agrees to
# within tolerance on the same DGP.  All comparisons happen on the
# *unstandardised* (raw y / X) scale to catch any scale-unstandardisation
# bug.
#
# Skipped if KRLS isn't installed (declared as Suggests in DESCRIPTION).

skip_if_no_krls <- function() {
  testthat::skip_if_not_installed("KRLS")
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
# Parity at matched (sigma, lambda): coefs / fitted / Looe.
# ----------------------------------------------------------------------
test_that("krls coefs / fitted / Looe match KRLS at matched (sigma, lambda)", {
  skip_if_no_krls()
  d <- make_dgp(n = 60L, p = 3L)
  sigma <- ncol(d$X)
  lambda <- 0.5
  ref <- KRLS::krls(d$X, d$y, sigma = sigma, lambda = lambda,
                    derivative = FALSE, vcov = TRUE, print.level = 0)
  fit <- roadrunner::krls(d$X, d$y, sigma = sigma, lambda = lambda,
                          derivative = FALSE, vcov = TRUE,
                          print.level = 0)
  expect_equal(as.numeric(fit$coeffs), as.numeric(ref$coeffs),
               tolerance = 1e-9)
  expect_equal(as.numeric(fit$fitted), as.numeric(ref$fitted),
               tolerance = 1e-9)
  expect_equal(as.numeric(fit$Looe), as.numeric(ref$Looe),
               tolerance = 1e-9)
  expect_equal(fit$sigma, ref$sigma)
  expect_equal(fit$lambda, ref$lambda)
})

# ----------------------------------------------------------------------
# Parity of automatic lambda search.
# ----------------------------------------------------------------------
test_that("auto lambda search agrees with KRLS::krls", {
  skip_if_no_krls()
  d <- make_dgp(n = 60L, p = 3L, seed = 2L)
  ref <- KRLS::krls(d$X, d$y, derivative = FALSE, vcov = TRUE,
                    print.level = 0)
  fit <- roadrunner::krls(d$X, d$y, derivative = FALSE, vcov = TRUE,
                          print.level = 0)
  ## lambda search uses tol = 1e-3 * n which is fairly generous; allow
  ## differences within tol.
  expect_equal(fit$lambda, ref$lambda, tolerance = 1e-3 * nrow(d$X))
  expect_equal(as.numeric(fit$coeffs), as.numeric(ref$coeffs),
               tolerance = 1e-6)
  expect_equal(as.numeric(fit$fitted), as.numeric(ref$fitted),
               tolerance = 1e-6)
})

# ----------------------------------------------------------------------
# Derivatives parity (pointwise and average).
# ----------------------------------------------------------------------
test_that("marginal effects match KRLS at matched lambda", {
  skip_if_no_krls()
  d <- make_dgp(n = 50L, p = 3L, seed = 3L)
  sigma <- ncol(d$X); lambda <- 0.3
  ref <- KRLS::krls(d$X, d$y, sigma = sigma, lambda = lambda,
                    derivative = TRUE, vcov = TRUE, binary = FALSE,
                    print.level = 0)
  fit <- roadrunner::krls(d$X, d$y, sigma = sigma, lambda = lambda,
                          derivative = TRUE, vcov = TRUE,
                          binary = FALSE, print.level = 0)
  expect_equal(as.numeric(fit$derivatives),
               as.numeric(ref$derivatives), tolerance = 1e-8)
  expect_equal(as.numeric(fit$avgderivatives),
               as.numeric(ref$avgderivatives), tolerance = 1e-8)
  expect_equal(as.numeric(fit$var.avgderivatives),
               as.numeric(ref$var.avgderivatives), tolerance = 1e-8)
})

# ----------------------------------------------------------------------
# predict parity (test set).
# ----------------------------------------------------------------------
test_that("predict matches KRLS on holdout data", {
  skip_if_no_krls()
  d <- make_dgp(n = 60L, p = 3L, seed = 4L)
  d_new <- make_dgp(n = 25L, p = 3L, seed = 5L)
  sigma <- ncol(d$X); lambda <- 0.4
  ref <- KRLS::krls(d$X, d$y, sigma = sigma, lambda = lambda,
                    derivative = FALSE, vcov = TRUE, print.level = 0)
  fit <- roadrunner::krls(d$X, d$y, sigma = sigma, lambda = lambda,
                          derivative = FALSE, vcov = TRUE,
                          print.level = 0)
  ref_p <- predict(ref, newdata = d_new$X, se.fit = TRUE)
  fit_p <- predict(fit, newdata = d_new$X, se.fit = TRUE)
  expect_equal(as.numeric(fit_p$fit), as.numeric(ref_p$fit),
               tolerance = 1e-9)
  expect_equal(as.numeric(fit_p$se.fit), as.numeric(ref_p$se.fit),
               tolerance = 1e-9)
})

# ----------------------------------------------------------------------
# Binary first-difference handling (fdskrls parity).
# ----------------------------------------------------------------------
test_that("binary = TRUE produces same first differences as KRLS::fdskrls", {
  skip_if_no_krls()
  set.seed(11L)
  n <- 80L
  X <- cbind(rnorm(n), sample(0:1, n, replace = TRUE),
             rnorm(n))
  colnames(X) <- c("x1", "x2_bin", "x3")
  y <- 0.5 * X[, 1L] + 1.2 * X[, 2L] + 0.3 * X[, 3L]^2 +
       rnorm(n, sd = 0.1)
  sigma <- ncol(X); lambda <- 0.2
  ref <- KRLS::krls(X, y, sigma = sigma, lambda = lambda,
                    derivative = TRUE, binary = TRUE, vcov = TRUE,
                    print.level = 0)
  fit <- roadrunner::krls(X, y, sigma = sigma, lambda = lambda,
                          derivative = TRUE, binary = TRUE,
                          vcov = TRUE, print.level = 0)
  ## ref / fit derivatives column 2 should be finite differences,
  ## columns 1 and 3 remain continuous derivatives.
  expect_equal(as.numeric(fit$derivatives),
               as.numeric(ref$derivatives), tolerance = 1e-8)
  expect_equal(as.numeric(fit$avgderivatives),
               as.numeric(ref$avgderivatives), tolerance = 1e-8)
  expect_equal(as.numeric(fit$var.avgderivatives),
               as.numeric(ref$var.avgderivatives), tolerance = 1e-8)
  expect_true(fit$binaryindicator[1L, 2L])
  expect_false(fit$binaryindicator[1L, 1L])
})

# ----------------------------------------------------------------------
# Structural / smoke tests (independent of KRLS).
# ----------------------------------------------------------------------
test_that("krls() rejects bad inputs", {
  set.seed(7L)
  X <- matrix(rnorm(40), 20, 2); y <- rnorm(20)
  expect_error(roadrunner::krls(X, rep(0, 20)), "constant")
  expect_error(roadrunner::krls(cbind(X, 1), y), "constant")
  X_na <- X; X_na[1L, 1L] <- NA
  # Phase 1 changed default na.action to "impute"; NA in X now warns + imputes
  expect_warning(roadrunner::krls(X_na, y), "median-imputed")
  # explicit omit also works (warning, not error)
  expect_warning(roadrunner::krls(X_na, y, na.action = "omit"), "dropped")
  expect_error(roadrunner::krls(X, y, sigma = -1), "sigma")
  expect_error(roadrunner::krls(X, y, lambda = -0.1), "lambda")
})

test_that("predict.krls input checks", {
  d <- make_dgp(n = 30L, p = 2L)
  fit <- roadrunner::krls(d$X, d$y, sigma = 2, lambda = 0.2,
                          derivative = FALSE, vcov = TRUE,
                          print.level = 0)
  ## column count mismatch
  expect_error(predict(fit, newdata = d$X[, 1L, drop = FALSE]),
               "ncol")
  ## se.fit without vcov fit
  fit_novcov <- roadrunner::krls(d$X, d$y, sigma = 2, lambda = 0.2,
                                 derivative = FALSE, vcov = FALSE,
                                 print.level = 0)
  expect_error(predict(fit_novcov, newdata = d$X, se.fit = TRUE),
               "vcov")
})

test_that("R2 is sensible on a clean DGP", {
  d <- make_dgp(n = 200L, p = 3L, seed = 9L)
  fit <- roadrunner::krls(d$X, d$y, derivative = FALSE, vcov = FALSE,
                          print.level = 0)
  expect_gt(fit$R2, 0.85)
  expect_lt(fit$R2, 1.0)
})

test_that("predict on training data recovers fitted values", {
  d <- make_dgp(n = 40L, p = 2L, seed = 13L)
  fit <- roadrunner::krls(d$X, d$y, sigma = 2, lambda = 0.3,
                          derivative = FALSE, vcov = FALSE,
                          print.level = 0)
  p <- predict(fit, newdata = d$X)
  expect_equal(as.numeric(p$fit), as.numeric(fit$fitted),
               tolerance = 1e-10)
})
