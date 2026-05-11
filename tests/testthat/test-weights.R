# Tests for the `weights` argument (observation weights / WLS).
#
# Quick tests:
#   - weights = rep(1, n) byte-identical to unweighted (NULL) path.
#   - weights = NULL byte-identical to omitting the argument.
#   - validation errors (length, NA, negative).
#   - basic fit succeeds with non-uniform weights and produces sensible
#     fitted values (R^2 against signal stays high on a clean DGP).
#   - bag fits compose with weights (n.boot > 0) under quick suite.
#
# Heavy tests (skip_if_quick): CV + weights, autotune + weights, determinism
# (nthreads = 1 vs N) under weights.

# Small Friedman-1-like DGP. Returns x, y, plus a smoothly-varying weight
# pattern (rows in the middle x[, 1] range get higher weight). Keeps n small
# so the quick tests stay fast.
.w_dgp <- function(n = 200, p = 5, seed = 1) {
  set.seed(seed)
  x <- matrix(stats::runif(n * p), n, p)
  f <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
       10 * x[, 4] + 5 * x[, 5]
  y <- f + stats::rnorm(n)
  # Weights: heteroscedastic — variance proportional to (0.5 + x1).
  w <- 1 / (0.5 + x[, 1])
  list(x = x, y = y, w = w, signal = f)
}

test_that("weights = rep(1, n) reproduces unweighted fit byte-for-byte", {
  d <- .w_dgp()
  f_unw <- ares(d$x, d$y, nthreads = 2)
  f_w1  <- ares(d$x, d$y, weights = rep(1, length(d$y)), nthreads = 2)
  expect_identical(f_w1$rss, f_unw$rss)
  expect_identical(f_w1$gcv, f_unw$gcv)
  expect_identical(f_w1$coefficients, f_unw$coefficients,
                   ignore_attr = TRUE)
  expect_identical(f_w1$selected.terms, f_unw$selected.terms)
  expect_identical(f_w1$dirs, f_unw$dirs)
  expect_identical(f_w1$cuts, f_unw$cuts)
})

test_that("weights = NULL identical to omitting weights", {
  d <- .w_dgp(seed = 2)
  f1 <- ares(d$x, d$y, nthreads = 2)
  f2 <- ares(d$x, d$y, weights = NULL, nthreads = 2)
  expect_equal(f1$rss, f2$rss, tolerance = 0)
  expect_equal(f1$coefficients, f2$coefficients, tolerance = 0,
               ignore_attr = TRUE)
})

test_that("weights validation: length, finiteness, negativity", {
  d <- .w_dgp()
  expect_error(ares(d$x, d$y, weights = c(1, 1), nthreads = 2),
               "length\\(weights\\)")
  bad_na <- d$w; bad_na[1] <- NA
  expect_error(ares(d$x, d$y, weights = bad_na, nthreads = 2),
               "finite")
  bad_neg <- d$w; bad_neg[1] <- -1
  expect_error(ares(d$x, d$y, weights = bad_neg, nthreads = 2),
               "non-negative")
  expect_error(ares(d$x, d$y, weights = rep(0, length(d$y)), nthreads = 2),
               "zero")
})

test_that("weighted fit produces sensible predictions on a clean DGP", {
  d <- .w_dgp(n = 400, seed = 3)
  fw <- ares(d$x, d$y, weights = d$w, degree = 2, nthreads = 2)
  # Should still fit the underlying signal well — weighted MSE shouldn't
  # explode and R^2 vs the true signal should be high.
  yhat <- predict(fw, d$x)
  ss_res <- sum((d$signal - yhat)^2)
  ss_tot <- sum((d$signal - mean(d$signal))^2)
  r2 <- 1 - ss_res / ss_tot
  expect_gt(r2, 0.85)
})

test_that("weighted fit differs from unweighted on the same DGP", {
  d <- .w_dgp(n = 300, seed = 4)
  fu <- ares(d$x, d$y, nthreads = 2)
  fw <- ares(d$x, d$y, weights = d$w, nthreads = 2)
  # Coefficients should differ on a non-trivial weight pattern.
  # (Defensive: occasionally the term-set ends up identical with very close
  # coefs; check at least one number changes by a non-trivial margin.)
  expect_false(isTRUE(all.equal(fu$coefficients, fw$coefficients,
                                tolerance = 1e-6)))
})

test_that("bagging composes with weights (quick)", {
  d <- .w_dgp(n = 200, seed = 5)
  fb <- ares(d$x, d$y, weights = d$w, n.boot = 3,
             seed.cv = 1L, nthreads = 2)
  expect_equal(length(fb$boot$fits), 3)
  yhat <- predict(fb, d$x)
  expect_true(all(is.finite(yhat)))
})

# ---- Heavy: CV + weights, autotune + weights, determinism ----
test_that("CV pruning composes with weights (heavy)", {
  skip_if_quick()
  d <- .w_dgp(n = 300, seed = 6)
  fc <- ares(d$x, d$y, weights = d$w, nfold = 3L, seed.cv = 1L,
             nthreads = 2)
  expect_identical(fc$pmethod, "cv")
  expect_true(is.finite(fc$rss))
  # cv.mse should be finite for at least one size.
  expect_true(any(is.finite(fc$cv$cv.mse)))
})

test_that("autotune composes with weights (heavy)", {
  skip_if_quick()
  d <- .w_dgp(n = 300, seed = 7)
  fa <- ares(d$x, d$y, weights = d$w, autotune = TRUE,
             seed.cv = 1L, nthreads = 2)
  expect_true(!is.null(fa$autotune))
  expect_true(is.finite(fa$rss))
})

test_that("weighted determinism: nthreads=1 == nthreads=N byte-identical (heavy)", {
  skip_if_quick()
  d <- .w_dgp(n = 250, seed = 8)
  f1 <- ares(d$x, d$y, weights = d$w, degree = 2, nthreads = 1)
  fN <- ares(d$x, d$y, weights = d$w, degree = 2, nthreads = 4)
  expect_identical(f1$rss, fN$rss)
  expect_identical(f1$coefficients, fN$coefficients, ignore_attr = TRUE)
  expect_identical(f1$selected.terms, fN$selected.terms)
})
