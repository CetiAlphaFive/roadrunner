# BUG-009 regression test (v0.0.0.9032).
#
# For bagged non-gaussian fits, predict(..., type = "link") used to return
# g(mean(g^{-1}(eta_b))) -- i.e. it averaged the per-replicate inverse-
# link predictions on the response scale, then re-applied the link.
# Jensen's inequality means that's NOT mean(eta_b). Simulator showed
# median |diff| of 1.8-51 log-odds for binomial bags at moderate signal
# (p90 200-708), making predict(type="link") numerically unreliable for
# any downstream use.
#
# Fix: track per-replicate eta_b and per-replicate response separately;
# `type = "link"` returns rowMeans(etas), `type = "response"` returns
# rowMeans(resps) (unchanged from prior behaviour). Bag SD is computed
# on whichever scale was returned.

skip_if_quick <- function() {
  if (!nzchar(Sys.getenv("ARES_FULL_TESTS"))) skip("ARES_FULL_TESTS not set")
}

test_that("bagged binomial type='link' = rowMeans of per-rep linear predictors", {
  skip_if_quick()
  set.seed(11)
  n <- 300
  x <- matrix(rnorm(n * 4), n, 4)
  eta_true <- 0.6 * x[, 1] - 0.5 * x[, 2]
  y <- rbinom(n, 1, plogis(eta_true))
  fit <- ares(x, y, family = "binomial", n.boot = 5, nthreads = 2, seed.cv = 7)
  # Recompute the per-replicate eta manually from the stored boot fits.
  eta_central <- drop(mars_basis_cpp(x, fit$dirs, fit$cuts,
                                     as.integer(fit$selected.terms)) %*%
                       fit$coefficients)
  etas <- matrix(NA_real_, nrow = n, ncol = length(fit$boot$fits) + 1L)
  etas[, 1] <- eta_central
  for (b in seq_along(fit$boot$fits)) {
    fb <- fit$boot$fits[[b]]
    bx_b <- mars_basis_cpp(x, fb$dirs, fb$cuts, as.integer(fb$selected.terms))
    etas[, b + 1L] <- drop(bx_b %*% fb$coefficients)
  }
  link_expected <- rowMeans(etas)
  link_actual <- as.numeric(predict(fit, x, type = "link"))
  expect_equal(link_actual, link_expected, tolerance = 1e-10)
})

test_that("bagged poisson type='link' = rowMeans of per-rep linear predictors", {
  skip_if_quick()
  set.seed(12)
  n <- 250
  x <- matrix(rnorm(n * 3), n, 3)
  eta_true <- 0.3 * x[, 1] + 0.2 * x[, 2] + 1
  y <- rpois(n, exp(eta_true))
  fit <- ares(x, y, family = "poisson", n.boot = 4, nthreads = 2, seed.cv = 8)
  eta_central <- drop(mars_basis_cpp(x, fit$dirs, fit$cuts,
                                     as.integer(fit$selected.terms)) %*%
                       fit$coefficients)
  etas <- matrix(NA_real_, nrow = n, ncol = length(fit$boot$fits) + 1L)
  etas[, 1] <- eta_central
  for (b in seq_along(fit$boot$fits)) {
    fb <- fit$boot$fits[[b]]
    bx_b <- mars_basis_cpp(x, fb$dirs, fb$cuts, as.integer(fb$selected.terms))
    etas[, b + 1L] <- drop(bx_b %*% fb$coefficients)
  }
  link_expected <- rowMeans(etas)
  link_actual <- as.numeric(predict(fit, x, type = "link"))
  expect_equal(link_actual, link_expected, tolerance = 1e-10)
})

test_that("bagged binomial type='response' unchanged (no regression)", {
  skip_if_quick()
  set.seed(13)
  n <- 200
  x <- matrix(rnorm(n * 3), n, 3)
  eta_true <- 0.4 * x[, 1] - 0.3 * x[, 2]
  y <- rbinom(n, 1, plogis(eta_true))
  fit <- ares(x, y, family = "binomial", n.boot = 4, nthreads = 2, seed.cv = 9)
  resp <- as.numeric(predict(fit, x, type = "response"))
  expect_true(all(resp >= 0 & resp <= 1))
  # Manual: rowMeans of plogis(eta_b).
  eta_central <- drop(mars_basis_cpp(x, fit$dirs, fit$cuts,
                                     as.integer(fit$selected.terms)) %*%
                       fit$coefficients)
  ps <- matrix(NA_real_, nrow = n, ncol = length(fit$boot$fits) + 1L)
  ps[, 1] <- plogis(pmin(pmax(eta_central, -30), 30))
  for (b in seq_along(fit$boot$fits)) {
    fb <- fit$boot$fits[[b]]
    bx_b <- mars_basis_cpp(x, fb$dirs, fb$cuts, as.integer(fb$selected.terms))
    e_b <- drop(bx_b %*% fb$coefficients)
    ps[, b + 1L] <- plogis(pmin(pmax(e_b, -30), 30))
  }
  expect_equal(resp, rowMeans(ps), tolerance = 1e-10)
})

test_that("bagged se.fit on type='link' uses link-scale SD", {
  skip_if_quick()
  set.seed(14)
  n <- 150
  x <- matrix(rnorm(n * 2), n, 2)
  eta_true <- 0.5 * x[, 1]
  y <- rbinom(n, 1, plogis(eta_true))
  fit <- ares(x, y, family = "binomial", n.boot = 5, nthreads = 2, seed.cv = 10)
  out <- predict(fit, x, type = "link", se.fit = TRUE)
  sdv <- attr(out, "sd")
  expect_true(!is.null(sdv))
  expect_equal(length(sdv), n)
  expect_true(all(is.finite(sdv)))
  expect_true(all(sdv >= 0))
})
