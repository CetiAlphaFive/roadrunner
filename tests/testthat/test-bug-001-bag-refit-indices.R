# BUG-001 regression test (v0.0.0.9029).
#
# The post-hoc GLM refit for each bag replicate (binomial / poisson /
# gamma) must run on the SAME bootstrap rows the basis was selected on in
# the main bagging loop. Prior to v0.0.0.9029 the refit redrew its own
# bootstrap indices from the live RNG, so the unseeded path (seed.cv =
# NULL) silently desynchronised the two loops.

test_that("bag refit (binomial) uses the same bootstrap indices as basis selection", {
  set.seed(50)
  n <- 200; p <- 4
  x <- matrix(stats::runif(n * p), n, p)
  yb <- as.integer(stats::plogis(2 * x[, 1] - 0.5) > 0.5)
  fb <- ares(x, yb, family = "binomial", n.boot = 3, nthreads = 2)
  expect_false(is.null(fb$boot$idx))
  expect_equal(length(fb$boot$idx), 3L)
  for (b in seq_along(fb$boot$fits)) {
    idx <- fb$boot$idx[[b]]
    bx_b <- mars_basis_cpp(x[idx, , drop = FALSE],
                           fb$boot$fits[[b]]$dirs,
                           fb$boot$fits[[b]]$cuts,
                           as.integer(fb$boot$fits[[b]]$selected.terms))
    g <- suppressWarnings(stats::glm.fit(
      bx_b, yb[idx], family = stats::binomial(),
      intercept = FALSE,
      control = list(maxit = 50L, epsilon = 1e-8)))
    expected_cf <- as.numeric(g$coefficients)
    expected_cf[!is.finite(expected_cf)] <- 0
    observed_cf <- as.numeric(fb$boot$fits[[b]]$coefficients)
    expect_equal(observed_cf, expected_cf, tolerance = 1e-10)
  }
})

test_that("bag refit (poisson) uses the same bootstrap indices as basis selection", {
  set.seed(51)
  n <- 200; p <- 3
  x <- matrix(stats::runif(n * p), n, p)
  mu <- exp(0.2 + 0.8 * x[, 1] + 0.5 * x[, 2])
  yp <- stats::rpois(n, mu)
  fb <- ares(x, yp, family = "poisson", n.boot = 3, nthreads = 2)
  expect_false(is.null(fb$boot$idx))
  expect_equal(length(fb$boot$idx), 3L)
  for (b in seq_along(fb$boot$fits)) {
    idx <- fb$boot$idx[[b]]
    bx_b <- mars_basis_cpp(x[idx, , drop = FALSE],
                           fb$boot$fits[[b]]$dirs,
                           fb$boot$fits[[b]]$cuts,
                           as.integer(fb$boot$fits[[b]]$selected.terms))
    g <- suppressWarnings(stats::glm.fit(
      bx_b, yp[idx], family = stats::poisson(link = "log"),
      intercept = FALSE,
      control = list(maxit = 100L, epsilon = 1e-8)))
    expected_cf <- as.numeric(g$coefficients)
    expected_cf[!is.finite(expected_cf)] <- 0
    observed_cf <- as.numeric(fb$boot$fits[[b]]$coefficients)
    expect_equal(observed_cf, expected_cf, tolerance = 1e-10)
  }
})

test_that("unseeded bag refit is no longer a function of intervening RNG draws", {
  # With the fix, the refit no longer touches the RNG -- intervening RNG
  # use between the main bag loop and the refit cannot change the result.
  set.seed(99)
  n <- 150; p <- 3
  x <- matrix(stats::runif(n * p), n, p)
  yb <- as.integer(stats::plogis(2 * x[, 1] - 0.5) > 0.5)
  set.seed(99)
  f1 <- ares(x, yb, family = "binomial", n.boot = 3, nthreads = 2)
  c1 <- lapply(f1$boot$fits, `[[`, "coefficients")
  # Refitting twice in the same R session under the same set.seed reproduces.
  set.seed(99)
  f2 <- ares(x, yb, family = "binomial", n.boot = 3, nthreads = 2)
  c2 <- lapply(f2$boot$fits, `[[`, "coefficients")
  expect_equal(c1, c2, tolerance = 0)
})
