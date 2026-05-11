# Tests for family = "poisson" and family = "gamma".
#
# Both families share the same plumbing: forward + backward run as gaussian
# on the original y, then the selected basis is refit via stats::glm.fit()
# under the appropriate log-link GLM. Linear predictor is clamped before
# exp() to keep predictions finite under extrapolation. Tests mirror the
# existing test-family-binomial structure.

# Helper: simple poisson DGP (log-mean is linear in x1, x2). Expected counts
# kept moderate so glm.fit converges fast.
.pois_dgp <- function(n = 300, p = 4, seed = 1) {
  set.seed(seed)
  x <- matrix(stats::runif(n * p), n, p)
  eta <- 0.5 + 0.7 * x[, 1] - 0.5 * x[, 2] + 0.3 * x[, 3]
  mu  <- exp(eta)
  y   <- stats::rpois(n, lambda = mu)
  list(x = x, y = y, eta = eta, mu = mu)
}

# Helper: simple gamma DGP (log-mean is linear in x). Gamma with fixed shape.
.gamma_dgp <- function(n = 300, p = 4, seed = 1) {
  set.seed(seed)
  x <- matrix(stats::runif(n * p), n, p)
  eta <- 0.7 + 0.5 * x[, 1] - 0.4 * x[, 2] + 0.3 * x[, 3]
  mu  <- exp(eta)
  shape <- 4
  y <- stats::rgamma(n, shape = shape, rate = shape / mu)
  list(x = x, y = y, eta = eta, mu = mu)
}

# ============================================================================
# Poisson
# ============================================================================

test_that("poisson family fits and stores GLM bookkeeping", {
  d <- .pois_dgp()
  fp <- ares(d$x, d$y, family = "poisson", nthreads = 2)
  expect_identical(fp$family, "poisson")
  expect_true(all(c("deviance", "null.deviance", "aic", "converged",
                    "family", "link") %in% names(fp$glm)))
  expect_identical(fp$glm$family, "poisson")
  expect_identical(fp$glm$link, "log")
  # Fitted values should be positive (it's a count rate).
  expect_true(all(fp$fitted.values > 0))
  expect_true(fp$glm$deviance <= fp$glm$null.deviance)
})

test_that("poisson rejects non-integer / negative y", {
  d <- .pois_dgp()
  bad <- as.numeric(d$y); bad[1] <- 1.5
  expect_error(ares(d$x, bad, family = "poisson", nthreads = 2),
               "integer")
  neg <- d$y; neg[1] <- -1L
  expect_error(ares(d$x, neg, family = "poisson", nthreads = 2),
               "negative")
})

test_that("poisson predict(type=response/link) round-trip", {
  d <- .pois_dgp(seed = 2)
  fp <- ares(d$x, d$y, family = "poisson", nthreads = 2)
  p_resp <- predict(fp, d$x, type = "response")
  p_link <- predict(fp, d$x, type = "link")
  expect_equal(p_resp, exp(p_link), tolerance = 1e-10)
  expect_true(all(p_resp > 0))
})

test_that("poisson coefs == stats::glm.fit on the selected basis", {
  d <- .pois_dgp(seed = 3)
  fp <- ares(d$x, d$y, family = "poisson", nthreads = 2)
  g <- stats::glm.fit(fp$bx, d$y, family = stats::poisson(link = "log"),
                       intercept = FALSE)
  expect_equal(unname(fp$coefficients), unname(g$coefficients),
               tolerance = 1e-6)
})

# ============================================================================
# Gamma
# ============================================================================

test_that("gamma family fits and stores GLM bookkeeping", {
  d <- .gamma_dgp()
  fg <- ares(d$x, d$y, family = "gamma", nthreads = 2)
  expect_identical(fg$family, "gamma")
  expect_identical(fg$glm$family, "gamma")
  expect_identical(fg$glm$link, "log")
  expect_true(all(fg$fitted.values > 0))
})

test_that("gamma rejects y <= 0", {
  d <- .gamma_dgp()
  bad <- d$y; bad[1] <- 0
  expect_error(ares(d$x, bad, family = "gamma", nthreads = 2),
               "y > 0")
  bad2 <- d$y; bad2[1] <- -0.5
  expect_error(ares(d$x, bad2, family = "gamma", nthreads = 2),
               "y > 0")
})

test_that("gamma predict(type=response/link) round-trip", {
  d <- .gamma_dgp(seed = 4)
  fg <- ares(d$x, d$y, family = "gamma", nthreads = 2)
  p_resp <- predict(fg, d$x, type = "response")
  p_link <- predict(fg, d$x, type = "link")
  expect_equal(p_resp, exp(p_link), tolerance = 1e-8)
  expect_true(all(p_resp > 0))
})

# ============================================================================
# Compose with weights
# ============================================================================

test_that("poisson + weights composes (quick)", {
  d <- .pois_dgp(n = 200, seed = 5)
  w <- runif(200, 0.5, 2)
  fp <- ares(d$x, d$y, family = "poisson", weights = w, nthreads = 2)
  expect_identical(fp$family, "poisson")
  expect_true(all(is.finite(fp$fitted.values)))
})

test_that("gamma + weights composes (quick)", {
  d <- .gamma_dgp(n = 200, seed = 6)
  w <- runif(200, 0.5, 2)
  fg <- ares(d$x, d$y, family = "gamma", weights = w, nthreads = 2)
  expect_identical(fg$family, "gamma")
  expect_true(all(is.finite(fg$fitted.values)))
})

# ============================================================================
# Heavy: CV, autotune, bagging compose
# ============================================================================

test_that("poisson + autotune composes (heavy)", {
  skip_if_quick()
  d <- .pois_dgp(n = 300, seed = 7)
  fa <- ares(d$x, d$y, family = "poisson", autotune = TRUE,
             seed.cv = 1L, nthreads = 2)
  expect_identical(fa$family, "poisson")
  expect_true(!is.null(fa$autotune))
})

test_that("gamma + bagging composes (heavy)", {
  skip_if_quick()
  d <- .gamma_dgp(n = 250, seed = 8)
  fb <- ares(d$x, d$y, family = "gamma", n.boot = 3,
             seed.cv = 1L, nthreads = 2)
  expect_equal(length(fb$boot$fits), 3)
  p <- predict(fb, d$x)
  expect_true(all(p > 0))
})
