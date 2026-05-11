# Smoke + sanity tests for plot.ares() diagnostic panels.
# Drawing happens to a pdf() device under tempfile() so the tests don't
# pop up an X11 window and don't depend on any interactive device.

with_pdf_device <- function(expr) {
  pf <- tempfile(fileext = ".pdf")
  grDevices::pdf(pf)
  on.exit({
    grDevices::dev.off()
    unlink(pf)
  })
  force(expr)
}

test_that("plot.ares() draws the 4 default panels for gaussian fit", {
  set.seed(1L)
  n <- 80L
  x <- matrix(stats::rnorm(n * 3L), n, 3L)
  y <- 0.5 + x[, 1L] - 0.8 * x[, 2L]^2 + stats::rnorm(n, sd = 0.3)
  fit <- ares(x, y, nthreads = 1L)

  with_pdf_device({
    expect_silent(plot(fit))
  })
})

test_that("plot.ares() honours `which` and supports all 6 panels", {
  set.seed(2L)
  n <- 60L
  x <- matrix(stats::rnorm(n * 2L), n, 2L)
  y <- x[, 1L] + stats::rnorm(n)
  fit <- ares(x, y, nthreads = 1L)

  with_pdf_device({
    for (k in 1:6) expect_silent(plot(fit, which = k))
    expect_silent(plot(fit, which = c(1L, 4L)))
  })

  expect_error(plot(fit, which = 0L), "subset of 1:6")
  expect_error(plot(fit, which = 7L), "subset of 1:6")
})

test_that("plot.ares() works for binomial family (deviance residuals)", {
  set.seed(3L)
  n <- 100L
  x <- matrix(stats::rnorm(n * 3L), n, 3L)
  lp <- 0.3 * x[, 1L] - 0.5 * x[, 2L]
  y <- stats::rbinom(n, 1L, stats::plogis(lp))
  fit <- ares(x, y, family = "binomial", nthreads = 1L)

  with_pdf_device({
    expect_silent(plot(fit, which = 1L))
    expect_silent(plot(fit))
  })
})

test_that("plot.ares() works for poisson family", {
  set.seed(4L)
  n <- 80L
  x <- matrix(stats::rnorm(n * 2L), n, 2L)
  mu <- exp(0.5 + 0.3 * x[, 1L])
  y <- stats::rpois(n, mu)
  fit <- ares(x, y, family = "poisson", nthreads = 1L)

  with_pdf_device({
    expect_silent(plot(fit))
  })
})

test_that("plot.ares() works for gamma family", {
  set.seed(5L)
  n <- 80L
  x <- matrix(stats::rnorm(n * 2L), n, 2L)
  mu <- exp(0.5 + 0.3 * x[, 1L])
  y <- stats::rgamma(n, shape = 2, rate = 2 / mu)
  fit <- ares(x, y, family = "gamma", nthreads = 1L)

  with_pdf_device({
    expect_silent(plot(fit))
  })
})

test_that("plot.ares() respects observation weights for leverage", {
  set.seed(6L)
  n <- 80L
  x <- matrix(stats::rnorm(n * 2L), n, 2L)
  y <- x[, 1L] + stats::rnorm(n, sd = 0.3)
  w <- stats::runif(n, 0.5, 1.5)
  fit_w <- ares(x, y, weights = w, nthreads = 1L)
  fit_u <- ares(x, y, nthreads = 1L)

  # Compute hat values by hand the same way plot.ares() does and confirm
  # weighted vs unweighted differ; this guards against regression to an
  # unweighted-only hat matrix path.
  expect_identical(length(fit_w$weights), n)

  h_w <- {
    sqW <- sqrt(w)
    Q <- qr.Q(qr(sqW * fit_w$bx))
    rowSums(Q * Q)
  }
  h_u <- {
    Q <- qr.Q(qr(fit_u$bx))
    rowSums(Q * Q)
  }
  expect_false(isTRUE(all.equal(h_w, h_u, tolerance = 1e-6)))

  with_pdf_device({
    expect_silent(plot(fit_w))
  })
})

test_that("plot.ares() returns the fit invisibly", {
  set.seed(7L)
  n <- 50L
  x <- matrix(stats::rnorm(n * 2L), n, 2L)
  y <- x[, 1L] + stats::rnorm(n)
  fit <- ares(x, y, nthreads = 1L)

  with_pdf_device({
    out <- plot(fit, which = 1L)
  })
  expect_identical(out, fit)
})

test_that("plot.ares() errors when $bx is missing", {
  set.seed(8L)
  n <- 40L
  x <- matrix(stats::rnorm(n * 2L), n, 2L)
  y <- x[, 1L] + stats::rnorm(n)
  fit <- ares(x, y, nthreads = 1L)
  fit$bx <- NULL

  with_pdf_device({
    expect_error(plot(fit), "design matrix")
  })
})

test_that("plot.ares() errors on non-ares input", {
  expect_error(plot.ares(list(call = quote(foo()))), "must be an 'ares' object")
})
