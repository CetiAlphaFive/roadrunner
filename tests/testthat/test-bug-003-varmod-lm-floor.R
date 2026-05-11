# BUG-003 regression test (v0.0.0.9029).
#
# When `varmod = "lm"` extrapolates past the training yhat range and the
# predicted MAD goes non-positive, predict() now (a) warns and (b) floors
# the PI sigma at a meaningful value (max of 5% of the const-equivalent
# sigma and half the smallest positive in-sample MAD), instead of
# silently collapsing to ~1e-12 widths.

test_that("varmod='lm' extrapolation triggers a warning and finite-width PI", {
  set.seed(60)
  n <- 300; p <- 2
  x <- matrix(stats::runif(n * p, 0, 1), n, p)
  f <- 10 * x[, 1]
  # sigma DECREASES with f -> slope of |resid| ~ yhat ends up negative
  y <- f + stats::rnorm(n, sd = pmax(0.2, 5 - 5 * x[, 1]))
  fit <- ares(x, y, varmod = "lm", nthreads = 2)
  # vm stash: floor info is captured at fit time.
  expect_false(is.null(fit$varmod$sigma_floor))
  expect_gt(fit$varmod$sigma_floor, 0)
  expect_lt(fit$varmod$sigma_floor, fit$varmod$sigma_const_floor)

  # Extrapolation grid: yhat goes far past the in-sample max
  x_big <- matrix(c(seq(0, 3, length.out = 10), rep(0.5, 10)), 10, 2)

  # Warns
  expect_warning(
    pm <- predict(fit, x_big, interval = "pint"),
    "non-positive predicted"
  )

  widths <- pm[, "upr"] - pm[, "lwr"]
  # All widths are finite and positive; the floored rows are MUCH wider
  # than the legacy ~4e-12 collapse.
  expect_true(all(is.finite(widths)))
  expect_gt(min(widths), 1e-6)
  # The floor produces an exact PI width of 2*qt(.975,df)*sigma_floor on
  # extrapolated rows.
  qq <- stats::qt(0.975, df = fit$varmod$df)
  expected_floor_width <- 2 * qq * fit$varmod$sigma_floor
  # At least the rightmost extrapolated rows should sit exactly on the floor.
  expect_true(any(abs(widths - expected_floor_width) < 1e-10))
})

test_that("in-sample varmod='lm' predictions are unchanged when MAD stays positive", {
  set.seed(61)
  n <- 300; p <- 2
  x <- matrix(stats::runif(n * p), n, p)
  f <- 2 * x[, 1] + x[, 2]
  y <- f + stats::rnorm(n, sd = 0.2)
  fit <- ares(x, y, varmod = "lm", nthreads = 2)
  # In-sample prediction: all rows should have raw MAD > 0 (the fit is
  # the same data) so no warning is emitted, no flooring kicks in.
  expect_silent(pm <- predict(fit, x, interval = "pint"))
  widths <- pm[, "upr"] - pm[, "lwr"]
  expect_true(all(widths > 0))
})
