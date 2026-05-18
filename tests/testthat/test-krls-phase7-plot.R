# Phase 7 tests: plot.krls_rr 4-panel diagnostics.

# --- 1. plot() returns invisibly the object -------------------------
test_that("Phase 7: plot() returns the fit invisibly", {
  set.seed(71L)
  X <- matrix(rnorm(40 * 2), 40, 2); y <- X[, 1] + rnorm(40, sd = 0.1)
  fit <- krls(X, y, sigma = 2, lambda = 0.3,
              derivative = FALSE, vcov = TRUE)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  res <- plot(fit)
  expect_identical(res, fit)
})

# --- 2. default plot draws 4 panels --------------------------------
test_that("Phase 7: default plot draws panels 1,2,3,5", {
  set.seed(72L)
  X <- matrix(rnorm(40 * 2), 40, 2); y <- X[, 1] + rnorm(40, sd = 0.1)
  fit <- krls(X, y, sigma = 2, lambda = 0.3,
              derivative = FALSE, vcov = TRUE)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(fit))
})

# --- 3. which = 1:6 draws all panels --------------------------------
test_that("Phase 7: which = 1:6 draws all six panels", {
  set.seed(73L)
  X <- matrix(rnorm(40 * 2), 40, 2); y <- X[, 1] + rnorm(40, sd = 0.1)
  fit <- krls(X, y, sigma = 2, lambda = 0.3,
              derivative = FALSE, vcov = TRUE)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(fit, which = 1:6))
})

# --- 4. plot uses varmod sigma_hat when present --------------------
test_that("Phase 7: plot honours stored sigma_hat from varmod", {
  set.seed(74L)
  X <- matrix(rnorm(60 * 2), 60, 2); y <- X[, 1] + rnorm(60, sd = 0.3)
  fit <- krls(X, y, sigma = 2, lambda = 0.3,
              derivative = FALSE, vcov = TRUE,
              varmod = "const")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(fit))
})

# --- 5. plot.krls_rr rejects non-krls object -----------------------
test_that("Phase 7: plot.krls_rr rejects non-krls input", {
  expect_error(plot.krls_rr(list()), "krls")
})

# --- 6. invalid `which` errors --------------------------------------
test_that("Phase 7: invalid which errors", {
  set.seed(75L)
  X <- matrix(rnorm(40 * 2), 40, 2); y <- X[, 1] + rnorm(40, sd = 0.1)
  fit <- krls(X, y, sigma = 2, lambda = 0.3,
              derivative = FALSE, vcov = TRUE)
  expect_error(plot(fit, which = 7), "1:6")
})
