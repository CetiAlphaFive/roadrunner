# Tests for causal isotonic calibration of cross-fitted binomial nuisances in
# meep() (van der Laan, Carone, Luedtke, van der Laan 2023). Calibration is ON
# by default; binomial nuisances get pooled isotonic calibration with truncation.

# Mann-Whitney AUC of scores `p` against binary labels `y` (1 = positive).
.auc_mw <- function(p, y) {
  y <- as.numeric(y)
  pos <- p[y == 1]
  neg <- p[y == 0]
  if (length(pos) == 0L || length(neg) == 0L) return(NA_real_)
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

# ---------------------------------------------------------------------------
#  Low-level calibrator helpers
# ---------------------------------------------------------------------------

test_that(".meep_isotonic_calibrator produces a monotone non-decreasing map", {
  set.seed(1)
  p <- runif(200)
  y <- rbinom(200, 1, plogis(-1 + 3 * p))
  cal <- .meep_isotonic_calibrator(p, y, c(1e-3, 1 - 1e-3))
  expect_true(all(diff(cal$yf) >= -1e-12))      # isotonic => non-decreasing
  expect_equal(cal$bounds, c(1e-3, 1 - 1e-3))
})

test_that(".meep_apply_calibrator respects bounds and is monotone in p", {
  set.seed(2)
  p <- runif(300)
  y <- rbinom(300, 1, plogis(-1 + 4 * p))
  cal <- .meep_isotonic_calibrator(p, y, c(1e-3, 1 - 1e-3))
  grid <- seq(-0.5, 1.5, length.out = 400)          # includes out-of-range
  out <- .meep_apply_calibrator(cal, grid)
  expect_true(all(out >= 1e-3 - 1e-12))
  expect_true(all(out <= 1 - 1e-3 + 1e-12))
  # monotone non-decreasing in the raw input
  expect_true(all(diff(out) >= -1e-12))
})

test_that("a perfectly-ordered (p, y) yields a sensible calibrator", {
  # y identical to a step in p: low p -> 0, high p -> 1
  p <- seq(0, 1, length.out = 100)
  y <- as.numeric(p > 0.5)
  cal <- .meep_isotonic_calibrator(p, y, c(1e-3, 1 - 1e-3))
  lo <- .meep_apply_calibrator(cal, 0.1)
  hi <- .meep_apply_calibrator(cal, 0.9)
  expect_lt(lo, hi)
  expect_gte(lo, 1e-3)
  expect_lte(hi, 1 - 1e-3)
})

# ---------------------------------------------------------------------------
#  Default meep behaviour on a binary-treatment fit
# ---------------------------------------------------------------------------

test_that("meep() calibrates the propensity by default", {
  set.seed(11)
  n <- 300; p <- 4
  X <- matrix(runif(n * p), n, p)
  D <- rbinom(n, 1, plogis(-0.5 + 1.5 * X[, 1] - X[, 2]))
  Y <- D + X[, 1] + rnorm(n, sd = 0.5)
  fit <- meep(X, Y, treatment = D, folds = 3, seed = 1,
              learners = c("ares", "krls"))
  expect_identical(fit$calibrate, "isotonic")
  expect_false(is.null(fit$calibrators$treatment))
  # d_hat_oof strictly within the calibration bounds
  b <- fit$calibrate_bounds
  d <- fit$d_hat_oof
  expect_true(all(d >= b[1] - 1e-12))
  expect_true(all(d <= b[2] + 1e-12))
  expect_gte(min(d), b[1])
  expect_lte(max(d), b[2])
})

test_that("calibration is monotone and (near-)AUC-preserving", {
  # Isotonic calibration is a monotone non-decreasing transform of the raw
  # combined prediction, so the calibrated d_hat_oof is a non-decreasing
  # function of the raw d_hat_oof. The pooled-adjacent-violators step ties
  # adjacent points, so AUC is only weakly preserved (PAV removes rank-order
  # violations against the response, never introduces them -- so the
  # calibrated AUC is >= the raw AUC up to tie effects).
  set.seed(12)
  n <- 350; p <- 4
  X <- matrix(runif(n * p), n, p)
  D <- rbinom(n, 1, plogis(-0.5 + 2 * X[, 1] - X[, 2]))
  Y <- D + X[, 1] + rnorm(n, sd = 0.5)
  fit_cal <- meep(X, Y, treatment = D, folds = 3, seed = 7,
                  learners = c("ares", "krls"), calibrate = "isotonic")
  fit_raw <- meep(X, Y, treatment = D, folds = 3, seed = 7,
                  learners = c("ares", "krls"), calibrate = "none")
  # calibrated propensity is a monotone non-decreasing function of the raw one
  o <- order(fit_raw$d_hat_oof)
  expect_true(all(diff(fit_cal$d_hat_oof[o]) >= -1e-12))
  auc_cal <- .auc_mw(fit_cal$d_hat_oof, D)
  auc_raw <- .auc_mw(fit_raw$d_hat_oof, D)
  expect_gte(auc_cal, auc_raw - 1e-8)
})

test_that("calibrate = 'none' reproduces uncalibrated behaviour", {
  set.seed(13)
  n <- 250; p <- 4
  X <- matrix(runif(n * p), n, p)
  D <- rbinom(n, 1, plogis(-0.5 + 1.5 * X[, 1]))
  Y <- D + X[, 1] + rnorm(n, sd = 0.5)
  fit <- meep(X, Y, treatment = D, folds = 3, seed = 3,
              learners = c("ares", "krls"), calibrate = "none")
  expect_identical(fit$calibrate, "none")
  expect_null(fit$calibrators$treatment)
  expect_true(all(vapply(fit$calibrators, is.null, logical(1L))) ||
              length(fit$calibrators) == 0L)
})

test_that("gaussian-only fit creates no calibrators and is unchanged", {
  set.seed(14)
  n <- 200; p <- 4
  X <- matrix(runif(n * p), n, p)
  Y <- X[, 1] + sin(3 * X[, 2]) + rnorm(n, sd = 0.4)
  fit_def <- meep(X, Y, folds = 3, seed = 5, learners = c("ares", "krls"))
  fit_non <- meep(X, Y, folds = 3, seed = 5, learners = c("ares", "krls"),
                  calibrate = "none")
  expect_null(fit_def$calibrators$outcome)
  expect_equal(fit_def$y_hat_oof, fit_non$y_hat_oof, tolerance = 1e-10)
})

test_that("binary outcome gets a calibrated outcome nuisance", {
  set.seed(15)
  n <- 300; p <- 4
  X <- matrix(runif(n * p), n, p)
  Y <- rbinom(n, 1, plogis(-0.5 + 2 * X[, 1] - X[, 2]))
  fit <- meep(X, Y, folds = 3, seed = 9, learners = c("ares", "krls"))
  expect_false(is.null(fit$calibrators$outcome))
  b <- fit$calibrate_bounds
  expect_true(all(fit$y_hat_oof >= b[1] - 1e-12))
  expect_true(all(fit$y_hat_oof <= b[2] + 1e-12))
})

# ---------------------------------------------------------------------------
#  Determinism
# ---------------------------------------------------------------------------

test_that("calibration is deterministic across thread counts", {
  skip_if_quick()
  set.seed(21)
  n <- 300; p <- 4
  X <- matrix(runif(n * p), n, p)
  D <- rbinom(n, 1, plogis(-0.5 + 1.5 * X[, 1] - X[, 2]))
  Y <- D + X[, 1] + rnorm(n, sd = 0.5)
  RcppParallel::setThreadOptions(numThreads = 1)
  f1 <- meep(X, Y, treatment = D, folds = 3, seed = 1,
             learners = c("ares", "krls"))
  RcppParallel::setThreadOptions(numThreads = 4)
  f4 <- meep(X, Y, treatment = D, folds = 3, seed = 1,
             learners = c("ares", "krls"))
  RcppParallel::setThreadOptions(numThreads = 1)
  expect_equal(f1$d_hat_oof, f4$d_hat_oof, tolerance = 0)
  expect_equal(f1$y_hat_oof, f4$y_hat_oof, tolerance = 0)
})

# ---------------------------------------------------------------------------
#  predict.meep consistency with the stored calibrator
# ---------------------------------------------------------------------------

test_that("predict.meep applies the stored calibrator on new data", {
  set.seed(31)
  n <- 300; p <- 4
  X <- matrix(runif(n * p), n, p)
  D <- rbinom(n, 1, plogis(-0.5 + 1.5 * X[, 1] - X[, 2]))
  Y <- D + X[, 1] + rnorm(n, sd = 0.5)
  fit <- meep(X, Y, treatment = D, folds = 3, seed = 2,
              learners = c("ares", "krls"))
  pr <- predict(fit, newdata = X, nuisance = "treatment")
  b <- fit$calibrate_bounds
  expect_true(all(pr >= b[1] - 1e-12))
  expect_true(all(pr <= b[2] + 1e-12))
  # predicted calibrated propensities are monotone in the raw combined pred
  cal <- fit$calibrators$treatment
  raw_grid <- seq(min(cal$x), max(cal$x), length.out = 50)
  cal_grid <- .meep_apply_calibrator(cal, raw_grid)
  expect_true(all(diff(cal_grid) >= -1e-12))
})
