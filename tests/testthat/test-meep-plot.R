# ----------------------------------------------------------------------------
#  Tests for plot.meep() -- S3 plot method for the cross-fitted ensemble
# ----------------------------------------------------------------------------

# Helper: run a plotting expression with output discarded (no file/window).
.with_null_device <- function(expr) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

# ----------------------------------------------------------------------------
#  .meep_roc unit tests
# ----------------------------------------------------------------------------

test_that(".meep_roc gives auc == 1 under perfect separation", {
  score <- c(0.1, 0.2, 0.8, 0.9)
  label <- c(0, 0, 1, 1)
  r <- .meep_roc(score, label)
  expect_equal(r$auc, 1)
})

test_that(".meep_roc matches a hand-computed Mann-Whitney AUC", {
  score <- c(0.9, 0.8, 0.3, 0.2, 0.6)
  label <- c(1, 1, 0, 0, 1)
  # positives: 0.9, 0.8, 0.6 ; negatives: 0.3, 0.2 -- all positives outrank
  # all negatives, so AUC == 1.
  r <- .meep_roc(score, label)
  # Recompute expected via Mann-Whitney rank formula directly.
  rk <- rank(score)
  n1 <- sum(label == 1)
  n0 <- sum(label == 0)
  expected <- (sum(rk[label == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  expect_equal(r$auc, expected, tolerance = 1e-12)
  expect_equal(r$auc, 1, tolerance = 1e-12)
})

test_that(".meep_roc gives auc == 0 for reversed scores", {
  score <- c(0.9, 0.8, 0.3, 0.2)
  label <- c(0, 0, 1, 1)
  r <- .meep_roc(score, label)
  expect_equal(r$auc, 0)
})

test_that(".meep_roc returns NA auc for a degenerate single-class label", {
  r <- .meep_roc(c(0.1, 0.5, 0.9), c(1, 1, 1))
  expect_true(is.na(r$auc))
})

# ----------------------------------------------------------------------------
#  .meep_r2 unit test
# ----------------------------------------------------------------------------

test_that(".meep_r2 is 1 for perfect predictions", {
  y <- c(1, 2, 3, 4, 5)
  expect_equal(.meep_r2(y, y), 1)
})

# ----------------------------------------------------------------------------
#  plot.meep -- gaussian outcome (no treatment): 2-panel path
# ----------------------------------------------------------------------------

test_that("plot.meep runs clean for a gaussian outcome-only fit", {
  set.seed(11)
  n <- 120; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] + 0.5 * X[, 2]^2 + rnorm(n, sd = 0.4)
  fit <- meep(X, y, folds = 3L, tune = "none", seed = 11)
  .with_null_device(expect_invisible(plot(fit)))
})

# ----------------------------------------------------------------------------
#  plot.meep -- binary outcome + binary treatment: ROC path (both nuisances)
# ----------------------------------------------------------------------------

test_that("plot.meep runs clean for a binary outcome + binary treatment fit", {
  set.seed(12)
  n <- 160; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  lp_d <- 0.8 * X[, 1]
  D <- rbinom(n, 1L, plogis(lp_d))
  lp_y <- 0.6 * X[, 1] - 0.5 * X[, 2] + 0.7 * D
  y <- rbinom(n, 1L, plogis(lp_y))
  fit <- suppressWarnings(
    meep(X, y, treatment = D, folds = 3L, tune = "none", seed = 12))
  .with_null_device(expect_invisible(plot(fit)))
})

# ----------------------------------------------------------------------------
#  plot.meep -- mixed families: gaussian outcome + binary treatment (3 panels)
# ----------------------------------------------------------------------------

test_that("plot.meep runs clean for gaussian outcome + binary treatment", {
  set.seed(13)
  n <- 160; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  D <- rbinom(n, 1L, plogis(0.8 * X[, 1]))
  y <- X[, 1] + 0.5 * D + rnorm(n, sd = 0.4)
  fit <- suppressWarnings(
    meep(X, y, treatment = D, folds = 3L, tune = "none", seed = 13))
  expect_equal(fit$family, "gaussian")
  expect_equal(fit$treatment_family, "binomial")
  .with_null_device(expect_invisible(plot(fit)))
})

# ----------------------------------------------------------------------------
#  plot.meep -- which = filtering
# ----------------------------------------------------------------------------

test_that("plot.meep silently drops absent nuisances and honors which=mu0", {
  set.seed(14)
  n <- 180; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  D <- rbinom(n, 1L, plogis(0.8 * X[, 1]))
  y <- X[, 1] + 0.5 * D + rnorm(n, sd = 0.4)
  fit <- suppressWarnings(
    meep(X, y, treatment = D, folds = 3L, arm_models = "auto",
         tune = "none", seed = 14))
  skip_if(is.null(fit$mu0_hat_oof), "arm models not fitted")
  # which = "mu0" alone -> arm-model nuisance, runs clean
  .with_null_device(expect_invisible(plot(fit, which = "mu0")))
  # default which includes "treatment"; for an outcome-only fit it is silently
  # dropped (no error).
  fit_oo <- meep(X, y, folds = 3L, tune = "none", seed = 14)
  .with_null_device(expect_invisible(plot(fit_oo)))
})

test_that("plot.meep errors when only absent nuisances are selected", {
  set.seed(15)
  n <- 120; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  y <- X[, 1] + rnorm(n, sd = 0.4)
  fit <- meep(X, y, folds = 3L, tune = "none", seed = 15)
  .with_null_device(
    expect_error(plot(fit, which = "treatment"), "nuisance"))
})
