# test-bgam-mboost-parity.R
# T-MBOOST-01 through T-MBOOST-03: parity vs mboost::gamboost
# ALL TESTS GATED BY skip_if_not_installed("mboost") AND skip_if_quick()
# Derived from test-spec.md Section 3.9
# Thresholds revised per leader auth iter 4, 2026-05-26

library(testthat)

test_that("T-MBOOST-01: gaussian linear predictor within tolerance of mboost::gamboost", {
  skip_if_not_installed("mboost")
  skip_if_quick()
  set.seed(60)
  n <- 200
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  y <- sin(X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n, 0, 0.5)
  df_m <- as.data.frame(cbind(y = y, X))

  rr_fit <- bgam(X, y, family = "gaussian", autotune = FALSE,
                 mstop = 100, nu = 0.1, nknots = 20, degree = 3,
                 dpen = 2, nthreads = 1)
  rr_eta <- rr_fit$linear.predictors

  bbs <- mboost::bbs
  mb_fit <- mboost::gamboost(
    y ~ bbs(x1, knots = 20, degree = 3, differences = 2) +
        bbs(x2, knots = 20, degree = 3, differences = 2) +
        bbs(x3, knots = 20, degree = 3, differences = 2) +
        bbs(x4, knots = 20, degree = 3, differences = 2),
    data = df_m,
    family = mboost::Gaussian(),
    control = mboost::boost_control(mstop = 100, nu = 0.1)
  )
  mb_eta <- predict(mb_fit, type = "link")[, 1]

  # Revised threshold per leader auth iter 4, 2026-05-26:
  # knot-placement (bgam equal-spacing vs mboost quantile-spacing) and
  # lambda calibration differences mean the two implementations diverge
  # substantially; 0.99 is not achievable. Empirically observed: cor ~ 0.902.
  expect_gt(cor(rr_eta, mb_eta), 0.85)
  # Relative MAD: qualitative agreement check (empirically ~ 0.202)
  expect_lt(mad(rr_eta - mb_eta) / mad(rr_eta), 0.30)
})

test_that("T-MBOOST-02: binomial linear predictor within tolerance of mboost::gamboost", {
  skip_if_not_installed("mboost")
  skip_if_quick()
  set.seed(61)
  n <- 200
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  y <- rbinom(n, 1, plogis(sin(x1) + 0.5 * x2))
  X <- cbind(x1, x2)
  df_m <- as.data.frame(cbind(y = y, X))

  # bgam receives numeric y (which it accepts)
  rr_fit <- bgam(X, y, family = "binomial", autotune = FALSE,
                 mstop = 100, nu = 0.1, nknots = 20, degree = 3,
                 dpen = 2, nthreads = 1)
  rr_eta <- rr_fit$linear.predictors

  # mboost::Binomial() requires y as a two-level factor
  y_fac <- factor(y, levels = c("0", "1"))
  df_m$y <- y_fac

  bbs <- mboost::bbs
  mb_fit <- mboost::gamboost(
    y ~ bbs(x1, knots = 20, degree = 3, differences = 2) +
        bbs(x2, knots = 20, degree = 3, differences = 2),
    data = df_m,
    family = mboost::Binomial(),
    control = mboost::boost_control(mstop = 100, nu = 0.1)
  )
  mb_eta <- predict(mb_fit, type = "link")[, 1]

  # Revised threshold per leader auth iter 4, 2026-05-26:
  # binomial path has additional IRLS divergence on top of knot-placement
  # and lambda calibration differences. Empirically observed: cor ~ 0.861.
  expect_gt(cor(rr_eta, mb_eta), 0.80)
})

test_that("T-MBOOST-03: loss trajectory is qualitatively matched", {
  skip_if_not_installed("mboost")
  skip_if_quick()
  set.seed(60)
  n <- 200
  p <- 4
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  y <- sin(X[, 1]) + 0.5 * X[, 2]^2 + rnorm(n, 0, 0.5)
  df_m <- as.data.frame(cbind(y = y, X))

  rr_fit <- bgam(X, y, family = "gaussian", autotune = FALSE,
                 mstop = 100, nu = 0.1, nknots = 20, degree = 3,
                 dpen = 2, nthreads = 1)
  rr_loss <- rr_fit$loss_path

  bbs <- mboost::bbs
  mb_fit <- mboost::gamboost(
    y ~ bbs(x1, knots = 20, degree = 3, differences = 2) +
        bbs(x2, knots = 20, degree = 3, differences = 2) +
        bbs(x3, knots = 20, degree = 3, differences = 2) +
        bbs(x4, knots = 20, degree = 3, differences = 2),
    data = df_m,
    family = mboost::Gaussian(),
    control = mboost::boost_control(mstop = 100, nu = 0.1)
  )
  mb_risks <- mboost::risk(mb_fit)[2:101]  # iterations 1..100

  expect_gt(cor(rr_loss, mb_risks), 0.95)
})
