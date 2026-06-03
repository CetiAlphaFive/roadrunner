# ============================================================================
#  Tests for meep() with the plda() learner (classification-only)
# ============================================================================
#
# plda() is a penalized LDA classifier. As a meep() learner it is
# binomial-only (like logreg): it is filtered out of gaussian nuisances and
# rejects observation weights (plda has no `weights` argument). plda is part
# of the DEFAULT learner set c("ares","krls","ols","logreg","plda"); per
# nuisance family-filtering then yields gaussian = {ares,krls,ols} and
# binomial = {ares,krls,logreg,plda}.

# ----------------------------------------------------------------------------
#  Helper DGPs
# ----------------------------------------------------------------------------

.meep_plda_dgp_binary <- function(n = 320, p = 4, seed = 1) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  y <- rbinom(n, 1, plogis(1.2 * X[, 1] - 0.9 * X[, 2] + 0.5 * X[, 3] - 0.4))
  list(X = X, y = y)
}

.meep_plda_dgp_gaussian <- function(n = 280, p = 4, seed = 1) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  y <- 1 + 0.8 * X[, 1] - 0.5 * X[, 2] + 0.3 * X[, 3] + rnorm(n, sd = 0.3)
  list(X = X, y = y)
}

# ----------------------------------------------------------------------------
#  Test 1 -- plda in default learners on a binary outcome
# ----------------------------------------------------------------------------

test_that("default learners include plda; binary OOF preds are finite in (0,1)", {
  d <- .meep_plda_dgp_binary(n = 320, p = 4, seed = 301)
  fit <- suppressWarnings(
    meep(d$X, d$y, folds = 4L, tune = "none", seed = 301))

  expect_s3_class(fit, "meep")
  expect_identical(fit$learners,
                   c("ares", "krls", "ols", "logreg", "plda"))
  expect_identical(fit$family, "binomial")

  # plda column present, populated, and a valid probability.
  oof <- fit$oof_matrix$outcome
  expect_true("plda" %in% colnames(oof))
  expect_true(all(is.finite(oof[, "plda"])))
  expect_true(all(oof[, "plda"] >= 0 & oof[, "plda"] <= 1))

  # plda is a binomial learner: it appears in the cv perf table.
  expect_true("plda" %in% fit$learner_cv_perf$learner)
  # ols (gaussian-only) is skipped on a binomial nuisance -> all-NA.
  expect_true(all(is.na(oof[, "ols"])))
})

# ----------------------------------------------------------------------------
#  Test 2 -- plda filtered out on a gaussian nuisance
# ----------------------------------------------------------------------------

test_that("plda is filtered out (all-NA) on a gaussian outcome", {
  d <- .meep_plda_dgp_gaussian(n = 280, p = 4, seed = 302)
  fit <- meep(d$X, d$y, folds = 4L, tune = "none", seed = 302)

  expect_identical(fit$family, "gaussian")
  oof <- fit$oof_matrix$outcome
  expect_true("plda" %in% colnames(oof))
  # plda does not apply to a gaussian nuisance -> all-NA column.
  expect_true(all(is.na(oof[, "plda"])))

  # an inapplicable learner is skipped, not failed: no fold_failures entry
  # naming plda, and plda is absent from learner_cv_perf.
  ff_learners <- vapply(fit$fold_failures,
                        function(z) z$learner, character(1L))
  expect_false("plda" %in% ff_learners)
  expect_false("plda" %in% fit$learner_cv_perf$learner)

  # the ensemble excludes the all-NA plda column entirely.
  expect_equal(unname(fit$ensemble_weights$outcome["plda"]), 0)
})

# ----------------------------------------------------------------------------
#  Test 3 -- plda + weights = error
# ----------------------------------------------------------------------------

test_that("plda in learners with weights errors clearly", {
  d <- .meep_plda_dgp_binary(n = 200, p = 3, seed = 303)
  w <- runif(length(d$y), 0.5, 1.5)

  expect_error(
    meep(d$X, d$y, folds = 4L,
         learners = c("logreg", "plda"), weights = w,
         tune = "none", seed = 303),
    "does not support observation weights")

  # the default learner set also triggers the error when weights are given
  # and the outcome is binomial (plda is in the default set).
  expect_error(
    suppressWarnings(meep(d$X, d$y, folds = 4L, weights = w,
                          tune = "none", seed = 303)),
    "does not support observation weights")
})

# ----------------------------------------------------------------------------
#  Test 4 -- determinism across thread counts
# ----------------------------------------------------------------------------

test_that("meep with plda is byte-identical across thread counts", {
  d <- .meep_plda_dgp_binary(n = 320, p = 4, seed = 304)

  RcppParallel::setThreadOptions(numThreads = 1L)
  fit1 <- suppressWarnings(
    meep(d$X, d$y, folds = 4L,
         learners = c("logreg", "plda"), ensemble = "stack",
         tune = "none", seed = 304))

  RcppParallel::setThreadOptions(numThreads = 4L)
  fit4 <- suppressWarnings(
    meep(d$X, d$y, folds = 4L,
         learners = c("logreg", "plda"), ensemble = "stack",
         tune = "none", seed = 304))
  RcppParallel::setThreadOptions(numThreads = 1L)

  expect_identical(fit1$y_hat_oof, fit4$y_hat_oof)
  expect_identical(fit1$oof_matrix$outcome[, "plda"],
                   fit4$oof_matrix$outcome[, "plda"])
  expect_identical(fit1$ensemble_weights$outcome,
                   fit4$ensemble_weights$outcome)
})

# ----------------------------------------------------------------------------
#  Test 5 -- plda probabilities are sane (correlate with logreg)
# ----------------------------------------------------------------------------

test_that("plda OOF probs are positively correlated with logreg OOF probs", {
  d <- .meep_plda_dgp_binary(n = 400, p = 4, seed = 305)
  fit <- suppressWarnings(
    meep(d$X, d$y, folds = 5L,
         learners = c("logreg", "plda"), tune = "none", seed = 305))

  oof <- fit$oof_matrix$outcome
  pl <- oof[, "plda"]
  lr <- oof[, "logreg"]
  ok <- is.finite(pl) & is.finite(lr)
  expect_gt(stats::cor(pl[ok], lr[ok]), 0.5)
})

# ----------------------------------------------------------------------------
#  Test 6 -- tune = "once" freezes plda hyperparameters (slow)
# ----------------------------------------------------------------------------

test_that("tune = 'once' freezes plda lambda/K and runs", {
  skip_if_quick()
  d <- .meep_plda_dgp_binary(n = 360, p = 4, seed = 306)
  fit <- suppressWarnings(
    meep(d$X, d$y, folds = 4L,
         learners = c("logreg", "plda"), tune = "once", seed = 306))

  fr <- fit$frozen_hyperparams$plda
  expect_true(is.list(fr))
  if (length(fr) > 0L) {
    expect_true("lambda" %in% names(fr))
    expect_false(isTRUE(fr$autotune))
  }
  oof <- fit$oof_matrix$outcome
  expect_true(all(is.finite(oof[, "plda"])))
  expect_true(all(oof[, "plda"] >= 0 & oof[, "plda"] <= 1))
})
