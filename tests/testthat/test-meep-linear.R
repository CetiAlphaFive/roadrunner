# ============================================================================
#  Tests for meep() with the linear learners ols() and logreg()
# ============================================================================
#
# Phase 3 of the ols()/logreg() feature: meep() integration. ols and logreg
# are opt-in learners (the default `learners` stays c("ares", "krls")). They
# carry NO hyperparameters, so `tune` ("once"/"per_fold"/"none") is a no-op:
# every fold is a plain refit. These tests verify registration, OOF honesty,
# stacking, determinism, and graceful degradation for the linear learners.

# ----------------------------------------------------------------------------
#  Helper DGPs
# ----------------------------------------------------------------------------

.meep_lin_dgp <- function(n = 300, p = 4, seed = 1) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  # genuinely linear mean so ols is well specified
  y <- 1 + 0.8 * X[, 1] - 0.5 * X[, 2] + 0.3 * X[, 3] + rnorm(n, sd = 0.3)
  list(X = X, y = y)
}

# ----------------------------------------------------------------------------
#  Test 1 -- learners = c("ols", "logreg"): registration + OOF basics
# ----------------------------------------------------------------------------

test_that("meep registers ols/logreg as learners and produces OOF predictions", {
  set.seed(101)
  n <- 320; p <- 4
  X <- matrix(runif(n * p), n, p)
  D <- rbinom(n, 1, plogis(X[, 1] - 0.5))
  y <- 1 + 0.7 * X[, 1] - 0.4 * X[, 2] + rnorm(n, sd = 0.3)

  # ols = gaussian outcome learner; logreg = binomial treatment learner.
  fit <- suppressWarnings(
    meep(X, y, treatment = D, folds = 4L,
         learners = c("ols", "logreg"), tune = "none", seed = 101))

  expect_s3_class(fit, "meep")
  expect_identical(fit$learners, c("ols", "logreg"))
  expect_length(fit$y_hat_oof, n)
  expect_true(all(is.finite(fit$y_hat_oof)))

  # outcome OOF is gaussian (ols); treatment OOF is a propensity in [0, 1].
  expect_true(all(is.finite(fit$d_hat_oof)))
  expect_true(all(fit$d_hat_oof >= 0 & fit$d_hat_oof <= 1))

  # every nuisance OOF matrix carries both learner columns.
  for (nm in names(fit$oof_matrix))
    expect_identical(colnames(fit$oof_matrix[[nm]]),
                     c("ols", "logreg"))

  # exact residual identity still holds.
  expect_equal(fit$y_resid, y - fit$y_hat_oof, tolerance = 0)
})

# ----------------------------------------------------------------------------
#  Test 2 -- default learners unchanged: ols/logreg are opt-in only
# ----------------------------------------------------------------------------

test_that("ols/logreg are opt-in: the default learners stay c('ares','krls')", {
  d <- .meep_lin_dgp(n = 200, p = 3, seed = 102)
  fit <- meep(d$X, d$y, folds = 4L, tune = "none", seed = 102)
  expect_identical(fit$learners, c("ares", "krls"))
  expect_false(any(c("ols", "logreg") %in% fit$learners))
})

# ----------------------------------------------------------------------------
#  Test 3 -- mixed learners c("ares", "ols")
# ----------------------------------------------------------------------------

test_that("meep mixes ares and ols learners and stacks them", {
  d <- .meep_lin_dgp(n = 360, p = 4, seed = 103)
  fit <- meep(d$X, d$y, folds = 5L,
              learners = c("ares", "ols"), ensemble = "stack",
              tune = "none", seed = 103)

  expect_identical(fit$learners, c("ares", "ols"))
  w <- fit$ensemble_weights$outcome
  expect_named(w, c("ares", "ols"))
  expect_true(all(w >= -1e-8))
  expect_equal(sum(w), 1, tolerance = 1e-6)

  # on a genuinely linear DGP ols should carry meaningful weight.
  expect_gt(w["ols"], 0.2)
})

# ----------------------------------------------------------------------------
#  Test 4 -- OOF honesty / leakage probe with the linear learners
# ----------------------------------------------------------------------------

test_that("ols/logreg OOF columns are honest: held-out y permutation is inert", {
  d <- .meep_lin_dgp(n = 320, p = 4, seed = 104)
  X <- d$X; y <- d$y

  fit0 <- meep(X, y, folds = 4L, learners = c("ols"),
               ensemble = "average", tune = "none", seed = 77)
  folds <- fit0$folds
  te <- which(folds == 1L)

  # permuting y inside the held-out fold cannot move that fold's OOF: the
  # model that predicts fold 1 was trained on folds 2..K.
  y_held <- y
  y_held[te] <- sample(y[te])
  fit_a <- meep(X, y_held, folds = folds, learners = c("ols"),
                ensemble = "average", tune = "none", seed = 77)
  expect_equal(fit_a$oof_matrix$outcome[te, "ols"],
               fit0$oof_matrix$outcome[te, "ols"], tolerance = 1e-10)

  # permuting y in a training fold MUST move fold 1's OOF.
  tr2 <- which(folds == 2L)
  y_tr <- y
  y_tr[tr2] <- sample(y[tr2])
  fit_b <- meep(X, y_tr, folds = folds, learners = c("ols"),
                ensemble = "average", tune = "none", seed = 77)
  expect_false(isTRUE(all.equal(fit_b$oof_matrix$outcome[te, "ols"],
                                fit0$oof_matrix$outcome[te, "ols"],
                                tolerance = 1e-6)))
})

# ----------------------------------------------------------------------------
#  Test 5 -- tune is a no-op for ols/logreg (no hyperparameters)
# ----------------------------------------------------------------------------

test_that("tune = 'once' froze nothing for ols/logreg and matches tune = 'none'", {
  set.seed(105)
  n <- 320; p <- 4
  X <- matrix(runif(n * p), n, p)
  D <- rbinom(n, 1, plogis(X[, 1] - 0.5))
  y <- 1 + 0.7 * X[, 1] - 0.4 * X[, 2] + rnorm(n, sd = 0.3)

  fit_once <- suppressWarnings(
    meep(X, y, treatment = D, folds = 4L,
         learners = c("ols", "logreg"), tune = "once", seed = 105))
  fit_none <- suppressWarnings(
    meep(X, y, treatment = D, folds = 4L,
         learners = c("ols", "logreg"), tune = "none", seed = 105))

  # ols/logreg have no hyperparameters: freezing must produce empty lists.
  expect_true(all(vapply(fit_once$frozen_hyperparams,
                         function(z) length(z) == 0L, logical(1L))))

  # with no hyperparameters to freeze, tune = "once" and tune = "none"
  # produce byte-identical OOF predictions.
  expect_equal(fit_once$y_hat_oof, fit_none$y_hat_oof, tolerance = 0)
  expect_equal(fit_once$d_hat_oof, fit_none$d_hat_oof, tolerance = 0)

  # per_fold is likewise a no-op for these learners.
  fit_pf <- suppressWarnings(
    meep(X, y, treatment = D, folds = 4L,
         learners = c("ols", "logreg"), tune = "per_fold", seed = 105))
  expect_equal(fit_pf$y_hat_oof, fit_none$y_hat_oof, tolerance = 0)
})

# ----------------------------------------------------------------------------
#  Test 6 -- determinism: nthreads = 1 byte-identical to nthreads = 4
# ----------------------------------------------------------------------------

test_that("meep with ols/logreg learners is byte-identical across thread counts", {
  set.seed(106)
  n <- 360; p <- 4
  X <- matrix(runif(n * p), n, p)
  D <- rbinom(n, 1, plogis(X[, 1] - 0.5))
  y <- 1 + 0.8 * X[, 1] - 0.5 * X[, 2] + rnorm(n, sd = 0.3)

  RcppParallel::setThreadOptions(numThreads = 1L)
  fit1 <- suppressWarnings(
    meep(X, y, treatment = D, folds = 4L,
         learners = c("ols", "logreg"), ensemble = "stack",
         tune = "none", seed = 106))

  RcppParallel::setThreadOptions(numThreads = 4L)
  fit4 <- suppressWarnings(
    meep(X, y, treatment = D, folds = 4L,
         learners = c("ols", "logreg"), ensemble = "stack",
         tune = "none", seed = 106))
  RcppParallel::setThreadOptions(numThreads = 1L)

  expect_identical(fit1$y_hat_oof, fit4$y_hat_oof)
  expect_identical(fit1$d_hat_oof, fit4$d_hat_oof)
  expect_identical(fit1$ensemble_weights$outcome,
                   fit4$ensemble_weights$outcome)
})

# ----------------------------------------------------------------------------
#  Test 7 -- predict.meep() works with the linear learners
# ----------------------------------------------------------------------------

test_that("predict.meep refits ols/logreg on full data for new rows", {
  set.seed(107)
  n <- 320; p <- 4
  X <- matrix(runif(n * p), n, p)
  D <- rbinom(n, 1, plogis(X[, 1] - 0.5))
  y <- 1 + 0.7 * X[, 1] - 0.4 * X[, 2] + rnorm(n, sd = 0.3)
  fit <- suppressWarnings(
    meep(X, y, treatment = D, folds = 4L,
         learners = c("ols", "logreg"), tune = "none", seed = 107))

  Xnew <- matrix(runif(20 * p), 20, p)
  py <- predict(fit, Xnew, nuisance = "outcome")
  pd <- predict(fit, Xnew, nuisance = "treatment")
  expect_length(py, 20L)
  expect_true(all(is.finite(py)))
  expect_true(all(pd >= 0 & pd <= 1))
})

# ----------------------------------------------------------------------------
#  Test 8 -- graceful degradation: a learner that errors on a fold -> NA col
# ----------------------------------------------------------------------------

# small ares adapter mirroring the package's built-in ares spec, used by the
# graceful-degradation test below (a list of learners bypasses the registry).
.meep_test_ares_spec <- function() {
  list(
    fit = function(X, y, family, weights, hp) ares(x = X, y = y),
    predict = function(model, newX, family)
      as.numeric(predict(model, newdata = newX, type = "response"))
  )
}

test_that("a failing ols fold yields an NA column and renormalized weights", {
  d <- .meep_lin_dgp(n = 240, p = 3, seed = 108)

  # adapter wrapping ols() that throws on fold 2 to exercise the
  # graceful-degradation path; ares is the surviving learner.
  flaky <- local({
    function() {
      list(
        fit = function(X, y, family, weights, hp) {
          if (abs(nrow(X) - 180L) <= 2L)  # ~3 of 4 folds -> ~180 train rows
            stop("synthetic ols fold failure")
          ols(x = X, y = y)
        },
        predict = function(model, newX, family)
          as.numeric(predict(model, newdata = newX, type = "response"))
      )
    }
  })()

  fit <- meep(d$X, d$y, folds = 4L,
              learners = list(ares = .meep_test_ares_spec(),
                              ols_flaky = flaky),
              ensemble = "average", tune = "none", seed = 108)

  # at least one fold's ols column is NA; ares still covers every row.
  oof <- fit$oof_matrix$outcome
  expect_true(anyNA(oof[, "ols_flaky"]))
  expect_true(all(is.finite(oof[, "ares"])))
  expect_true(all(is.finite(fit$y_hat_oof)))
  expect_gt(length(fit$fold_failures), 0L)
})
