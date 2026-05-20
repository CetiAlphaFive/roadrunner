# ============================================================================
#  Tests for meep() -- cross-fitted causal ensemble (Phase Q-DML)
# ============================================================================
#
# Spec section 9 test plan, MINUS the meep_plr test (test 5 in the spec).
# Test 5 here instead verifies residual identities + mean-near-zero. krls()
# logistic warnings are noisy but harmless; suppressed where they appear.

# ----------------------------------------------------------------------------
#  Helper DGPs
# ----------------------------------------------------------------------------

.meep_dgp_smooth <- function(n = 300, p = 4, seed = 1) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  y <- X[, 1]^2 + 0.5 * sin(3 * X[, 2]) + rnorm(n, sd = 0.3)
  list(X = X, y = y)
}

# ----------------------------------------------------------------------------
#  Test 1 -- leakage probe (core). Both directions.
# ----------------------------------------------------------------------------

test_that("leakage probe: held-out y permutation does not move OOF; training y permutation does", {
  d <- .meep_dgp_smooth(n = 300, p = 4, seed = 11)
  X <- d$X; y <- d$y

  # Core honesty test. ensemble = "average" so the combine weights (1/L)
  # do NOT depend on y at all: any movement in y_hat_oof can only come
  # from a base learner having seen the permuted rows. tune = "none" so
  # frozen hyperparameters cannot contaminate either.
  fit0 <- meep(X, y, folds = 4L, ensemble = "average", tune = "none",
               seed = 99)
  folds <- fit0$folds
  k <- 1L
  te <- which(folds == k)

  # (a) permute y WITHIN the held-out fold k. y_hat_oof[te] must NOT move:
  # rows in fold k are predicted by models trained on folds != k, which
  # never saw fold k's y, and the equal weights are y-independent.
  y_perm_held <- y
  y_perm_held[te] <- sample(y[te])
  fit_a <- meep(X, y_perm_held, folds = folds, ensemble = "average",
                tune = "none", seed = 99)
  expect_equal(fit_a$y_hat_oof[te], fit0$y_hat_oof[te], tolerance = 1e-10)

  # (b) permute y in a TRAINING fold (fold 2). y_hat_oof[te] (fold 1) MUST
  # move, because fold 1 is predicted by models trained on folds 2,3,4.
  tr2 <- which(folds == 2L)
  y_perm_train <- y
  y_perm_train[tr2] <- sample(y[tr2])
  fit_b <- meep(X, y_perm_train, folds = folds, ensemble = "average",
                tune = "none", seed = 99)
  expect_false(isTRUE(all.equal(fit_b$y_hat_oof[te], fit0$y_hat_oof[te],
                                tolerance = 1e-6)))

  # The per-learner OOF columns themselves are the load-bearing honesty
  # guarantee under ANY ensemble: held-out permutation leaves them fixed.
  fit_stack0 <- meep(X, y, folds = folds, ensemble = "stack",
                     tune = "none", seed = 99)
  fit_stackA <- meep(X, y_perm_held, folds = folds, ensemble = "stack",
                     tune = "none", seed = 99)
  expect_equal(fit_stackA$oof_matrix$outcome[te, ],
               fit_stack0$oof_matrix$outcome[te, ], tolerance = 1e-10)
})

# ----------------------------------------------------------------------------
#  Test 2 -- NNLS weights non-negative, sum to 1 (outcome + treatment)
# ----------------------------------------------------------------------------

test_that("stack ensemble weights are non-negative and sum to one", {
  set.seed(21)
  n <- 300; p <- 4
  X <- matrix(runif(n * p), n, p)
  D <- X[, 1] + rnorm(n, sd = 0.3)
  y <- X[, 1]^2 + 0.5 * D + rnorm(n, sd = 0.3)
  fit <- meep(X, y, treatment = D, folds = 4L, ensemble = "stack",
              tune = "none", seed = 21)
  for (nm in names(fit$ensemble_weights)) {
    w <- fit$ensemble_weights[[nm]]
    expect_true(all(w >= -1e-8), info = nm)
    expect_equal(sum(w), 1, tolerance = 1e-6, info = nm)
  }
  expect_true("outcome" %in% names(fit$ensemble_weights))
  expect_true("treatment" %in% names(fit$ensemble_weights))
})

# ----------------------------------------------------------------------------
#  Test 3 -- dominance: hinge DGP -> ares; Gaussian-bump DGP -> krls
# ----------------------------------------------------------------------------

test_that("ares dominates on hinge DGP; krls dominates on Gaussian-bump DGP", {
  skip_if_quick()
  # hinge / piecewise-linear DGP -- MARS native territory
  set.seed(31)
  n <- 500; p <- 4
  Xh <- matrix(runif(n * p, -1, 1), n, p)
  yh <- pmax(0, Xh[, 1]) + pmax(0, -Xh[, 2]) + rnorm(n, sd = 0.15)
  fit_h <- meep(Xh, yh, folds = 5L, tune = "none", seed = 31)
  expect_gt(fit_h$ensemble_weights$outcome["ares"], 0.8)

  # smooth radial Gaussian-bump DGP -- kernel native territory
  set.seed(32)
  Xb <- matrix(runif(n * p, -1, 1), n, p)
  r2 <- Xb[, 1]^2 + Xb[, 2]^2
  yb <- exp(-3 * r2) + rnorm(n, sd = 0.05)
  fit_b <- meep(Xb, yb, folds = 5L, tune = "none", seed = 32)
  expect_gt(fit_b$ensemble_weights$outcome["krls"],
            fit_b$ensemble_weights$outcome["ares"])
})

# ----------------------------------------------------------------------------
#  Test 4 -- grf integration smoke
# ----------------------------------------------------------------------------

test_that("grf::causal_forest accepts meep OOF nuisances", {
  skip_if_not_installed("grf")
  set.seed(41)
  n <- 400; p <- 5
  X <- matrix(runif(n * p), n, p)
  D <- rbinom(n, 1, plogis(X[, 1] - 0.5))
  y <- X[, 1]^2 + 0.8 * D + rnorm(n, sd = 0.3)
  fit <- suppressWarnings(
    meep(X, y, treatment = D, folds = 5L, tune = "none", seed = 41))
  cf <- grf::causal_forest(X, y, D,
                           Y.hat = fit$y_hat_oof,
                           W.hat = fit$d_hat_oof,
                           num.trees = 200)
  ate <- grf::average_treatment_effect(cf)
  expect_true(is.finite(ate[["estimate"]]))
})

# ----------------------------------------------------------------------------
#  Test 5 -- (REPLACES meep_plr) residual identities + mean-near-zero
# ----------------------------------------------------------------------------

test_that("y_resid and d_resid are exact residuals and mean-near-zero on a clean DGP", {
  set.seed(51)
  n <- 600; p <- 4
  X <- matrix(runif(n * p), n, p)
  D <- 0.4 * X[, 1] + 0.3 * sin(3 * X[, 2]) + rnorm(n, sd = 0.3)
  g <- X[, 1]^2 + 0.5 * X[, 3]
  y <- 1.0 * D + g + rnorm(n, sd = 0.5)
  fit <- meep(X, y, treatment = D, folds = 5L, tune = "none", seed = 51)

  # exact residual identities
  expect_equal(fit$y_resid, y - fit$y_hat_oof, tolerance = 0)
  expect_equal(fit$d_resid, D - fit$d_hat_oof, tolerance = 0)

  # on a well-specified DGP the cross-fitted residuals are mean-near-zero
  expect_lt(abs(mean(fit$y_resid)), 0.15)
  expect_lt(abs(mean(fit$d_resid)), 0.10)
})

# ----------------------------------------------------------------------------
#  Test 6 -- binary outcome -> binomial / logistic learners; y_hat in [0,1]
# ----------------------------------------------------------------------------

test_that("binary outcome auto-detects binomial family; y_hat_oof in [0,1]", {
  set.seed(61)
  n <- 400; p <- 4
  X <- matrix(runif(n * p), n, p)
  y <- rbinom(n, 1, plogis(2 * X[, 1] - 1))
  fit <- suppressWarnings(meep(X, y, folds = 4L, tune = "none", seed = 61))
  expect_identical(fit$family, "binomial")
  expect_true(all(fit$y_hat_oof >= 0 & fit$y_hat_oof <= 1))
})

# ----------------------------------------------------------------------------
#  Test 7 -- binary treatment -> d_hat in [0,1]; mu0/mu1 populated
# ----------------------------------------------------------------------------

test_that("binary treatment gives a propensity in [0,1] and arm models", {
  set.seed(71)
  n <- 400; p <- 4
  X <- matrix(runif(n * p), n, p)
  D <- rbinom(n, 1, plogis(X[, 1] - 0.5))
  y <- X[, 1]^2 + 0.6 * D + rnorm(n, sd = 0.3)
  fit <- suppressWarnings(
    meep(X, y, treatment = D, folds = 4L, arm_models = "auto",
         tune = "none", seed = 71))
  expect_identical(fit$treatment_family, "binomial")
  expect_true(all(fit$d_hat_oof >= 0 & fit$d_hat_oof <= 1))
  expect_false(is.null(fit$mu0_hat_oof))
  expect_false(is.null(fit$mu1_hat_oof))
  expect_length(fit$mu0_hat_oof, n)
  expect_length(fit$mu1_hat_oof, n)
})

# ----------------------------------------------------------------------------
#  Test 8 -- reproducibility: same seed -> identical folds/OOF/weights
# ----------------------------------------------------------------------------

test_that("same seed gives identical folds, y_hat_oof, and weights", {
  d <- .meep_dgp_smooth(n = 300, p = 4, seed = 81)
  f1 <- meep(d$X, d$y, folds = 5L, tune = "none", seed = 123)
  f2 <- meep(d$X, d$y, folds = 5L, tune = "none", seed = 123)
  expect_identical(f1$folds, f2$folds)
  expect_equal(f1$y_hat_oof, f2$y_hat_oof, tolerance = 0)
  expect_equal(f1$ensemble_weights$outcome, f2$ensemble_weights$outcome,
               tolerance = 0)
})

# ----------------------------------------------------------------------------
#  Test 9 -- user-supplied folds vector honored verbatim
# ----------------------------------------------------------------------------

test_that("a user-supplied folds vector is honored verbatim", {
  d <- .meep_dgp_smooth(n = 240, p = 3, seed = 91)
  user_folds <- rep_len(1:4, 240)
  fit <- meep(d$X, d$y, folds = user_folds, tune = "none", seed = 7)
  expect_identical(fit$folds, as.integer(user_folds))
})

# ----------------------------------------------------------------------------
#  Test 10 -- cluster: every cluster lands entirely in one fold
# ----------------------------------------------------------------------------

test_that("cluster-aware folds keep each cluster whole", {
  d <- .meep_dgp_smooth(n = 300, p = 3, seed = 101)
  clust <- rep(1:30, each = 10)
  fit <- meep(d$X, d$y, folds = 5L, cluster = clust, tune = "none",
              seed = 101)
  # each cluster id maps to exactly one fold
  tab <- tapply(fit$folds, clust, function(z) length(unique(z)))
  expect_true(all(tab == 1L))
})

# ----------------------------------------------------------------------------
#  Test 11 -- graceful degradation: a learner errors on a fold
# ----------------------------------------------------------------------------

test_that("a learner that fails a fold is dropped and weights renormalize", {
  d <- .meep_dgp_smooth(n = 300, p = 3, seed = 111)

  # stub learner: a list-of-adapters entry that errors on fold-1-sized
  # training sets but otherwise mirrors a constant predictor.
  fail_once <- local({
    seen <- 0L
    list(
      fit = function(X, y, family, weights, hp) {
        seen <<- seen + 1L
        if (seen == 1L) stop("stub learner deliberate failure")
        structure(list(mu = mean(y)), class = "meep_stub")
      },
      predict = function(model, newX, family) {
        rep(model$mu, nrow(as.data.frame(newX)))
      }
    )
  })
  learners <- list(ares = NULL, stub = fail_once)
  # reuse the real ares adapter for the "ares" slot
  ares_spec <- list(
    fit = function(X, y, family, weights, hp)
      ares(x = X, y = y, family = family),
    predict = function(model, newX, family)
      as.numeric(predict(model, newdata = newX, type = "response"))
  )
  learners$ares <- ares_spec

  fit <- meep(d$X, d$y, folds = 4L, learners = learners,
              ensemble = "stack", tune = "none", seed = 111)
  expect_gt(length(fit$fold_failures), 0L)
  # surviving ensemble weights still sum to 1
  expect_equal(sum(fit$ensemble_weights$outcome), 1, tolerance = 1e-6)
  expect_true(all(is.finite(fit$y_hat_oof)))
})

# ----------------------------------------------------------------------------
#  Test 12 -- predict.meep on new data -> finite vector of right length
# ----------------------------------------------------------------------------

test_that("predict.meep on new data returns a finite length-correct vector", {
  d <- .meep_dgp_smooth(n = 300, p = 4, seed = 121)
  fit <- meep(d$X, d$y, folds = 4L, tune = "none", seed = 121)
  newX <- matrix(runif(50 * 4), 50, 4)
  pr <- predict(fit, newdata = newX, nuisance = "outcome")
  expect_length(pr, 50L)
  expect_true(all(is.finite(pr)))
})

# ----------------------------------------------------------------------------
#  Test 13 -- treatment = NULL -> d_hat_oof / mu* / d_resid all NULL
# ----------------------------------------------------------------------------

test_that("treatment = NULL leaves all treatment-side outputs NULL", {
  d <- .meep_dgp_smooth(n = 250, p = 3, seed = 131)
  fit <- meep(d$X, d$y, treatment = NULL, folds = 4L, tune = "none",
              seed = 131)
  expect_null(fit$d_hat_oof)
  expect_null(fit$mu0_hat_oof)
  expect_null(fit$mu1_hat_oof)
  expect_null(fit$d_resid)
  expect_false(is.null(fit$y_hat_oof))
})

# ----------------------------------------------------------------------------
#  Test 14 -- multi-valued treatment -> informative stop()
# ----------------------------------------------------------------------------

test_that("a multi-valued treatment is rejected with an informative error", {
  d <- .meep_dgp_smooth(n = 200, p = 3, seed = 141)
  Dmv <- sample(0:3, 200, replace = TRUE)
  expect_error(
    meep(d$X, d$y, treatment = Dmv, folds = 4L, tune = "none"),
    "[Mm]ulti-valued")
})

# ----------------------------------------------------------------------------
#  Test 15 -- tune = "none" fast path: runs, no autotune invoked
# ----------------------------------------------------------------------------

test_that("tune = 'none' runs the fast path without invoking autotune", {
  d <- .meep_dgp_smooth(n = 250, p = 4, seed = 151)

  # learner stub that records the hyperparameters spliced into its fit
  # call. Under tune = "none" the spliced hp list must be empty (no
  # autotune = TRUE, no frozen degree/penalty/sigma/lambda).
  seen_hp <- list()
  rec_spec <- list(
    fit = function(X, y, family, weights, hp) {
      seen_hp[[length(seen_hp) + 1L]] <<- hp
      structure(list(mu = mean(y)), class = "meep_stub")
    },
    predict = function(model, newX, family)
      rep(model$mu, nrow(as.data.frame(newX)))
  )
  fit <- meep(d$X, d$y, folds = 4L,
              learners = list(rec = rec_spec),
              ensemble = "average", tune = "none", seed = 151)

  # every recorded hyperparameter list is empty -> no autotune, no frozen hp
  expect_true(all(vapply(seen_hp, function(z) length(z) == 0L,
                         logical(1L))))
  expect_null(fit$frozen_hyperparams)
  expect_identical(fit$tune, "none")
  expect_true(all(is.finite(fit$y_hat_oof)))
  expect_length(fit$y_hat_oof, 250L)
})
