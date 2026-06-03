# Tests for meep() opt-in external-package learners: "forest" (ranger) and
# "BART" (dbarts). These packages are Suggests-only; tests are skipped when
# the package is not installed. Token-normalization and availability-check
# logic (case-insensitivity, unknown-token rejection) is tested via the
# package-internal helper and does NOT require ranger/dbarts to be present.

# ---- normalization / availability helper (no external pkgs needed) --------

test_that(".meep_normalize_extra_learners canonicalizes case-insensitively", {
  f <- roadrunner:::.meep_normalize_extra_learners
  expect_equal(f("forest"), "forest")
  expect_equal(f("Forest"), "forest")
  expect_equal(f("FOREST"), "forest")
  expect_equal(f("BART"), "bart")
  expect_equal(f("bart"), "bart")
  expect_equal(f("Bart"), "bart")
})

test_that(".meep_normalize_extra_learners dedupes and accepts vectors", {
  f <- roadrunner:::.meep_normalize_extra_learners
  expect_equal(f(c("forest", "BART")), c("forest", "bart"))
  expect_equal(f(c("BART", "bart")), "bart")
  expect_equal(f(c("Forest", "forest")), "forest")
})

test_that(".meep_normalize_extra_learners returns character(0) for NULL/empty", {
  f <- roadrunner:::.meep_normalize_extra_learners
  expect_equal(f(NULL), character(0))
  expect_equal(f(character(0)), character(0))
})

test_that(".meep_normalize_extra_learners rejects unknown tokens", {
  f <- roadrunner:::.meep_normalize_extra_learners
  expect_error(f("xgboost"), "forest")
  expect_error(f("xgboost"), "BART")
  expect_error(f(c("forest", "lightgbm")), "lightgbm")
})

test_that("invalid extra.learners token errors from meep() listing valid values", {
  set.seed(1)
  n <- 40; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1] + rnorm(n)
  expect_error(
    meep(X, Y, folds = 2L, learners = "ols", extra.learners = "xgboost"),
    "forest")
})

# ---- forest learner: gaussian outcome -------------------------------------

test_that("extra.learners='forest' adds a finite OOF column on a gaussian meep", {
  skip_if_not_installed("ranger")
  set.seed(2)
  n <- 120; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1] + 0.5 * X[, 2]^2 + rnorm(n)
  fit <- meep(X, Y, folds = 3L, learners = "ols",
              extra.learners = "forest", seed = 7,
              forest_args = list(num.trees = 100))
  expect_true("forest" %in% fit$learners)
  oof <- fit$oof_matrix$outcome[, "forest"]
  expect_true(all(is.finite(oof)))
  expect_true("forest" %in% names(fit$ensemble_weights$outcome))
})

# ---- BART learner: binary outcome + treatment -----------------------------

test_that("extra.learners='BART' adds in-range probabilities on a binary meep", {
  skip_if_not_installed("dbarts")
  set.seed(3)
  n <- 120; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  lin <- X[, 1] - 0.5 * X[, 2]
  Y <- rbinom(n, 1, plogis(lin))
  D <- rbinom(n, 1, plogis(0.5 * X[, 3]))
  fit <- meep(X, Y, treatment = D, folds = 3L, learners = "logreg",
              extra.learners = "BART", seed = 11,
              bart_args = list(nskip = 50, ndpost = 100))
  expect_true("bart" %in% fit$learners)
  oof <- fit$oof_matrix$outcome[, "bart"]
  ok <- is.finite(oof)
  expect_true(any(ok))
  expect_true(all(oof[ok] > 0 & oof[ok] < 1))
})

# ---- both together --------------------------------------------------------

test_that("extra.learners=c('forest','BART') adds both learners", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("dbarts")
  set.seed(4)
  n <- 120; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1] + rnorm(n)
  fit <- meep(X, Y, folds = 3L, learners = "ols",
              extra.learners = c("forest", "BART"), seed = 5,
              forest_args = list(num.trees = 100),
              bart_args = list(nskip = 50, ndpost = 100))
  expect_true(all(c("forest", "bart") %in% fit$learners))
})

# ---- mix of built-in and extra learners -----------------------------------

test_that("mixing built-in + extra learners yields a stack over all names", {
  skip_if_not_installed("ranger")
  set.seed(6)
  n <- 120; p <- 4
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1] + 0.3 * X[, 3] + rnorm(n)
  fit <- meep(X, Y, folds = 3L,
              learners = c("ares", "ols"),
              extra.learners = "forest", seed = 9,
              forest_args = list(num.trees = 100))
  expect_true(all(c("ares", "ols", "forest") %in% fit$learners))
  w <- fit$ensemble_weights$outcome
  expect_true(all(c("ares", "ols", "forest") %in% names(w)))
  expect_equal(sum(w), 1, tolerance = 1e-8)
})

# ---- plot smoke test with extra learners present --------------------------

test_that("plot.meep works when extra learners are present", {
  skip_if_not_installed("ranger")
  set.seed(8)
  n <- 120; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  D <- rbinom(n, 1, plogis(X[, 1]))
  Y <- X[, 1] + 0.7 * D + rnorm(n)
  fit <- meep(X, Y, treatment = D, folds = 3L,
              learners = c("ares", "ols"),
              extra.learners = "forest", seed = 3,
              forest_args = list(num.trees = 100))
  pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_silent(plot(fit))
})

# ---- predict.meep works with an extra learner -----------------------------

test_that("predict.meep produces finite predictions with extra learners", {
  skip_if_not_installed("ranger")
  set.seed(10)
  n <- 120; p <- 3
  X <- matrix(rnorm(n * p), n, p)
  Y <- X[, 1] + rnorm(n)
  fit <- meep(X, Y, folds = 3L, learners = "ols",
              extra.learners = "forest", seed = 4,
              forest_args = list(num.trees = 100))
  newX <- matrix(rnorm(20 * p), 20, p)
  pr <- predict(fit, newdata = newX, nuisance = "outcome")
  expect_length(pr, 20L)
  expect_true(all(is.finite(pr)))
})
