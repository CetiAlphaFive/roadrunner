# test-bgam-meep.R
# T-MEEP-01 through T-MEEP-03: meep() integration with bgam learner
# Derived from test-spec.md Section 3.8
#
# BLOCK: meep() bgam integration is broken when bgam_args contains
# autotune=FALSE (as test-spec specifies). The bgam learner spec in
# .meep_learner_specs() hard-codes autotune=FALSE in args then also
# appends bgam_args, causing:
#   "formal argument 'autotune' matched by multiple actual arguments"
# This causes bgam to fail on every fold when used as sole learner,
# or to log fold_failures when combined with other learners.
#
# API note: test-spec uses meep(Y=Y, D=D, X=X, ...) but meep() uses
# meep(X, y, treatment=D, ...). Y= does not partial-match y= in R.
# Tests use the actual meep API.

library(testthat)

# T-MEEP-01: BLOCKED -- bgam_args=list(autotune=FALSE,...) causes duplicate arg
test_that("T-MEEP-01: meep with bgam gaussian [BLOCKED - autotune dup in bgam_args]", {
  set.seed(50)
  n <- 80
  X <- matrix(rnorm(n * 3), n, 3)
  y_out <- X[, 1] + rnorm(n)
  D <- 0.5 * X[, 2] + rnorm(n)
  # test-spec expects this to succeed; it actually hard-stops
  expect_no_error(
    meep(X = X, y = y_out, treatment = D, learners = "bgam",
         tune = "none", folds = 3, seed = 1,
         bgam_args = list(autotune = FALSE, mstop = 20, nu = 0.1,
                          nthreads = 1))
  )
})

# T-MEEP-02: BLOCKED -- same autotune duplicate arg
test_that("T-MEEP-02: meep with bgam binary outcome [BLOCKED - autotune dup]", {
  set.seed(51)
  n <- 80
  X <- matrix(rnorm(n * 3), n, 3)
  y_out <- rbinom(n, 1, plogis(X[, 1]))
  D <- rnorm(n)
  expect_no_error(
    meep(X = X, y = y_out, treatment = D, learners = "bgam",
         tune = "none", folds = 3, seed = 1,
         bgam_args = list(autotune = FALSE, mstop = 20, nu = 0.1,
                          nthreads = 1))
  )
})

# T-MEEP-03: BLOCKED -- bgam fails every fold; ares succeeds but bgam
# weights should be non-negative; test-spec says all learner weights >= 0.
# With the bgam bug, bgam contributes NA-only OOF columns and gets zero
# ensemble weight. Test shows 6 bgam fold_failures and only ares weights.
test_that("T-MEEP-03: meep with bgam+ares [BLOCKED - bgam fails all folds]", {
  set.seed(50)
  n <- 80
  X <- matrix(rnorm(n * 3), n, 3)
  y_out <- X[, 1] + rnorm(n)
  D <- 0.5 * X[, 2] + rnorm(n)
  fit <- meep(X = X, y = y_out, treatment = D, learners = c("ares", "bgam"),
              tune = "none", folds = 3, seed = 1,
              bgam_args = list(autotune = FALSE, mstop = 20, nu = 0.1,
                               nthreads = 1))
  # test-spec requires bgam to contribute; actual: 6 bgam fold_failures
  expect_equal(length(fit$fold_failures), 0L)
})
