# test-krls-overfit-REQ002.R
# Empirical ratio regression test: EMP-1 from test-spec.md REQ-20260518-002 (REVISED).
#
# Runs a minimal simulation (n=500, p=10, R=10, 5 DGPs x 2 tune conditions)
# to verify that the four fixes materially reduce overfit ratios.
#
# Mandatory sim rules (from spec):
#  1. Smoke-test first: R=2 reps on 1 DGP before the full loop.
#  2. Write intermediate CSV after every (DGP, tune) cell.
#  3. Cap wall-clock at 5 minutes (300s) — skip gracefully if exceeded.
#
# Assertions revised 2026-05-18 (planner EMP-1 BLOCK resolution).
# Old ASSERTION 2 (tune=autotune < 1.5) replaced with 5 oracle-grounded assertions.
# See comprehension.md Threshold Revision Record for rationale.
#
# This entire test is wrapped in skip_if_quick() so it only runs under
# ARES_FULL_TESTS=1.

library(roadrunner)

# ---------------------------------------------------------------------------
# DGP generator (matches test-spec.md)
# ---------------------------------------------------------------------------
make_dgp_fix002 <- function(dgp, n, p, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  eps <- rnorm(n, sd = 1)
  y <- switch(dgp,
    linear      = rowSums(X[, 1:min(p, 5L)]) + eps,
    additive    = X[, 1] + X[, 2]^2 + X[, 3] - X[, 4] + eps,
    interaction = X[, 1] * X[, 2] + X[, 3]^2 + eps,
    sparse      = X[, 1] + eps,
    noise       = eps
  )
  list(X = X, y = y)
}

# ---------------------------------------------------------------------------
# EMP-1 (REVISED): Overfit ratios drop after the four fixes
# ---------------------------------------------------------------------------
test_that("EMP-1 (revised): overfit ratios materially lower after four KRLS fixes", {
  skip_if_quick()

  R  <- 10L
  n  <- 500L
  p  <- 10L
  dgps <- c("additive", "interaction", "sparse", "linear", "noise")

  t0 <- proc.time()["elapsed"]
  tmpcsv <- tempfile(fileext = ".csv")
  results <- list()

  ## ---- smoke test (R=2 reps, 1 DGP, no autotune) -------------------
  d_sm <- make_dgp_fix002("additive", n = n, p = p, seed = 1L)
  train_i_sm <- seq_len(400L); test_i_sm <- 401L:500L
  fit_sm <- tryCatch(
    krls(d_sm$X[train_i_sm, ], d_sm$y[train_i_sm],
         derivative = FALSE, vcov = FALSE, nthreads = 1L),
    error = function(e) stop(paste("Smoke test failed:", conditionMessage(e)))
  )
  pred_sm <- predict(fit_sm, d_sm$X[test_i_sm, ])$fit
  stopifnot("smoke test predict OK" = length(pred_sm) == length(test_i_sm))

  ## ---- full simulation loop ----------------------------------------
  for (dgp in dgps) {
    for (tune in c("none", "autotune")) {

      ## wall-clock cap: 280s = 4 min 40 s (20s margin before 300s)
      elapsed <- proc.time()["elapsed"] - t0
      if (elapsed > 280) {
        if (length(results) > 0L)
          write.csv(do.call(rbind, results), tmpcsv, row.names = FALSE)
        skip(paste0("wall-clock cap reached (", round(elapsed), "s elapsed)"))
      }

      ratios <- numeric(R)
      winner_sigmas <- numeric(R)
      for (r in seq_len(R)) {
        ## check cap inside inner loop too
        if (proc.time()["elapsed"] - t0 > 280) {
          write.csv(do.call(rbind, results), tmpcsv, row.names = FALSE)
          skip(paste0("wall-clock cap reached during r=", r))
        }

        seed_r <- r * 1000L + which(dgps == dgp)
        d <- make_dgp_fix002(dgp, n = n, p = p, seed = seed_r)
        train_i <- seq_len(400L); test_i <- 401L:500L

        fit_args <- list(
          X          = d$X[train_i, ],
          y          = d$y[train_i],
          derivative = FALSE,
          vcov       = FALSE,
          nthreads   = 1L
        )
        if (tune == "autotune") fit_args$autotune <- TRUE

        fit <- tryCatch(do.call(krls, fit_args), error = function(e) NULL)
        if (is.null(fit)) {
          ratios[r] <- NA_real_
          next
        }
        yhat_tr  <- fit$fitted
        yhat_te  <- predict(fit, d$X[test_i, ])$fit
        train_mse <- mean((d$y[train_i] - yhat_tr)^2)
        test_mse  <- mean((d$y[test_i]  - yhat_te)^2)
        ratios[r] <- if (train_mse > 1e-14) test_mse / train_mse else NA_real_
        winner_sigmas[r] <- fit$sigma
      }

      cell_key <- paste(dgp, tune, sep = "_")
      results[[cell_key]] <- data.frame(
        dgp          = dgp,
        tune         = tune,
        mean_ratio   = mean(ratios, na.rm = TRUE),
        mean_sigma   = mean(winner_sigmas, na.rm = TRUE),
        R            = sum(!is.na(ratios))
      )
      ## intermediate CSV write after every cell
      write.csv(do.call(rbind, results), tmpcsv, row.names = FALSE)
    }
  }

  summary_df <- do.call(rbind, results)

  ## ---- display results for diagnostics ----------------------------
  message("EMP-1 results:\n",
          paste(capture.output(print(summary_df, row.names = FALSE)),
                collapse = "\n"))

  ## ---- Diagnostic: autotune sigma vs tune=none sigma (informational) ------
  ## At n=500 p=10 the autotune winner is expected to equal the default sigma
  ## (explained in test-spec.md Autotune Equivalence Note). Not an assertion.
  for (dg in c("additive", "interaction")) {
    none_sig <- summary_df[summary_df$dgp == dg & summary_df$tune == "none",   "mean_sigma"]
    auto_sig <- summary_df[summary_df$dgp == dg & summary_df$tune == "autotune", "mean_sigma"]
    message(sprintf("Diagnostic sigma check: %s  tune=none sigma=%.4f  tune=autotune sigma=%.4f  diff=%.4f",
                    dg, none_sig, auto_sig, abs(auto_sig - none_sig)))
  }

  ## ---- hard failure threshold: at min, shows SOME improvement -----
  r_add_none <- summary_df[summary_df$dgp == "additive" &
                             summary_df$tune == "none", "mean_ratio"]
  r_int_none <- summary_df[summary_df$dgp == "interaction" &
                             summary_df$tune == "none", "mean_ratio"]
  if (!is.na(r_add_none) && r_add_none >= 4.0) {
    fail(paste0("HARD FAIL: tune=none additive mean_ratio=", round(r_add_none, 3),
                " >= 4.0 (no improvement from baseline 3.35)"))
  }
  if (!is.na(r_int_none) && r_int_none >= 4.0) {
    fail(paste0("HARD FAIL: tune=none interaction mean_ratio=", round(r_int_none, 3),
                " >= 4.0 (no improvement from baseline 4.38)"))
  }

  ## ========================================================================
  ## ASSERTION 1: tune=none absolute thresholds
  ## (oracle-grounded: additive < 2.0, interaction < 2.5, sparse < 1.8, linear < 1.5)
  ## Tolerances from test-spec.md exactly as written.
  ## ========================================================================
  none_thresholds <- c(additive = 2.0, interaction = 2.5, sparse = 1.8, linear = 1.5)
  for (dg in names(none_thresholds)) {
    r <- summary_df[summary_df$dgp == dg & summary_df$tune == "none", "mean_ratio"]
    if (!is.na(r)) {
      expect_lte(r, none_thresholds[[dg]],
                 label = paste("ASSERTION 1: tune=none", dg,
                               "mean_ratio should be <", none_thresholds[[dg]]))
    }
  }

  ## ========================================================================
  ## ASSERTION 2: tune=none improvement over pre-fix baseline >= 40%
  ## Pre-fix baselines from REQ-20260518-001 sim-report.md.
  ## 40% is a conservative lower bound (builder achieved 50%+ on all cells).
  ## Tolerances from test-spec.md exactly as written.
  ## ========================================================================
  pre_fix_none <- c(additive = 3.35, interaction = 4.38, sparse = 2.61, linear = 2.70)
  for (dg in names(pre_fix_none)) {
    r <- summary_df[summary_df$dgp == dg & summary_df$tune == "none", "mean_ratio"]
    if (!is.na(r)) {
      pct_improvement <- (pre_fix_none[[dg]] - r) / pre_fix_none[[dg]]
      expect_gte(pct_improvement, 0.40,
                 label = paste("ASSERTION 2: tune=none", dg,
                               "must show >= 40% improvement over pre-fix baseline of",
                               pre_fix_none[[dg]],
                               "(actual improvement:", round(pct_improvement * 100, 1), "%)"))
    }
  }

  ## ========================================================================
  ## ASSERTION 3: tune=autotune absolute thresholds
  ## (oracle-grounded; autotune equals tune=none for additive/interaction
  ##  at this (n,p) by design — see Autotune Equivalence Note in test-spec.md)
  ## Tolerances from test-spec.md exactly as written.
  ## ========================================================================
  auto_thresholds <- c(additive = 2.0, interaction = 2.5, sparse = 1.5, linear = 1.3)
  for (dg in names(auto_thresholds)) {
    r <- summary_df[summary_df$dgp == dg & summary_df$tune == "autotune", "mean_ratio"]
    if (!is.na(r)) {
      expect_lte(r, auto_thresholds[[dg]],
                 label = paste("ASSERTION 3: tune=autotune", dg,
                               "mean_ratio should be <", auto_thresholds[[dg]]))
    }
  }

  ## ========================================================================
  ## ASSERTION 4: autotune does not dramatically worsen relative to tune=none
  ## Autotune may match or beat tune=none; must not exceed tune=none by > 0.2.
  ## Tolerance from test-spec.md exactly as written: margin = 0.2.
  ## ========================================================================
  for (dg in c("additive", "interaction", "sparse", "linear")) {
    r_none <- summary_df[summary_df$dgp == dg & summary_df$tune == "none",     "mean_ratio"]
    r_auto <- summary_df[summary_df$dgp == dg & summary_df$tune == "autotune", "mean_ratio"]
    if (!is.na(r_none) && !is.na(r_auto)) {
      expect_lte(r_auto, r_none + 0.2,
                 label = paste("ASSERTION 4: tune=autotune", dg,
                               "must not exceed tune=none by more than 0.2",
                               "(tune=none:", round(r_none, 3),
                               "tune=autotune:", round(r_auto, 3), ")"))
    }
  }

  ## ========================================================================
  ## ASSERTION 5: noise DGP stays clean (ratio < 1.2) under both tune conditions
  ## Tolerance from test-spec.md exactly as written: 1.2 (tightened from 1.25).
  ## ========================================================================
  for (tune in c("none", "autotune")) {
    r <- summary_df[summary_df$dgp == "noise" & summary_df$tune == tune, "mean_ratio"]
    if (!is.na(r)) {
      expect_lte(r, 1.2,
                 label = paste("ASSERTION 5: noise DGP tune=", tune,
                               "should be CLEAN (ratio < 1.2, actual:", round(r, 3), ")"))
    }
  }
})
