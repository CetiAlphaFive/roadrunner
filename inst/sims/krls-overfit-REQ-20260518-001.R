## =============================================================================
## krls-overfit-REQ-20260518-001.R
##
## Overfitting Audit: roadrunner::krls()
## Simulation Study — REQ-20260518-001 (2026-05-18)
##
## Purpose:
##   Monte Carlo study measuring train vs test MSE/R2 for roadrunner::krls()
##   across a grid of DGPs, sample sizes, predictor counts, and tuning
##   conditions ("none" = default LOO lambda, "autotune" = inner-CV sigma grid).
##   Primary diagnostic: overfit_ratio = test_MSE / train_MSE.
##
## Usage:
##   Rscript inst/sims/krls-overfit-REQ-20260518-001.R          # full run
##   Rscript inst/sims/krls-overfit-REQ-20260518-001.R quick    # smoke test
##
## Quick mode: n={200,500}, p={5,10}, R=5, tune="none" only.
##   Produces a valid CSV with 5 DGP x 2 n x 2 p x 1 tune x 5 rep = 100 rows.
##
## Outputs:
##   inst/sims/results/krls-overfit-REQ-20260518-001.csv   (per-replicate)
##   inst/sims/results/krls-overfit-REQ-20260518-001.rds   (list: results+spec)
##   inst/sims/results/krls-overfit-plots/                 (5 DGP PNGs + heatmap)
##   [run dir copies are made externally by the simulation harness caller]
##
## Parallelism:
##   Outer loop parallelised via parallel::mclapply over cells; mc.cores = 4.
##   Inside each cell, RcppParallel::setThreadOptions(numThreads = 1) is set
##   to prevent TBB nesting. Sequential fallback on Windows (mc.cores = 1).
##
## Determinism:
##   Per-replicate seed is computed from (dgp_idx, n_idx, p_idx, tune_idx, rep)
##   per the formula in sim-spec.md Section 5. Every row of the output CSV
##   includes the seed so any replicate can be replayed exactly.
##
## Warnings:
##   Captured per replicate via withCallingHandlers; counted in n_warnings.
##   NOT suppressed globally. LOO search boundary warnings on noise-only DGP
##   are expected and are informative data.
## =============================================================================

suppressPackageStartupMessages({
  devtools::load_all("/home/jack/Dropbox/roadrunner", quiet = TRUE)
})

## ---- Command-line args ------------------------------------------------------
args  <- commandArgs(trailingOnly = TRUE)
quick <- length(args) >= 1L && tolower(trimws(args[1])) == "quick"

## ---- Output paths -----------------------------------------------------------
repo_dir  <- "/home/jack/Dropbox/roadrunner"
res_dir   <- file.path(repo_dir, "inst/sims/results")
plot_repo <- file.path(res_dir,  "krls-overfit-plots")
run_dir   <- "/home/jack/.claude/plugins/data/statsclaw-statsclaw/workspace/roadrunner/runs/REQ-20260518-001"
plot_run  <- file.path(run_dir, "sim-plots")

dir.create(res_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(plot_repo, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_run,  showWarnings = FALSE, recursive = TRUE)

csv_repo  <- file.path(res_dir,  "krls-overfit-REQ-20260518-001.csv")
rds_repo  <- file.path(res_dir,  "krls-overfit-REQ-20260518-001.rds")
csv_run   <- file.path(run_dir,  "sim-results.csv")
rds_run   <- file.path(run_dir,  "sim-results.rds")
sum_run   <- file.path(run_dir,  "sim-summary.csv")
sum_rds   <- file.path(run_dir,  "sim-summary.rds")

## ---- DGP definitions --------------------------------------------------------

## All DGPs:
##   Input:  n (integer), p (integer), n_test (integer), seed (integer)
##   Output: list(X_train, y_train, X_test, y_test, sigma_eps)
##   Seeds: set at the start; draws are sequential:
##     1. X_train (n x p)
##     2. eps_train -> y_train
##     3. X_test (n_test x p)
##     4. eps_test -> y_test

## Helper: set column names
.namex <- function(X) {
  colnames(X) <- paste0("V", seq_len(ncol(X)))
  X
}

## DGP 1 — linear-additive
dgp_linear <- function(n, p, n_test, seed) {
  set.seed(seed)
  beta <- c(1.5, -1.0, 0.8, -0.6, 0.4, rep(0, max(0, p - 5L)))
  k    <- min(p, 5L)
  X_tr <- matrix(rnorm(n * p), n, p)
  f_tr <- as.numeric(X_tr[, seq_len(k), drop = FALSE] %*% beta[seq_len(k)])
  sigma_eps <- sd(f_tr) / 2
  y_tr  <- f_tr + rnorm(n, sd = sigma_eps)
  X_te  <- matrix(rnorm(n_test * p), n_test, p)
  f_te  <- as.numeric(X_te[, seq_len(k), drop = FALSE] %*% beta[seq_len(k)])
  y_te  <- f_te + rnorm(n_test, sd = sigma_eps)
  list(X_train = .namex(X_tr), y_train = y_tr,
       X_test  = .namex(X_te), y_test  = y_te,
       sigma_eps = sigma_eps)
}

## DGP 2 — additive-smooth
dgp_additive <- function(n, p, n_test, seed) {
  set.seed(seed)
  .f <- function(X) {
    sin(X[, 1]) +
      0.5 * X[, 2]^2 -
      0.3 * exp(-X[, 3]^2) +
      0.8 * X[, 4] -
      0.5 * X[, 5]
  }
  X_tr <- matrix(rnorm(n * p), n, p)
  f_tr <- .f(X_tr)
  sigma_eps <- sd(f_tr) / 2
  y_tr  <- f_tr + rnorm(n, sd = sigma_eps)
  X_te  <- matrix(rnorm(n_test * p), n_test, p)
  f_te  <- .f(X_te)
  y_te  <- f_te + rnorm(n_test, sd = sigma_eps)
  list(X_train = .namex(X_tr), y_train = y_tr,
       X_test  = .namex(X_te), y_test  = y_te,
       sigma_eps = sigma_eps)
}

## DGP 3 — low-rank interaction
dgp_interaction <- function(n, p, n_test, seed) {
  set.seed(seed)
  .f <- function(X) X[, 1] * X[, 2] + X[, 3] * X[, 4]
  X_tr <- matrix(rnorm(n * p), n, p)
  f_tr <- .f(X_tr)
  sigma_eps <- sd(f_tr) / 2
  y_tr  <- f_tr + rnorm(n, sd = sigma_eps)
  X_te  <- matrix(rnorm(n_test * p), n_test, p)
  f_te  <- .f(X_te)
  y_te  <- f_te + rnorm(n_test, sd = sigma_eps)
  list(X_train = .namex(X_tr), y_train = y_tr,
       X_test  = .namex(X_te), y_test  = y_te,
       sigma_eps = sigma_eps)
}

## DGP 4 — sparse-strong
dgp_sparse <- function(n, p, n_test, seed) {
  set.seed(seed)
  .f <- function(X) 3 * X[, 1] + 2 * X[, 2] - 1.5 * X[, 3]
  X_tr <- matrix(rnorm(n * p), n, p)
  f_tr <- .f(X_tr)
  sigma_eps <- sd(f_tr) / 2
  y_tr  <- f_tr + rnorm(n, sd = sigma_eps)
  X_te  <- matrix(rnorm(n_test * p), n_test, p)
  f_te  <- .f(X_te)
  y_te  <- f_te + rnorm(n_test, sd = sigma_eps)
  list(X_train = .namex(X_tr), y_train = y_tr,
       X_test  = .namex(X_te), y_test  = y_te,
       sigma_eps = sigma_eps)
}

## DGP 5 — noise-only
dgp_noise <- function(n, p, n_test, seed) {
  set.seed(seed)
  X_tr <- matrix(rnorm(n * p), n, p)
  y_tr <- rnorm(n)
  X_te <- matrix(rnorm(n_test * p), n_test, p)
  y_te <- rnorm(n_test)
  list(X_train = .namex(X_tr), y_train = y_tr,
       X_test  = .namex(X_te), y_test  = y_te,
       sigma_eps = 1.0)
}

DGP_FNS <- list(
  linear      = dgp_linear,
  additive    = dgp_additive,
  interaction = dgp_interaction,
  sparse      = dgp_sparse,
  noise       = dgp_noise
)

## ---- Seed formula -----------------------------------------------------------
## Based on sim-spec.md Section 5 formula; index arrays extended to support
## the trimmed-full grid values n in {300,800} and p in {15,30} which are not
## in the original spec arrays. Extended arrays preserve original indices for
## original values and add new positions for new values. This ensures seeds
## are unique, deterministic, and fully documented in output CSV for replay.
##
## Extended n index: c(200,300,500,800,1000,1500) -> n_idx 1..6
## Extended p index: c(5,10,15,25,30,50)          -> p_idx 1..6
make_seed <- function(dgp, n, p, tune, rep) {
  dgp_idx  <- match(dgp,  c("linear","additive","interaction","sparse","noise"))
  tune_idx <- match(tune, c("none","autotune"))
  n_idx    <- match(n,    c(200L, 300L, 500L, 800L, 1000L, 1500L))
  p_idx    <- match(p,    c(5L, 10L, 15L, 25L, 30L, 50L))
  seed <- (dgp_idx  * 10000000L +
           n_idx    *  1000000L +
           p_idx    *   100000L +
           tune_idx *    10000L +
           rep      *        1L) %% (2L^31L - 1L)
  as.integer(seed)
}

## ---- Cell grid --------------------------------------------------------------
## Full-spec dimensions: 5 DGP x 4 n x 4 p x 2 tune = 160 cells
## Trimmed-full dimensions (wall-clock budget):
##   n in {300, 800, 1500}, p in {5, 15, 30}
##   tune in {"none", "autotune"} (both required)
##   R = 25 for tune="none", R = 15 for tune="autotune"
##   Total cells: 5 x 3 x 3 x 2 = 90 cells
##   Total fits: 5 x 3 x 3 x (25+15) = ~1800

ALL_DGPS  <- c("linear", "additive", "interaction", "sparse", "noise")
N_TEST    <- 1000L

if (quick) {
  message("[quick mode] n={200,500}, p={5,10}, tune='none', R=5")
  ALL_N     <- c(200L, 500L)
  ALL_P     <- c(5L, 10L)
  ALL_TUNES <- "none"
} else {
  ## Trimmed full grid (wall-clock-budgeted version of full spec)
  ## Deviates from sim-spec Section 4 as documented in sim-report.md:
  ##   n: {300, 800, 1500} instead of {200, 500, 1000, 1500}
  ##   p: {5, 15, 30}      instead of {5, 10, 25, 50}
  ##   R: 25 (none) / 15 (autotune) instead of 50/40/30
  message("[full-trimmed mode] n={300,800,1500}, p={5,15,30}, tune={none,autotune}")
  ALL_N     <- c(300L, 800L, 1500L)
  ALL_P     <- c(5L, 15L, 30L)
  ALL_TUNES <- c("none", "autotune")
}

cells <- expand.grid(
  dgp  = ALL_DGPS,
  n    = ALL_N,
  p    = ALL_P,
  tune = ALL_TUNES,
  stringsAsFactors = FALSE
)

## Replicate counts per cell
## Trimmed-full: R=25 for tune="none", R=15 for tune="autotune"
## Quick mode:   R=5 (all)
get_R <- function(n, p, tune) {
  if (quick) return(5L)
  ## Trimmed-full fixed R counts
  if (tune == "none") return(25L)
  return(15L)
}
cells$R <- mapply(get_R, cells$n, cells$p, cells$tune)

## ---- Single replicate function ----------------------------------------------

run_one_rep <- function(dgp, n, p, tune, rep, n_test) {
  seed <- make_seed(dgp, n, p, tune, rep)

  ## 1. Draw data
  dat <- DGP_FNS[[dgp]](n, p, n_test, seed)
  X_train <- dat$X_train
  y_train <- dat$y_train
  X_test  <- dat$X_test
  y_test  <- dat$y_test

  ## 2. Fit model — capture warnings, catch errors
  warns <- character(0L)
  t0 <- proc.time()[["elapsed"]]

  fit_result <- withCallingHandlers(
    tryCatch({
      ## Inside each krls call, set nthreads=1 to avoid TBB nesting
      ## when the outer loop is parallel.
      RcppParallel::setThreadOptions(numThreads = 1L)
      if (tune == "none") {
        roadrunner::krls(
          X          = X_train,
          y          = y_train,
          derivative = FALSE,
          vcov       = FALSE,
          trace      = 0L
        )
      } else {
        ## autotune = TRUE uses default sigma grid p * c(0.25,.5,1,2,4,8)
        roadrunner::krls(
          X          = X_train,
          y          = y_train,
          derivative = FALSE,
          vcov       = FALSE,
          autotune   = TRUE,
          trace      = 0L
        )
      }
    }, error = function(e) {
      structure(list(msg = conditionMessage(e)), class = "krls_error")
    }),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  wall_clock_s <- proc.time()[["elapsed"]] - t0
  n_warnings   <- length(warns)
  n_errors     <- if (inherits(fit_result, "krls_error")) 1L else 0L

  ## 3. Compute metrics
  if (n_errors == 1L) {
    return(data.frame(
      dgp           = dgp,
      n             = n,
      p             = p,
      tune          = tune,
      rep           = rep,
      seed          = seed,
      train_R2      = NA_real_,
      test_R2       = NA_real_,
      train_MSE     = NA_real_,
      test_MSE      = NA_real_,
      overfit_ratio = NA_real_,
      overfit_gap   = NA_real_,
      sigma_hat     = NA_real_,
      lambda_hat    = NA_real_,
      df_eff        = NA_real_,
      wall_clock_s  = wall_clock_s,
      n_warnings    = n_warnings,
      n_errors      = n_errors,
      stringsAsFactors = FALSE
    ))
  }

  fit <- fit_result

  ## Training predictions: fit$fitted is a matrix column
  yhat_train <- as.numeric(fit$fitted)

  ## Test predictions: predict returns list(fit=matrix, ...)
  pred_test  <- predict(fit, newdata = X_test)
  yhat_test  <- as.numeric(pred_test$fit)

  ## MSE and R2
  train_MSE <- mean((y_train - yhat_train)^2)
  test_MSE  <- mean((y_test  - yhat_test)^2)
  train_R2  <- 1 - train_MSE / var(y_train)
  test_R2   <- 1 - test_MSE  / var(y_test)

  overfit_ratio <- test_MSE / train_MSE
  overfit_gap   <- train_R2 - test_R2

  sigma_hat  <- fit$sigma
  lambda_hat <- fit$lambda
  df_eff     <- NA_real_

  data.frame(
    dgp           = dgp,
    n             = n,
    p             = p,
    tune          = tune,
    rep           = rep,
    seed          = seed,
    train_R2      = train_R2,
    test_R2       = test_R2,
    train_MSE     = train_MSE,
    test_MSE      = test_MSE,
    overfit_ratio = overfit_ratio,
    overfit_gap   = overfit_gap,
    sigma_hat     = sigma_hat,
    lambda_hat    = lambda_hat,
    df_eff        = df_eff,
    wall_clock_s  = wall_clock_s,
    n_warnings    = n_warnings,
    n_errors      = n_errors,
    stringsAsFactors = FALSE
  )
}

## ---- Cell runner (all reps for one cell) ------------------------------------

run_cell <- function(i, cells, n_test) {
  row  <- cells[i, ]
  dgp  <- row$dgp
  n    <- as.integer(row$n)
  p    <- as.integer(row$p)
  tune <- row$tune
  R    <- as.integer(row$R)

  n_cells <- nrow(cells)
  t_cell0 <- proc.time()[["elapsed"]]
  message(sprintf("[Cell %d/%d] dgp=%-12s n=%4d p=%2d tune=%-9s R=%d ...",
                  i, n_cells, dgp, n, p, tune, R))

  reps_list <- vector("list", R)
  for (r in seq_len(R)) {
    reps_list[[r]] <- tryCatch(
      run_one_rep(dgp, n, p, tune, r, n_test),
      error = function(e) {
        message("  !! rep ", r, " outer error: ", conditionMessage(e))
        data.frame(
          dgp = dgp, n = n, p = p, tune = tune, rep = r,
          seed          = make_seed(dgp, n, p, tune, r),
          train_R2      = NA_real_, test_R2      = NA_real_,
          train_MSE     = NA_real_, test_MSE     = NA_real_,
          overfit_ratio = NA_real_, overfit_gap  = NA_real_,
          sigma_hat     = NA_real_, lambda_hat   = NA_real_,
          df_eff        = NA_real_, wall_clock_s = NA_real_,
          n_warnings    = NA_integer_, n_errors  = 1L,
          stringsAsFactors = FALSE
        )
      }
    )
  }

  cell_df    <- do.call(rbind, reps_list)
  t_cell_el  <- proc.time()[["elapsed"]] - t_cell0
  good_rows  <- cell_df[!is.na(cell_df$overfit_ratio), ]
  mean_ratio <- if (nrow(good_rows) > 0) mean(good_rows$overfit_ratio) else NA_real_
  message(sprintf("  -> mean_ratio=%.3f  completed=%d/%d  (wall %.1fs)",
                  ifelse(is.na(mean_ratio), NaN, mean_ratio),
                  nrow(good_rows), R, t_cell_el))
  cell_df
}

## ---- Main simulation loop ---------------------------------------------------

message("=============================================================")
message("krls overfitting audit — REQ-20260518-001")
message(sprintf("Mode: %s", if (quick) "QUICK" else "FULL-TRIMMED"))
message(sprintf("Cells: %d  Total reps (approx): %d",
                nrow(cells), sum(cells$R)))
if (!quick) {
  message("Grid: n={300,800,1500} x p={5,15,30} x tune={none,autotune} x 5 DGPs")
  message("R: 25 (tune=none) / 15 (tune=autotune)")
}
message(sprintf("Parallelism: outer mclapply, mc.cores = 4; nthreads=1 per call"))
message("=============================================================")

## Determine parallelism
use_parallel <- .Platform$OS.type != "windows"
n_cores      <- if (use_parallel) 4L else 1L

t_total0 <- proc.time()[["elapsed"]]

if (use_parallel) {
  ## Parallel outer loop — each cell is one mclapply job.
  ## Each job runs all R replicates for that cell sequentially.
  ## Dispatch in two batches so we can checkpoint timing after 10 cells.
  n_cells_total <- nrow(cells)
  probe_batch   <- min(10L, n_cells_total)

  ## First batch: cells 1..probe_batch (timing probe)
  message(sprintf("[timing probe] Dispatching first %d cells ...", probe_batch))
  t_probe0 <- proc.time()[["elapsed"]]
  cell_results_probe <- parallel::mclapply(
    seq_len(probe_batch),
    function(i) run_cell(i, cells, N_TEST),
    mc.cores = n_cores,
    mc.set.seed = FALSE
  )
  t_probe_el <- proc.time()[["elapsed"]] - t_probe0
  message(sprintf("[timing probe] %d cells done in %.1f s (%.1f min)",
                  probe_batch, t_probe_el, t_probe_el / 60))

  ## Project total wall-clock
  if (n_cells_total > probe_batch) {
    ## Reps per cell in probe vs rest; use a rep-weighted projection
    probe_reps <- sum(cells$R[seq_len(probe_batch)])
    remain_reps <- sum(cells$R[(probe_batch + 1L):n_cells_total])
    t_per_rep   <- t_probe_el / probe_reps
    t_remain    <- t_per_rep * remain_reps / n_cores   ## parallel speedup
    t_project   <- t_probe_el + t_remain
    message(sprintf("[timing probe] Projected total: %.1f min (remaining: %.1f min)",
                    t_project / 60, t_remain / 60))

    ## If projected total > 45 min, reduce R uniformly by a factor to stay in budget
    if (t_project / 60 > 45) {
      scale_factor <- 45 * 60 / t_project
      old_R <- cells$R
      cells$R <- pmax(5L, as.integer(floor(cells$R * scale_factor)))
      message(sprintf("[wall-clock guard] Projected > 45 min. Scaling R by %.2f",
                      scale_factor))
      message(sprintf("  R range before: %d-%d; after: %d-%d",
                      min(old_R), max(old_R), min(cells$R), max(cells$R)))
    }
  }

  ## Remaining cells
  if (n_cells_total > probe_batch) {
    cell_results_rest <- parallel::mclapply(
      (probe_batch + 1L):n_cells_total,
      function(i) run_cell(i, cells, N_TEST),
      mc.cores = n_cores,
      mc.set.seed = FALSE
    )
    cell_results_list <- c(cell_results_probe, cell_results_rest)
  } else {
    cell_results_list <- cell_results_probe
  }
} else {
  cell_results_list <- lapply(
    seq_len(nrow(cells)),
    function(i) run_cell(i, cells, N_TEST)
  )
}

results <- do.call(rbind, cell_results_list)
rownames(results) <- NULL

t_total_el <- proc.time()[["elapsed"]] - t_total0
message(sprintf("\nTotal wall-clock: %.1f seconds (%.1f min)",
                t_total_el, t_total_el / 60))

## ---- Save per-replicate results ---------------------------------------------

write.csv(results, csv_repo, row.names = FALSE, quote = FALSE)
write.csv(results, csv_run,  row.names = FALSE, quote = FALSE)
message("Wrote per-replicate CSV: ", csv_repo)

rds_obj <- list(
  results = results,
  spec    = paste0("sim-spec.md REQ-20260518-001 (2026-05-18); ",
                   "roadrunner::krls overfitting audit"),
  session = sessionInfo()
)
saveRDS(rds_obj, rds_repo)
saveRDS(rds_obj, rds_run)
message("Wrote per-replicate RDS: ", rds_repo)

## ---- Summary aggregation ----------------------------------------------------

## Verdict function — applied per (dgp, n, p, tune) cell
verdict_fn <- function(mean_ratio, ci_lo) {
  if (!is.finite(mean_ratio)) return("NO_DATA")
  if (mean_ratio > 1.5 && ci_lo > 1.0) "OVERFIT"
  else if (mean_ratio > 1.2)            "BORDERLINE"
  else                                   "CLEAN"
}

## Bootstrap CI on mean overfit_ratio (B=2000 draws)
boot_ci_mean <- function(x, B = 2000L) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(c(ci_lo = NA_real_, ci_hi = NA_real_))
  boot_means <- replicate(B, mean(sample(x, length(x), replace = TRUE)))
  c(ci_lo = as.numeric(quantile(boot_means, 0.025)),
    ci_hi = as.numeric(quantile(boot_means, 0.975)))
}

## Compute summary per cell
cell_keys <- unique(results[, c("dgp", "n", "p", "tune")])
summary_rows <- vector("list", nrow(cell_keys))

for (i in seq_len(nrow(cell_keys))) {
  key   <- cell_keys[i, ]
  sub   <- results[results$dgp  == key$dgp &
                   results$n    == key$n   &
                   results$p    == key$p   &
                   results$tune == key$tune, ]
  good  <- sub[!is.na(sub$overfit_ratio), ]
  R_act <- nrow(sub)

  ratio_vec   <- good$overfit_ratio
  mean_ratio  <- if (length(ratio_vec) > 0) mean(ratio_vec) else NA_real_
  se_ratio    <- if (length(ratio_vec) > 1) sd(ratio_vec) / sqrt(length(ratio_vec)) else NA_real_

  ci_vals     <- boot_ci_mean(ratio_vec)
  ci_lo       <- ci_vals["ci_lo"]
  ci_hi       <- ci_vals["ci_hi"]

  verd        <- verdict_fn(mean_ratio, ci_lo)

  summary_rows[[i]] <- data.frame(
    dgp           = key$dgp,
    n             = key$n,
    p             = key$p,
    tune          = key$tune,
    R_actual      = R_act,
    mean_ratio    = mean_ratio,
    se_ratio      = se_ratio,
    ci_lo         = as.numeric(ci_lo),
    ci_hi         = as.numeric(ci_hi),
    mean_gap      = if (nrow(good) > 0) mean(good$overfit_gap,  na.rm = TRUE) else NA_real_,
    mean_train_R2 = if (nrow(good) > 0) mean(good$train_R2,     na.rm = TRUE) else NA_real_,
    mean_test_R2  = if (nrow(good) > 0) mean(good$test_R2,      na.rm = TRUE) else NA_real_,
    mean_sigma    = if (nrow(good) > 0) mean(good$sigma_hat,    na.rm = TRUE) else NA_real_,
    mean_lambda   = if (nrow(good) > 0) mean(good$lambda_hat,   na.rm = TRUE) else NA_real_,
    mean_wall_s   = if (nrow(good) > 0) mean(good$wall_clock_s, na.rm = TRUE) else NA_real_,
    verdict       = verd,
    stringsAsFactors = FALSE
  )
}

summary_df <- do.call(rbind, summary_rows)
rownames(summary_df) <- NULL

write.csv(summary_df, sum_run, row.names = FALSE, quote = FALSE)
saveRDS(summary_df, sum_rds)
message("Wrote summary CSV: ", sum_run)
message("Wrote summary RDS: ", sum_rds)

## ---- Print summary table ----------------------------------------------------

message("\n======= VERDICT SUMMARY =======")
for (tun in unique(summary_df$tune)) {
  message(sprintf("\ntune = %s", tun))
  message(sprintf("%-12s  %4s  %2s  %8s  %8s  %8s  %10s",
                  "dgp", "n", "p", "mean_ratio", "ci_lo", "ci_hi", "verdict"))
  sub <- summary_df[summary_df$tune == tun, ]
  sub <- sub[order(sub$dgp, sub$n, sub$p), ]
  for (i in seq_len(nrow(sub))) {
    r <- sub[i, ]
    message(sprintf("%-12s  %4d  %2d  %8.3f  %8.3f  %8.3f  %10s",
                    r$dgp, r$n, r$p,
                    ifelse(is.na(r$mean_ratio), NaN, r$mean_ratio),
                    ifelse(is.na(r$ci_lo), NaN, r$ci_lo),
                    ifelse(is.na(r$ci_hi), NaN, r$ci_hi),
                    r$verdict))
  }
}

## ---- Plots ------------------------------------------------------------------

## Helper: draw one DGP panel (overfit_ratio vs n, lines by p, 2 tune sub-panels)
.plot_dgp <- function(dgp_name, summary_df, plot_path) {
  sub   <- summary_df[summary_df$dgp == dgp_name, ]
  tunes <- unique(sub$tune)
  ps    <- sort(unique(sub$p))
  ns    <- sort(unique(sub$n))

  n_tune <- length(tunes)
  ## 2 rows x 2 cols: left col = overfit_ratio, right col = overfit_gap
  ## rows  = tune conditions
  grDevices::png(plot_path, width = 1200, height = 600 * n_tune, res = 120)
  graphics::par(mfrow = c(n_tune, 2),
                mar = c(4, 4, 3, 1) + 0.1,
                oma = c(0, 0, 3, 0))

  p_cols <- grDevices::rainbow(length(ps))

  for (ti in seq_along(tunes)) {
    tun <- tunes[ti]
    st  <- sub[sub$tune == tun, ]

    ## ----- Panel left: overfit_ratio vs n ---------------------------------
    y_rat_all <- c(st$ci_lo, st$ci_hi)
    y_rat_all <- y_rat_all[is.finite(y_rat_all)]
    ylim_rat <- if (length(y_rat_all) > 0)
      range(c(0.8, 2.0, y_rat_all)) else c(0.8, 2.0)

    graphics::plot(NA, xlim = range(ns), ylim = ylim_rat,
                   xlab = "n (training)", ylab = "overfit_ratio (test_MSE / train_MSE)",
                   main = sprintf("tune=%s | overfit ratio", tun),
                   log = "x", xaxt = "n")
    graphics::axis(1, at = ns, labels = ns)
    graphics::abline(h = 1.0, col = "black", lty = 1, lwd = 1.5)
    graphics::abline(h = 1.5, col = "red",   lty = 2, lwd = 1.5)
    graphics::abline(h = 1.2, col = "orange",lty = 3, lwd = 1.2)

    for (pi_idx in seq_along(ps)) {
      pv <- ps[pi_idx]
      sp <- st[st$p == pv, ]
      sp <- sp[order(sp$n), ]
      if (nrow(sp) < 1) next
      graphics::lines(sp$n, sp$mean_ratio, col = p_cols[pi_idx], lwd = 2, pch = 16, type = "b")
      ## CI ribbon (approximate via lines)
      graphics::lines(sp$n, sp$ci_lo, col = p_cols[pi_idx], lty = 3, lwd = 1)
      graphics::lines(sp$n, sp$ci_hi, col = p_cols[pi_idx], lty = 3, lwd = 1)
    }
    graphics::legend("topright", legend = paste0("p=", ps), col = p_cols,
                     lwd = 2, lty = 1, cex = 0.8, bty = "n")

    ## ----- Panel right: overfit_gap vs n ----------------------------------
    y_gap_all <- c(st$mean_gap - 2 * st$se_ratio, st$mean_gap + 2 * st$se_ratio)
    y_gap_all <- y_gap_all[is.finite(y_gap_all)]
    ylim_gap  <- if (length(y_gap_all) > 0)
      range(c(-0.1, 0.5, y_gap_all)) else c(-0.1, 0.5)

    graphics::plot(NA, xlim = range(ns), ylim = ylim_gap,
                   xlab = "n (training)", ylab = "overfit_gap (train_R2 - test_R2)",
                   main = sprintf("tune=%s | R2 gap", tun),
                   log = "x", xaxt = "n")
    graphics::axis(1, at = ns, labels = ns)
    graphics::abline(h = 0, col = "black", lty = 1, lwd = 1.5)

    for (pi_idx in seq_along(ps)) {
      pv <- ps[pi_idx]
      sp <- st[st$p == pv, ]
      sp <- sp[order(sp$n), ]
      if (nrow(sp) < 1) next
      graphics::lines(sp$n, sp$mean_gap, col = p_cols[pi_idx], lwd = 2, pch = 16, type = "b")
    }
    graphics::legend("topright", legend = paste0("p=", ps), col = p_cols,
                     lwd = 2, lty = 1, cex = 0.8, bty = "n")
  }

  graphics::mtext(sprintf("DGP: %s — overfit_ratio and R2 gap vs n", dgp_name),
                  outer = TRUE, cex = 1.3, font = 2)
  grDevices::dev.off()
  invisible(plot_path)
}

## Generate per-DGP plots
for (dgp_nm in ALL_DGPS) {
  for (pdir in c(plot_repo, plot_run)) {
    ppath <- file.path(pdir, sprintf("overfit-ratio-%s.png", dgp_nm))
    .plot_dgp(dgp_nm, summary_df, ppath)
    message("Wrote plot: ", ppath)
  }
}

## Heatmap: verdict grid (dgp x p rows, n columns, faceted by tune)
.plot_heatmap <- function(summary_df, plot_path) {
  tunes  <- sort(unique(summary_df$tune))
  dgps   <- rev(c("linear","additive","interaction","sparse","noise"))
  ps     <- sort(unique(summary_df$p))
  ns     <- sort(unique(summary_df$n))
  n_tune <- length(tunes)

  verdict_col <- function(v) {
    switch(v,
      OVERFIT    = "tomato",
      BORDERLINE = "orange",
      CLEAN      = "seagreen3",
      NO_DATA    = "grey70",
      "white"
    )
  }

  ## rows: dgp x p combinations; cols: n
  rows  <- expand.grid(dgp = dgps, p = ps, stringsAsFactors = FALSE)
  rows  <- rows[order(rows$dgp, rows$p), ]
  nrows <- nrow(rows)
  ncols <- length(ns)

  ## y-axis labels
  row_labels <- sprintf("%s/p=%d", rows$dgp, rows$p)

  grDevices::png(plot_path, width = 200 * (ncols + 2),
                            height = 70 * nrows * n_tune + 120,
                            res = 96)

  graphics::par(mfrow = c(n_tune, 1),
                mar   = c(2, 8, 3, 1) + 0.1,
                oma   = c(2, 0, 3, 0))

  for (tun in tunes) {
    st <- summary_df[summary_df$tune == tun, ]

    graphics::plot.new()
    graphics::plot.window(xlim = c(0.5, ncols + 0.5), ylim = c(0.5, nrows + 0.5))
    graphics::title(main = sprintf("tune = %s", tun), font.main = 2)
    graphics::axis(2, at = seq_len(nrows), labels = row_labels, las = 1, cex.axis = 0.7)
    graphics::axis(1, at = seq_along(ns), labels = ns)
    graphics::mtext("n (training)", side = 1, line = 1.5, cex = 0.8)

    for (ri in seq_len(nrows)) {
      dg <- rows$dgp[ri]
      pv <- rows$p[ri]
      for (ci in seq_along(ns)) {
        nv  <- ns[ci]
        sub <- st[st$dgp == dg & st$p == pv & st$n == nv, ]
        vc  <- if (nrow(sub) == 1) verdict_col(sub$verdict) else "white"
        graphics::rect(ci - 0.45, ri - 0.45, ci + 0.45, ri + 0.45,
                       col = vc, border = "white")
        if (nrow(sub) == 1 && !is.na(sub$mean_ratio)) {
          graphics::text(ci, ri, sprintf("%.2f", sub$mean_ratio), cex = 0.55)
        }
      }
    }
    ## Legend
    graphics::legend("bottomright",
                     legend = c("OVERFIT", "BORDERLINE", "CLEAN"),
                     fill   = c("tomato", "orange", "seagreen3"),
                     bty = "n", cex = 0.75)
  }

  graphics::mtext("Verdict heatmap: overfit_ratio per (DGP, p, n, tune)",
                  outer = TRUE, cex = 1.1, font = 2)
  grDevices::dev.off()
  invisible(plot_path)
}

for (pdir in c(plot_repo, plot_run)) {
  hp <- file.path(pdir, "overfit-heatmap.png")
  .plot_heatmap(summary_df, hp)
  message("Wrote heatmap: ", hp)
}

message("\n=== DONE ===")
message(sprintf("Total wall-clock: %.1f s (%.1f min)", t_total_el, t_total_el / 60))
message(sprintf("Results: %s", csv_repo))
message(sprintf("Summary: %s", sum_run))
