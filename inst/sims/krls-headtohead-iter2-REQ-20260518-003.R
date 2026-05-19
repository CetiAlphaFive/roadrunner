## krls-headtohead-iter2-REQ-20260518-003.R
##
## Iter-2 verification: re-run the iter-0 head-to-head after the
## .krls_sigma_anchor patch (median(d2) -> sqrt(median(d2) * ncol(Xs))).
##
## Uses devtools::load_all() to pick up the in-source patch.
## Same seeds, n=500, p=10, R=5, paired, sequential.
##
## Outputs written to REQ-20260518-003/ as iter-2 variants.
## iter-0 artifacts are NOT overwritten.

suppressPackageStartupMessages({
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all("/home/jack/Dropbox/roadrunner", quiet = TRUE)
  } else {
    library(roadrunner)
  }
  library(KRLS)
  library(RcppParallel)
})

RcppParallel::setThreadOptions(numThreads = 1L)

## -----------------------------------------------------------------------
## Paths
## -----------------------------------------------------------------------
RUN_DIR   <- "/home/jack/.claude/plugins/data/statsclaw-statsclaw/workspace/roadrunner/runs/REQ-20260518-003"
PLOTS_DIR <- file.path(RUN_DIR, "headtohead-plots")        # reuse same dir
CSV_PATH  <- file.path(RUN_DIR, "headtohead-results-iter2.csv")
RDS_PATH  <- file.path(RUN_DIR, "headtohead-results-iter2.rds")
SUM_CSV   <- file.path(RUN_DIR, "headtohead-summary-iter2.csv")
REPORT_MD <- file.path(RUN_DIR, "headtohead-report-iter2.md")
ITER0_SUM <- file.path(RUN_DIR, "headtohead-summary.csv")     # for delta col

N_TRAIN    <- 500L
N_TEST     <- 1000L
P          <- 10L
R_REPS     <- 5L
WALL_CAP_S <- 300
MASTER_SEED <- 20260518L

dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
wall_start <- proc.time()[["elapsed"]]

## -----------------------------------------------------------------------
## Load iter-0 summary for delta computation
## -----------------------------------------------------------------------
iter0_sum <- read.csv(ITER0_SUM, stringsAsFactors = FALSE)

## -----------------------------------------------------------------------
## DGP definitions (identical to iter-0 script)
## -----------------------------------------------------------------------
dgp_linear <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  f <- X[, 1] + 0.5 * X[, 2] - 0.3 * X[, 3]
  eps_sd <- sd(f) / 2
  y <- f + rnorm(n, sd = eps_sd)
  list(X = X, y = y, f = f, eps_sd = eps_sd)
}

dgp_additive <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  f <- sin(X[, 1]) + X[, 2]^2 + exp(-abs(X[, 3]))
  eps_sd <- sd(f) / 2
  y <- f + rnorm(n, sd = eps_sd)
  list(X = X, y = y, f = f, eps_sd = eps_sd)
}

dgp_interaction <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  f <- X[, 1] * X[, 2] + sin(X[, 3]) * X[, 4]
  eps_sd <- sd(f) / 2
  y <- f + rnorm(n, sd = eps_sd)
  list(X = X, y = y, f = f, eps_sd = eps_sd)
}

dgp_sparse <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  f <- 2 * X[, 1] + 0.5 * X[, 2]
  eps_sd <- sd(f) / 2
  y <- f + rnorm(n, sd = eps_sd)
  list(X = X, y = y, f = f, eps_sd = eps_sd)
}

dgp_noise <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  f <- rep(0, n)
  y <- rnorm(n, sd = 1)
  list(X = X, y = y, f = f, eps_sd = 1)
}

dgp_friedman1 <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(runif(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  f <- 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 +
       10 * X[, 4] + 5 * X[, 5]
  eps_sd <- sd(f) / 2
  y <- f + rnorm(n, sd = eps_sd)
  list(X = X, y = y, f = f, eps_sd = eps_sd)
}

dgp_friedman2 <- function(n, p, seed) {
  set.seed(seed)
  x1 <- runif(n, 0, 100)
  x2 <- runif(n, 40 * pi, 560 * pi)
  x3 <- runif(n, 0, 1)
  x4 <- runif(n, 1, 11)
  X_extra <- matrix(runif(n * (p - 4)), n, p - 4)
  X <- cbind(x1, x2, x3, x4, X_extra)
  colnames(X) <- paste0("x", seq_len(p))
  f <- sqrt(x1^2 + (x2 * x3 - 1 / (x2 * x4))^2)
  eps_sd <- sd(f) / 2
  y <- f + rnorm(n, sd = eps_sd)
  list(X = X, y = y, f = f, eps_sd = eps_sd)
}

dgp_friedman3 <- function(n, p, seed) {
  set.seed(seed)
  x1 <- runif(n, 0, 100)
  x2 <- runif(n, 40 * pi, 560 * pi)
  x3 <- runif(n, 0, 1)
  x4 <- runif(n, 1, 11)
  X_extra <- matrix(runif(n * (p - 4)), n, p - 4)
  X <- cbind(x1, x2, x3, x4, X_extra)
  colnames(X) <- paste0("x", seq_len(p))
  f <- atan((x2 * x3 - 1 / (x2 * x4)) / x1)
  eps_sd <- sd(f) / 2
  y <- f + rnorm(n, sd = eps_sd)
  list(X = X, y = y, f = f, eps_sd = eps_sd)
}

dgp_sin_sum <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  f <- rowSums(sin(X[, 1:5]))
  eps_sd <- sd(f) / 2
  y <- f + rnorm(n, sd = eps_sd)
  list(X = X, y = y, f = f, eps_sd = eps_sd)
}

dgp_exp_decay <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  f <- exp(-rowSums(X^2) / p)
  eps_sd <- sd(f) / 2
  y <- f + rnorm(n, sd = eps_sd)
  list(X = X, y = y, f = f, eps_sd = eps_sd)
}

dgp_poly2 <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  f <- rowSums(X[, 1:5]^2)
  eps_sd <- sd(f) / 2
  y <- f + rnorm(n, sd = eps_sd)
  list(X = X, y = y, f = f, eps_sd = eps_sd)
}

dgp_tanh_interaction <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  f <- tanh(X[, 1] * X[, 2]) + tanh(X[, 3] * X[, 4])
  eps_sd <- sd(f) / 2
  y <- f + rnorm(n, sd = eps_sd)
  list(X = X, y = y, f = f, eps_sd = eps_sd)
}

dgp_mixture <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  f <- ifelse(X[, 1] > 0, sin(X[, 2]), cos(X[, 2]))
  eps_sd <- sd(f) / 2
  y <- f + rnorm(n, sd = eps_sd)
  list(X = X, y = y, f = f, eps_sd = eps_sd)
}

dgp_monotone <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  f <- log(1 + exp(X[, 1] + X[, 2]))
  eps_sd <- sd(f) / 2
  y <- f + rnorm(n, sd = eps_sd)
  list(X = X, y = y, f = f, eps_sd = eps_sd)
}

dgp_heterosked <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  f <- X[, 1] + X[, 2]
  eps_scale <- 1 + abs(X[, 1])
  eps_sd_base <- sd(f) / (2 * mean(eps_scale))
  y <- f + rnorm(n, sd = eps_sd_base) * eps_scale
  list(X = X, y = y, f = f, eps_sd = eps_sd_base)
}

## -----------------------------------------------------------------------
## DGP catalog (same order as iter-0)
## -----------------------------------------------------------------------
DGPS <- list(
  list(name = "linear",           fn = dgp_linear),
  list(name = "additive",         fn = dgp_additive),
  list(name = "interaction",      fn = dgp_interaction),
  list(name = "sparse",           fn = dgp_sparse),
  list(name = "noise",            fn = dgp_noise),
  list(name = "friedman1",        fn = dgp_friedman1),
  list(name = "friedman2",        fn = dgp_friedman2),
  list(name = "friedman3",        fn = dgp_friedman3),
  list(name = "sin-sum",          fn = dgp_sin_sum),
  list(name = "exp-decay",        fn = dgp_exp_decay),
  list(name = "poly2",            fn = dgp_poly2),
  list(name = "tanh-interaction", fn = dgp_tanh_interaction),
  list(name = "mixture",          fn = dgp_mixture),
  list(name = "monotone",         fn = dgp_monotone),
  list(name = "heterosked",       fn = dgp_heterosked)
)

## -----------------------------------------------------------------------
## Helpers
## -----------------------------------------------------------------------
predict_vec <- function(fit, X_test) {
  pred <- predict(fit, newdata = X_test)
  if (is.list(pred)) {
    if (!is.null(pred$fit))      return(as.numeric(pred$fit))
    if (!is.null(pred$yfitted))  return(as.numeric(pred$yfitted))
  }
  as.numeric(pred)
}

mse <- function(a, b) mean((a - b)^2)
r2  <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)

run_one <- function(pkg, X_train, y_train, X_test, y_test) {
  n_warn <- 0L
  n_err  <- 0L
  sigma_hat  <- NA_real_
  lambda_hat <- NA_real_
  train_MSE  <- NA_real_
  test_MSE   <- NA_real_
  train_R2   <- NA_real_
  test_R2    <- NA_real_

  t0 <- proc.time()[["elapsed"]]

  fit <- tryCatch(
    withCallingHandlers(
      {
        if (pkg == "rr") {
          roadrunner::krls(X_train, y_train,
                           derivative = FALSE, vcov = FALSE,
                           trace = 0L, nthreads = 1L)
        } else {
          KRLS::krls(X_train, y_train,
                     derivative = FALSE, vcov = FALSE,
                     print.level = 0L)
        }
      },
      warning = function(w) {
        n_warn <<- n_warn + 1L
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      n_err <<- n_err + 1L
      NULL
    }
  )

  wall_s <- proc.time()[["elapsed"]] - t0

  if (!is.null(fit)) {
    sigma_hat  <- as.numeric(fit$sigma)
    lambda_hat <- as.numeric(fit$lambda)
    y_fit      <- as.numeric(fit$fitted)
    train_MSE  <- mse(y_train, y_fit)
    train_R2   <- r2(y_train, y_fit)

    y_pred <- tryCatch(
      predict_vec(fit, X_test),
      error = function(e) { n_err <<- n_err + 1L; NULL }
    )

    if (!is.null(y_pred) && length(y_pred) == length(y_test)) {
      test_MSE <- mse(y_test, y_pred)
      test_R2  <- r2(y_test, y_pred)
    }
  }

  list(sigma_hat = sigma_hat, lambda_hat = lambda_hat,
       train_MSE = train_MSE, test_MSE = test_MSE,
       train_R2 = train_R2, test_R2 = test_R2,
       wall_s = wall_s, n_warnings = n_warn, n_errors = n_err)
}

make_row <- function(dgp_name, rep, seed, pkg, metrics) {
  data.frame(
    dgp        = dgp_name, rep = rep, seed = seed, pkg = pkg,
    sigma_hat  = metrics$sigma_hat, lambda_hat = metrics$lambda_hat,
    train_MSE  = metrics$train_MSE, test_MSE   = metrics$test_MSE,
    train_R2   = metrics$train_R2,  test_R2    = metrics$test_R2,
    overfit_ratio = metrics$test_MSE / metrics$train_MSE,
    wall_s     = metrics$wall_s,
    n_warnings = metrics$n_warnings, n_errors  = metrics$n_errors,
    stringsAsFactors = FALSE
  )
}

csv_header_written <- FALSE
append_row <- function(row_df) {
  if (!csv_header_written || !file.exists(CSV_PATH)) {
    write.csv(row_df, CSV_PATH, row.names = FALSE)
    csv_header_written <<- TRUE
  } else {
    write.table(row_df, CSV_PATH, row.names = FALSE, col.names = FALSE,
                sep = ",", append = TRUE)
  }
}

## -----------------------------------------------------------------------
## SMOKE TEST: 1 DGP x 2 reps, both packages
## -----------------------------------------------------------------------
cat("=== SMOKE TEST (iter-2) ===\n")
smoke_seed <- MASTER_SEED + 999999L
set.seed(smoke_seed)
X_sm      <- matrix(rnorm(200 * P), 200, P)
y_sm      <- X_sm[, 1] + rnorm(200, sd = 0.5)
X_sm_test <- matrix(rnorm(100 * P), 100, P)
y_sm_test <- X_sm_test[, 1] + rnorm(100, sd = 0.5)

sm_rr   <- run_one("rr",   X_sm, y_sm, X_sm_test, y_sm_test)
sm_krls <- run_one("krls", X_sm, y_sm, X_sm_test, y_sm_test)

stopifnot(
  "smoke: rr test_MSE not finite"   = is.finite(sm_rr$test_MSE),
  "smoke: krls test_MSE not finite" = is.finite(sm_krls$test_MSE),
  "smoke: rr sigma not positive"    = sm_rr$sigma_hat > 0,
  "smoke: krls sigma not positive"  = sm_krls$sigma_hat > 0
)
cat(sprintf("  rr   test_MSE=%.4f sigma=%.3f lambda=%.6f wall=%.2fs\n",
            sm_rr$test_MSE, sm_rr$sigma_hat, sm_rr$lambda_hat, sm_rr$wall_s))
cat(sprintf("  krls test_MSE=%.4f sigma=%.3f lambda=%.6f wall=%.2fs\n",
            sm_krls$test_MSE, sm_krls$sigma_hat, sm_krls$lambda_hat, sm_krls$wall_s))
cat("  SMOKE TEST PASSED\n\n")

## -----------------------------------------------------------------------
## FULL GRID
## -----------------------------------------------------------------------
cat(sprintf("Starting full grid: %d DGPs x %d reps x 2 pkgs\n",
            length(DGPS), R_REPS))
cat(sprintf("n_train=%d  n_test=%d  p=%d\n\n", N_TRAIN, N_TEST, P))

all_rows <- list()
row_idx  <- 0L
cell_done <- 0L

for (dgp_i in seq_along(DGPS)) {
  dgp_name <- DGPS[[dgp_i]]$name
  dgp_fn   <- DGPS[[dgp_i]]$fn

  cat(sprintf("[%d/%d] DGP: %s\n", dgp_i, length(DGPS), dgp_name))

  for (r in seq_len(R_REPS)) {
    elapsed <- proc.time()[["elapsed"]] - wall_start
    if (elapsed > WALL_CAP_S - 15) {
      cat(sprintf("  WALL-CLOCK CAP approaching (%.0fs). Aborting.\n", elapsed))
      break
    }

    ## Identical seed derivation as iter-0
    train_seed <- MASTER_SEED + (dgp_i - 1L) * R_REPS * 3L + (r - 1L) * 3L + 1L
    test_seed  <- MASTER_SEED + (dgp_i - 1L) * R_REPS * 3L + (r - 1L) * 3L + 2L

    d_train <- dgp_fn(N_TRAIN, P, train_seed)
    d_test  <- dgp_fn(N_TEST,  P, test_seed)

    ## ---- roadrunner ----
    m_rr   <- run_one("rr",   d_train$X, d_train$y, d_test$X, d_test$y)
    row_rr <- make_row(dgp_name, r, train_seed, "rr", m_rr)
    row_idx <- row_idx + 1L; all_rows[[row_idx]] <- row_rr
    append_row(row_rr)

    ## ---- KRLS ----
    m_krls   <- run_one("krls", d_train$X, d_train$y, d_test$X, d_test$y)
    row_krls <- make_row(dgp_name, r, train_seed, "krls", m_krls)
    row_idx <- row_idx + 1L; all_rows[[row_idx]] <- row_krls
    append_row(row_krls)

    mse_ratio <- m_rr$test_MSE / m_krls$test_MSE
    cat(sprintf("  rep%d: rr_MSE=%.4f krls_MSE=%.4f ratio=%.3f  [rr_sig=%.2f krls_sig=%.2f]\n",
                r, m_rr$test_MSE, m_krls$test_MSE, mse_ratio,
                m_rr$sigma_hat, m_krls$sigma_hat))

    cell_done <- cell_done + 1L
  }

  elapsed <- proc.time()[["elapsed"]] - wall_start
  if (elapsed > WALL_CAP_S - 15) break
}

## -----------------------------------------------------------------------
## Save RDS
## -----------------------------------------------------------------------
results_df <- do.call(rbind, all_rows)
saveRDS(results_df, RDS_PATH)
cat(sprintf("\nFull results: %d rows -> %s\n", nrow(results_df), RDS_PATH))

## -----------------------------------------------------------------------
## Build paired dataset and summary
## -----------------------------------------------------------------------
rr_df   <- results_df[results_df$pkg == "rr",   ]
krls_df <- results_df[results_df$pkg == "krls", ]

paired <- merge(
  rr_df[, c("dgp","rep","seed","test_MSE","train_MSE","sigma_hat","lambda_hat","wall_s","n_warnings","n_errors")],
  krls_df[, c("dgp","rep","seed","test_MSE","train_MSE","sigma_hat","lambda_hat","wall_s","n_warnings","n_errors")],
  by = c("dgp","rep","seed"), suffixes = c("_rr","_krls")
)
paired$mse_ratio_rr_over_krls <- paired$test_MSE_rr / paired$test_MSE_krls

## -----------------------------------------------------------------------
## Summary per DGP
## -----------------------------------------------------------------------
dgp_names_done <- unique(paired$dgp)

summary_rows <- lapply(dgp_names_done, function(dgp_nm) {
  sub  <- paired[paired$dgp == dgp_nm, ]
  n_rep <- nrow(sub)

  mean_rr    <- mean(sub$test_MSE_rr,   na.rm = TRUE)
  mean_krls  <- mean(sub$test_MSE_krls, na.rm = TRUE)
  mean_ratio <- mean(sub$mse_ratio_rr_over_krls, na.rm = TRUE)
  sd_ratio   <- sd(sub$mse_ratio_rr_over_krls,   na.rm = TRUE)

  log_ratio  <- log(sub$mse_ratio_rr_over_krls)
  if (n_rep >= 2 && !all(is.na(log_ratio))) {
    tt    <- t.test(log_ratio, mu = 0)
    pval  <- tt$p.value
    ci_lo <- exp(tt$conf.int[1])
    ci_hi <- exp(tt$conf.int[2])
  } else {
    pval  <- NA_real_; ci_lo <- NA_real_; ci_hi <- NA_real_
  }

  winner <- if (is.na(pval) || pval > 0.10) "tie" else
            if (mean_ratio < 1) "rr" else "krls"
  if (!is.na(mean_ratio) && abs(mean_ratio - 1) < 0.02) winner <- "tie"

  data.frame(
    dgp = dgp_nm, n_rep = n_rep,
    mean_test_MSE_rr = mean_rr, mean_test_MSE_krls = mean_krls,
    mean_ratio = mean_ratio, sd_ratio = sd_ratio,
    ci_lo = ci_lo, ci_hi = ci_hi,
    pval_log_ratio = pval, winner = winner,
    stringsAsFactors = FALSE
  )
})

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, SUM_CSV, row.names = FALSE)
cat(sprintf("Summary written: %s\n", SUM_CSV))

## -----------------------------------------------------------------------
## Delta vs iter-0
## -----------------------------------------------------------------------
## ratio_delta = iter2 mean_ratio - iter0 mean_ratio (negative = improved)
summary_df$ratio_delta_vs_iter0 <- NA_real_
for (i in seq_len(nrow(summary_df))) {
  dgp_nm <- summary_df$dgp[i]
  row0   <- iter0_sum[iter0_sum$dgp == dgp_nm, ]
  if (nrow(row0) == 1) {
    summary_df$ratio_delta_vs_iter0[i] <-
      summary_df$mean_ratio[i] - row0$mean_ratio
  }
}
# Re-write CSV with delta column
write.csv(summary_df, SUM_CSV, row.names = FALSE)

## -----------------------------------------------------------------------
## Plots: iter-2 paired scatter (per-DGP)
## -----------------------------------------------------------------------
png(file.path(PLOTS_DIR, "paired-scatter-iter2.png"),
    width = 1800, height = 1200, res = 150)
n_dgp <- length(dgp_names_done)
n_col <- ceiling(sqrt(n_dgp))
n_row <- ceiling(n_dgp / n_col)
par(mfrow = c(n_row, n_col), mar = c(3, 3, 2, 1), mgp = c(2, 0.7, 0))
for (dgp_nm in dgp_names_done) {
  sub <- paired[paired$dgp == dgp_nm, ]
  rng <- range(c(sub$test_MSE_rr, sub$test_MSE_krls), na.rm = TRUE)
  plot(sub$test_MSE_krls, sub$test_MSE_rr,
       xlab = "KRLS test MSE", ylab = "rr test MSE",
       main = dgp_nm, xlim = rng, ylim = rng,
       pch = 16, col = "steelblue", cex = 1.2)
  abline(0, 1, lty = 2, col = "grey50")
  sm <- summary_df[summary_df$dgp == dgp_nm, ]
  legend("topleft", bty = "n",
         legend = sprintf("ratio=%.3f\n%s", sm$mean_ratio, sm$winner),
         cex = 0.8)
}
dev.off()

## -----------------------------------------------------------------------
## Counts
## -----------------------------------------------------------------------
lose_dgps <- summary_df[!is.na(summary_df$winner) & summary_df$winner == "krls", ]
win_dgps  <- summary_df[!is.na(summary_df$winner) & summary_df$winner == "rr",   ]
tie_dgps  <- summary_df[!is.na(summary_df$winner) & summary_df$winner == "tie",  ]

if (nrow(lose_dgps) > 0) {
  worst_idx   <- which.max(lose_dgps$mean_ratio)
  worst_dgp   <- lose_dgps$dgp[worst_idx]
  worst_ratio <- lose_dgps$mean_ratio[worst_idx]
} else {
  worst_dgp   <- "(none)"
  worst_ratio <- NA_real_
}

## -----------------------------------------------------------------------
## Headline table (with delta col)
## -----------------------------------------------------------------------
## order by DGP name
sum_ord <- summary_df[order(summary_df$dgp), ]

headline_lines <- apply(sum_ord, 1, function(r) {
  delta <- as.numeric(r["ratio_delta_vs_iter0"])
  delta_str <- if (is.na(delta)) "  N/A  " else sprintf("%+.3f", delta)
  sprintf("| %-18s | %5s | %.4f | %.4f | %.3f [%.3f, %.3f] | %s | %s |",
          r["dgp"], r["winner"],
          as.numeric(r["mean_test_MSE_rr"]),
          as.numeric(r["mean_test_MSE_krls"]),
          as.numeric(r["mean_ratio"]),
          as.numeric(r["ci_lo"]),
          as.numeric(r["ci_hi"]),
          formatC(as.numeric(r["pval_log_ratio"]), digits = 3, format = "f"),
          delta_str)
})

## iter-0 headline for comparison
iter0_lines <- apply(iter0_sum[order(iter0_sum$dgp), ], 1, function(r) {
  sprintf("| %-18s | %5s | %.4f | %.4f | %.3f |",
          r["dgp"], r["winner"],
          as.numeric(r["mean_test_MSE_rr"]),
          as.numeric(r["mean_test_MSE_krls"]),
          as.numeric(r["mean_ratio"]))
})

## loss diagnostics
loss_hyp_lines <- if (nrow(lose_dgps) > 0) {
  sapply(seq_len(nrow(lose_dgps)), function(i) {
    dgp_nm  <- lose_dgps$dgp[i]
    ratio   <- lose_dgps$mean_ratio[i]
    sub_rr   <- rr_df[rr_df$dgp == dgp_nm, ]
    sub_krls <- krls_df[krls_df$dgp == dgp_nm, ]
    ms_rr    <- mean(sub_rr$sigma_hat,   na.rm = TRUE)
    ms_krls  <- mean(sub_krls$sigma_hat, na.rm = TRUE)

    hyp <- if (ms_rr < ms_krls * 0.8) {
      sprintf("rr sigma (%.2f) << KRLS sigma (%.2f): kernel still too narrow.", ms_rr, ms_krls)
    } else if (ms_rr > ms_krls * 1.5) {
      sprintf("rr sigma (%.2f) >> KRLS sigma (%.2f): kernel too wide.", ms_rr, ms_krls)
    } else {
      sprintf("rr sigma (%.2f) similar to KRLS sigma (%.2f); loss driven by lambda or noise DGP structure.", ms_rr, ms_krls)
    }
    sprintf("- **%s** (ratio=%.3f): %s", dgp_nm, ratio, hyp)
  })
} else {
  c("(none)")
}

## PASS/FAIL verdict
n_wins   <- nrow(win_dgps)
n_losses <- nrow(lose_dgps)
n_ties   <- nrow(tie_dgps)
verdict  <- if (n_wins + n_ties >= 11 && n_losses <= 4) "PASS" else "FAIL"

total_wall <- proc.time()[["elapsed"]] - wall_start

report_md <- c(
  "# KRLS Head-to-Head Iter-2: roadrunner (geomean_p sigma anchor) vs KRLS (CRAN)",
  "",
  sprintf("**Run**: REQ-20260518-003 iter-2  **Date**: 2026-05-18"),
  sprintf("**Patch**: `.krls_sigma_anchor` uses `sqrt(median(d2) * ncol(Xs))` (geomean_p)"),
  sprintf("**Config**: n=%d, p=%d, n_test=%d, R=%d reps, sequential (nthreads=1)",
          N_TRAIN, P, N_TEST, R_REPS),
  sprintf("**DGPs completed**: %d / %d", length(dgp_names_done), length(DGPS)),
  sprintf("**Total wall-clock**: %.1f s (cap: %.0f s)", total_wall, WALL_CAP_S),
  "",
  sprintf("## Verdict: %s", verdict),
  "",
  sprintf("roadrunner wins **%d / %d** DGPs  (losses: %d, ties: %d)",
          n_wins, length(dgp_names_done), n_losses, n_ties),
  sprintf("Criterion: PASS if wins + ties >= 11 and losses <= 4"),
  "",
  if (worst_dgp != "(none)") {
    sprintf("**Top remaining loss**: %s  mean_ratio=%.3f", worst_dgp, worst_ratio)
  } else {
    "**Top remaining loss**: (none)"
  },
  "",
  "---",
  "",
  "## Iter-2 Headline Table",
  "",
  "| DGP               | Winner | rr_MSE | krls_MSE | ratio [95% CI]         | p-val | delta_vs_iter0 |",
  "|-------------------|--------|--------|----------|------------------------|-------|----------------|",
  paste(headline_lines, collapse = "\n"),
  "",
  paste0("delta_vs_iter0 = iter2_mean_ratio - iter0_mean_ratio. ",
         "Negative = improved (rr relatively better than iter-0)."),
  "",
  "---",
  "",
  "## Iter-0 Reference Table (for comparison)",
  "",
  "| DGP               | Winner | rr_MSE | krls_MSE | ratio |",
  "|-------------------|--------|--------|----------|-------|",
  paste(iter0_lines, collapse = "\n"),
  "",
  "Iter-0 score: rr wins=4, KRLS wins=8, ties=3",
  "",
  "---",
  "",
  "## LOSE DGPs (roadrunner still worse after patch)",
  "",
  paste(loss_hyp_lines, collapse = "\n"),
  "",
  "---",
  "",
  "## Sigma Comparison (iter-2)",
  "",
  "| DGP               | mean_sigma_rr | mean_sigma_krls |",
  "|-------------------|--------------|-----------------|",
  paste(sapply(dgp_names_done, function(g) {
    s_rr   <- mean(rr_df$sigma_hat[rr_df$dgp == g],     na.rm = TRUE)
    s_krls <- mean(krls_df$sigma_hat[krls_df$dgp == g], na.rm = TRUE)
    sprintf("| %-18s | %13.3f | %15.3f |", g, s_rr, s_krls)
  }), collapse = "\n"),
  "",
  "---",
  "",
  "## Wall-Clock",
  "",
  sprintf("Total elapsed: %.1f s", total_wall),
  "",
  "| DGP               | mean_wall_rr_s | mean_wall_krls_s |",
  "|-------------------|---------------|-----------------|",
  paste(sapply(dgp_names_done, function(g) {
    w_rr   <- mean(rr_df$wall_s[rr_df$dgp == g],     na.rm = TRUE)
    w_krls <- mean(krls_df$wall_s[krls_df$dgp == g], na.rm = TRUE)
    sprintf("| %-18s | %14.2f | %16.2f |", g, w_rr, w_krls)
  }), collapse = "\n"),
  "",
  "---",
  "",
  "## Warnings and Errors",
  "",
  sprintf("Total rr warnings: %d  errors: %d",
          sum(rr_df$n_warnings, na.rm = TRUE), sum(rr_df$n_errors, na.rm = TRUE)),
  sprintf("Total krls warnings: %d  errors: %d",
          sum(krls_df$n_warnings, na.rm = TRUE), sum(krls_df$n_errors, na.rm = TRUE)),
  "",
  "---",
  "",
  "## Files",
  sprintf("- `%s`", CSV_PATH),
  sprintf("- `%s`", RDS_PATH),
  sprintf("- `%s`", SUM_CSV),
  sprintf("- `%s/paired-scatter-iter2.png`", PLOTS_DIR)
)

writeLines(report_md, REPORT_MD)
cat(sprintf("Report written: %s\n", REPORT_MD))

## -----------------------------------------------------------------------
## Console summary
## -----------------------------------------------------------------------
cat("\n========== ITER-2 FINAL SUMMARY ==========\n")
cat(sprintf("VERDICT: %s\n", verdict))
cat(sprintf("DGPs completed: %d / %d\n", length(dgp_names_done), length(DGPS)))
cat(sprintf("rr wins: %d  KRLS wins: %d  ties: %d\n",
            n_wins, n_losses, n_ties))
cat(sprintf("Worst loss: %s  ratio=%.3f\n",
            worst_dgp, if (is.na(worst_ratio)) 0 else worst_ratio))
cat(sprintf("Total wall-clock: %.1f s\n", total_wall))
cat("\nPer-DGP summary (iter-2 vs iter-0):\n")
print_df <- summary_df[order(summary_df$mean_ratio, decreasing = TRUE),
                        c("dgp","winner","mean_ratio","ci_lo","ci_hi","ratio_delta_vs_iter0")]
print(print_df, row.names = FALSE, digits = 4)
cat("==========================================\n")
