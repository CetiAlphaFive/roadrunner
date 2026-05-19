## krls-headtohead-REQ-20260518-003.R
##
## Head-to-head: roadrunner::krls() vs KRLS::krls() on 15 DGPs.
##
## Both packages called in DEFAULT mode (no autotune, no manual sigma/lambda).
##   roadrunner: roadrunner::krls(X, y)  -- median-heuristic sigma (v0.0.0.9040)
##   reference:  KRLS::krls(X, y, print.level=0) -- sigma = ncol(X)
##
## Configuration:
##   n = 500, p = 10, n_test = 1000, R = 5 reps per DGP per package
##   Sequential (mc.cores = 1), nthreads = 1
##   Paired: same (X, y) seed per rep
##   Wall-clock cap: 5 min
##
## REQ-20260518-003 -- simulator: statsclaw pipeline

suppressPackageStartupMessages({
  ## Use load_all so we pick up the source-tree version (including krls).
  ## Falls back to library(roadrunner) if devtools is unavailable.
  if (requireNamespace("devtools", quietly = TRUE)) {
    pkg_root <- "/home/jack/Dropbox/roadrunner"
    devtools::load_all(pkg_root, quiet = TRUE)
  } else {
    library(roadrunner)
  }
  library(KRLS)
  library(RcppParallel)
})

RcppParallel::setThreadOptions(numThreads = 1L)

## -----------------------------------------------------------------------
## Global config
## -----------------------------------------------------------------------
RUN_DIR    <- "/home/jack/.claude/plugins/data/statsclaw-statsclaw/workspace/roadrunner/runs/REQ-20260518-003"
PLOTS_DIR  <- file.path(RUN_DIR, "headtohead-plots")
CSV_PATH   <- file.path(RUN_DIR, "headtohead-results.csv")
RDS_PATH   <- file.path(RUN_DIR, "headtohead-results.rds")
SUM_CSV    <- file.path(RUN_DIR, "headtohead-summary.csv")
REPORT_MD  <- file.path(RUN_DIR, "headtohead-report.md")

N_TRAIN    <- 500L
N_TEST     <- 1000L
P          <- 10L
R_REPS     <- 5L
WALL_CAP_S <- 300      # 5 min hard cap
MASTER_SEED <- 20260518L

dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

wall_start <- proc.time()[["elapsed"]]

## -----------------------------------------------------------------------
## DGP definitions
## Each function takes (n, p, seed) and returns list(X, y, f_true)
## where f_true is the noiseless signal (used for SNR calibration).
## sigma_eps = sd(f_true) / 2  =>  SNR ~ 4.
## X is drawn from N(0,1)^p unless specified.
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
  ## Use Friedman's suggested scaling for x: x1~U(0,100), x2~U(40pi,560pi),
  ## x3~U(0,1), x4~U(1,11).  Scale to avoid numeric blow-up.
  x1 <- runif(n, 0, 100)
  x2 <- runif(n, 40 * pi, 560 * pi)
  x3 <- runif(n, 0, 1)
  x4 <- runif(n, 1, 11)
  ## Fill remaining cols with U(0,1)
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
  ## calibrate noise so SNR ~ 4 at typical X values
  eps_sd_base <- sd(f) / (2 * mean(eps_scale))
  y <- f + rnorm(n, sd = eps_sd_base) * eps_scale
  list(X = X, y = y, f = f, eps_sd = eps_sd_base)
}

## -----------------------------------------------------------------------
## DGP catalog
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
## Predict helper — unified for both packages
## -----------------------------------------------------------------------
predict_vec <- function(fit, X_test) {
  pred <- predict(fit, newdata = X_test)
  ## roadrunner returns list with $fit; KRLS predict.krls also returns list
  if (is.list(pred)) {
    if (!is.null(pred$fit)) return(as.numeric(pred$fit))
    if (!is.null(pred$yfitted)) return(as.numeric(pred$yfitted))
  }
  as.numeric(pred)
}

mse <- function(a, b) mean((a - b)^2)
r2  <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)

## -----------------------------------------------------------------------
## Single-fit runner for one package on one rep
## Returns named list of metrics or NA on error.
## -----------------------------------------------------------------------
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

    y_fit <- as.numeric(fit$fitted)
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

  list(
    sigma_hat  = sigma_hat,
    lambda_hat = lambda_hat,
    train_MSE  = train_MSE,
    test_MSE   = test_MSE,
    train_R2   = train_R2,
    test_R2    = test_R2,
    wall_s     = wall_s,
    n_warnings = n_warn,
    n_errors   = n_err
  )
}

## -----------------------------------------------------------------------
## Row builder
## -----------------------------------------------------------------------
make_row <- function(dgp_name, rep, seed, pkg, metrics) {
  data.frame(
    dgp        = dgp_name,
    rep        = rep,
    seed       = seed,
    pkg        = pkg,
    sigma_hat  = metrics$sigma_hat,
    lambda_hat = metrics$lambda_hat,
    train_MSE  = metrics$train_MSE,
    test_MSE   = metrics$test_MSE,
    train_R2   = metrics$train_R2,
    test_R2    = metrics$test_R2,
    overfit_ratio = metrics$test_MSE / metrics$train_MSE,
    wall_s     = metrics$wall_s,
    n_warnings = metrics$n_warnings,
    n_errors   = metrics$n_errors,
    stringsAsFactors = FALSE
  )
}

## -----------------------------------------------------------------------
## Append-to-CSV helper (writes header on first call)
## -----------------------------------------------------------------------
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
cat("=== SMOKE TEST ===\n")
smoke_seed <- MASTER_SEED + 999999L
set.seed(smoke_seed)
X_sm  <- matrix(rnorm(200 * P), 200, P)
y_sm  <- X_sm[, 1] + rnorm(200, sd = 0.5)
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
total_cells <- length(DGPS) * R_REPS
cell_done  <- 0L

for (dgp_i in seq_along(DGPS)) {
  dgp_name <- DGPS[[dgp_i]]$name
  dgp_fn   <- DGPS[[dgp_i]]$fn

  cat(sprintf("[%d/%d] DGP: %s\n", dgp_i, length(DGPS), dgp_name))

  for (r in seq_len(R_REPS)) {
    # check wall-clock budget
    elapsed <- proc.time()[["elapsed"]] - wall_start
    if (elapsed > WALL_CAP_S - 15) {
      cat(sprintf("  WALL-CLOCK CAP approaching (%.0fs elapsed). Aborting remaining DGPs.\n",
                  elapsed))
      break
    }

    ## Derive seeds deterministically
    train_seed <- MASTER_SEED + (dgp_i - 1L) * R_REPS * 3L + (r - 1L) * 3L + 1L
    test_seed  <- MASTER_SEED + (dgp_i - 1L) * R_REPS * 3L + (r - 1L) * 3L + 2L

    ## Generate TRAINING data
    d_train <- dgp_fn(N_TRAIN, P, train_seed)
    X_train <- d_train$X
    y_train <- d_train$y

    ## Generate TEST data (same DGP function, different seed)
    d_test  <- dgp_fn(N_TEST, P, test_seed)
    X_test  <- d_test$X
    y_test  <- d_test$y

    ## ---- roadrunner ----
    m_rr   <- run_one("rr",   X_train, y_train, X_test, y_test)
    row_rr <- make_row(dgp_name, r, train_seed, "rr", m_rr)
    row_idx <- row_idx + 1L
    all_rows[[row_idx]] <- row_rr
    append_row(row_rr)

    ## ---- KRLS ----
    m_krls   <- run_one("krls", X_train, y_train, X_test, y_test)
    row_krls <- make_row(dgp_name, r, train_seed, "krls", m_krls)
    row_idx <- row_idx + 1L
    all_rows[[row_idx]] <- row_krls
    append_row(row_krls)

    ## paired ratio
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

## match by (dgp, rep, seed)
paired <- merge(
  rr_df[, c("dgp","rep","seed","test_MSE","train_MSE","sigma_hat","lambda_hat","wall_s","n_warnings","n_errors")],
  krls_df[, c("dgp","rep","seed","test_MSE","train_MSE","sigma_hat","lambda_hat","wall_s","n_warnings","n_errors")],
  by = c("dgp","rep","seed"),
  suffixes = c("_rr","_krls")
)
paired$mse_ratio_rr_over_krls <- paired$test_MSE_rr / paired$test_MSE_krls

## -----------------------------------------------------------------------
## Summary per DGP
## -----------------------------------------------------------------------
dgp_names_done <- unique(paired$dgp)

summary_rows <- lapply(dgp_names_done, function(dgp_nm) {
  sub <- paired[paired$dgp == dgp_nm, ]
  n_rep <- nrow(sub)

  mean_rr   <- mean(sub$test_MSE_rr,   na.rm = TRUE)
  mean_krls <- mean(sub$test_MSE_krls, na.rm = TRUE)
  mean_ratio <- mean(sub$mse_ratio_rr_over_krls, na.rm = TRUE)
  sd_ratio   <- sd(sub$mse_ratio_rr_over_krls,   na.rm = TRUE)

  ## paired t-test on log MSE ratio (better distributional properties)
  log_ratio  <- log(sub$mse_ratio_rr_over_krls)
  if (n_rep >= 2 && !all(is.na(log_ratio))) {
    tt  <- t.test(log_ratio, mu = 0)
    pval <- tt$p.value
    ci_lo <- exp(tt$conf.int[1])
    ci_hi <- exp(tt$conf.int[2])
  } else {
    pval  <- NA_real_
    ci_lo <- NA_real_
    ci_hi <- NA_real_
  }

  ## winner: rr if mean_ratio < 1, krls if > 1, tie if within 2% and non-sig
  winner <- if (is.na(pval) || pval > 0.10) {
    "tie"
  } else if (mean_ratio < 1) {
    "rr"
  } else {
    "krls"
  }
  ## also call tie if ratio within 2% regardless of p-value (trivial diff)
  if (!is.na(mean_ratio) && abs(mean_ratio - 1) < 0.02) winner <- "tie"

  data.frame(
    dgp         = dgp_nm,
    n_rep       = n_rep,
    mean_test_MSE_rr   = mean_rr,
    mean_test_MSE_krls = mean_krls,
    mean_ratio   = mean_ratio,
    sd_ratio     = sd_ratio,
    ci_lo        = ci_lo,
    ci_hi        = ci_hi,
    pval_log_ratio = pval,
    winner       = winner,
    stringsAsFactors = FALSE
  )
})

summary_df <- do.call(rbind, summary_rows)
write.csv(summary_df, SUM_CSV, row.names = FALSE)
cat(sprintf("Summary written: %s\n", SUM_CSV))

## -----------------------------------------------------------------------
## Diagnostic plots
## -----------------------------------------------------------------------

## 1. Per-DGP paired scatter: test_MSE_rr vs test_MSE_krls
png(file.path(PLOTS_DIR, "paired-scatter.png"), width = 1800, height = 1200, res = 150)
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

## 2. Bar chart of mean ratio per DGP
png(file.path(PLOTS_DIR, "winner-bar.png"), width = 1600, height = 800, res = 150)
par(mar = c(7, 5, 3, 2))
cols <- ifelse(summary_df$mean_ratio < 0.98, "steelblue",
        ifelse(summary_df$mean_ratio > 1.02, "firebrick", "grey70"))
bp <- barplot(summary_df$mean_ratio,
              names.arg = summary_df$dgp,
              las = 2,
              col = cols,
              ylab = "Mean test_MSE ratio (rr/krls)",
              main = "roadrunner vs KRLS: test MSE ratio per DGP\n(blue < 1 = rr wins, red > 1 = rr loses)",
              ylim = c(0, max(summary_df$mean_ratio, na.rm = TRUE) * 1.15))
abline(h = 1, lty = 2, col = "black", lwd = 2)
## add error bars (sd / sqrt(n))
se_ratio <- summary_df$sd_ratio / sqrt(summary_df$n_rep)
arrows(bp, summary_df$mean_ratio - se_ratio,
       bp, summary_df$mean_ratio + se_ratio,
       angle = 90, code = 3, length = 0.05)
dev.off()

## 3. Sigma comparison plot
png(file.path(PLOTS_DIR, "sigma-comparison.png"), width = 1600, height = 800, res = 150)
par(mar = c(7, 5, 3, 2))
mean_sigma_rr   <- tapply(rr_df$sigma_hat,   rr_df$dgp,   mean, na.rm = TRUE)
mean_sigma_krls <- tapply(krls_df$sigma_hat, krls_df$dgp, mean, na.rm = TRUE)
dgp_ord <- names(mean_sigma_rr)[order(names(mean_sigma_rr))]
xx <- seq_along(dgp_ord)
plot(xx, mean_sigma_rr[dgp_ord], type = "b", pch = 16, col = "steelblue",
     xlab = "", ylab = "Mean sigma_hat", main = "Default sigma: rr (blue) vs KRLS (red)",
     xaxt = "n",
     ylim = range(c(mean_sigma_rr, mean_sigma_krls), na.rm = TRUE))
lines(xx, mean_sigma_krls[dgp_ord], type = "b", pch = 17, col = "firebrick")
axis(1, at = xx, labels = dgp_ord, las = 2)
legend("topright", legend = c("roadrunner", "KRLS"), col = c("steelblue","firebrick"),
       pch = c(16,17), lty = 1)
dev.off()

cat(sprintf("Plots written to %s\n", PLOTS_DIR))

## -----------------------------------------------------------------------
## Write report
## -----------------------------------------------------------------------
total_wall <- proc.time()[["elapsed"]] - wall_start

lose_dgps <- summary_df[!is.na(summary_df$winner) & summary_df$winner == "krls", ]
win_dgps  <- summary_df[!is.na(summary_df$winner) & summary_df$winner == "rr",   ]
tie_dgps  <- summary_df[!is.na(summary_df$winner) & summary_df$winner == "tie",  ]

## worst loss: highest mean_ratio
if (nrow(lose_dgps) > 0) {
  worst_idx <- which.max(lose_dgps$mean_ratio)
  worst_dgp <- lose_dgps$dgp[worst_idx]
  worst_ratio <- lose_dgps$mean_ratio[worst_idx]
} else {
  worst_dgp   <- "(none)"
  worst_ratio <- NA_real_
}

## headline table rows
headline_lines <- apply(summary_df[order(summary_df$dgp),], 1, function(r) {
  sprintf("| %-18s | %5s | %.4f | %.4f | %.3f [%.3f, %.3f] | %s |",
          r["dgp"], r["winner"],
          as.numeric(r["mean_test_MSE_rr"]),
          as.numeric(r["mean_test_MSE_krls"]),
          as.numeric(r["mean_ratio"]),
          as.numeric(r["ci_lo"]),
          as.numeric(r["ci_hi"]),
          formatC(as.numeric(r["pval_log_ratio"]), digits = 3, format = "f"))
})

## diagnostic hypotheses for losses
loss_hyp_lines <- if (nrow(lose_dgps) > 0) {
  sapply(seq_len(nrow(lose_dgps)), function(i) {
    dgp_nm <- lose_dgps$dgp[i]
    ratio  <- lose_dgps$mean_ratio[i]
    sub_rr   <- rr_df[rr_df$dgp == dgp_nm, ]
    sub_krls <- krls_df[krls_df$dgp == dgp_nm, ]
    ms_rr   <- mean(sub_rr$sigma_hat,   na.rm = TRUE)
    ms_krls <- mean(sub_krls$sigma_hat, na.rm = TRUE)

    hyp <- if (ms_rr < ms_krls * 0.8) {
      sprintf("rr sigma (%.2f) << KRLS sigma (%.2f): median heuristic selects too narrow a kernel for this DGP, producing underfitting.", ms_rr, ms_krls)
    } else if (ms_rr > ms_krls * 1.5) {
      sprintf("rr sigma (%.2f) >> KRLS sigma (%.2f): median heuristic selects too wide a kernel, over-smoothing the signal.", ms_rr, ms_krls)
    } else {
      sprintf("rr sigma (%.2f) similar to KRLS sigma (%.2f); loss likely driven by lambda selection (rr LOO with tighter tol vs KRLS).", ms_rr, ms_krls)
    }
    sprintf("- **%s** (ratio=%.3f): %s", dgp_nm, ratio, hyp)
  })
} else {
  "(none — roadrunner wins or ties all DGPs)"
}

report_md <- c(
  "# KRLS Head-to-Head: roadrunner (v0.0.0.9040) vs KRLS (CRAN)",
  "",
  sprintf("**Run**: REQ-20260518-003  **Date**: 2026-05-18"),
  sprintf("**Config**: n=%d, p=%d, n_test=%d, R=%d reps, sequential (nthreads=1)",
          N_TRAIN, P, N_TEST, R_REPS),
  sprintf("**DGPs completed**: %d / %d", length(dgp_names_done), length(DGPS)),
  sprintf("**Total wall-clock**: %.1f s (cap: %.0f s)", total_wall, WALL_CAP_S),
  "",
  "---",
  "",
  "## Headline Table",
  "",
  "| DGP               | Winner | rr_MSE | krls_MSE | ratio [95% CI]         | p-val |",
  "|-------------------|--------|--------|----------|------------------------|-------|",
  paste(headline_lines, collapse = "\n"),
  "",
  sprintf("**Winner counts**: roadrunner wins=%d, KRLS wins=%d, ties=%d",
          nrow(win_dgps), nrow(lose_dgps), nrow(tie_dgps)),
  "",
  "---",
  "",
  "## LOSE DGPs (roadrunner MSE > KRLS MSE, CI excludes 1 or strong trend)",
  "",
  paste(loss_hyp_lines, collapse = "\n"),
  "",
  "---",
  "",
  "## Worst Loss",
  "",
  sprintf("**DGP**: %s  **ratio**: %.3f", worst_dgp,
          if (is.na(worst_ratio)) 0 else worst_ratio),
  "",
  "---",
  "",
  "## Sigma Comparison",
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
  "## Wall-Clock Summary",
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
          sum(rr_df$n_warnings, na.rm=TRUE), sum(rr_df$n_errors, na.rm=TRUE)),
  sprintf("Total krls warnings: %d  errors: %d",
          sum(krls_df$n_warnings, na.rm=TRUE), sum(krls_df$n_errors, na.rm=TRUE)),
  "",
  "---",
  "",
  "## Files",
  sprintf("- `%s` (per-replicate results, both packages flat-stacked)", CSV_PATH),
  sprintf("- `%s` (RDS)", RDS_PATH),
  sprintf("- `%s` (per-DGP summary)", SUM_CSV),
  sprintf("- `%s/paired-scatter.png`", PLOTS_DIR),
  sprintf("- `%s/winner-bar.png`", PLOTS_DIR),
  sprintf("- `%s/sigma-comparison.png`", PLOTS_DIR)
)

writeLines(report_md, REPORT_MD)
cat(sprintf("Report written: %s\n", REPORT_MD))

## -----------------------------------------------------------------------
## Console summary
## -----------------------------------------------------------------------
cat("\n========== FINAL SUMMARY ==========\n")
cat(sprintf("DGPs completed: %d / %d\n", length(dgp_names_done), length(DGPS)))
cat(sprintf("rr wins: %d  KRLS wins (rr LOSES): %d  ties: %d\n",
            nrow(win_dgps), nrow(lose_dgps), nrow(tie_dgps)))
cat(sprintf("Worst loss: %s  ratio=%.3f\n",
            worst_dgp, if(is.na(worst_ratio)) 0 else worst_ratio))
cat(sprintf("Total wall-clock: %.1f s\n", total_wall))
cat("\nPer-DGP summary:\n")
print(summary_df[order(summary_df$mean_ratio, decreasing = TRUE),
                 c("dgp","winner","mean_ratio","ci_lo","ci_hi","mean_test_MSE_rr","mean_test_MSE_krls")],
      row.names = FALSE, digits = 4)
cat("=====================================\n")
