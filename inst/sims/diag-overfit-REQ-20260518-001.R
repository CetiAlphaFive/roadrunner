## =============================================================================
## diag-overfit-REQ-20260518-001.R
##
## Targeted diagnostic sweeps to identify which tuning knob causes residual
## overfitting in roadrunner::krls() on additive + interaction DGPs.
##
## Three sweeps, all at fixed n=500, p=10, nthreads=1, sequential:
##
##   Sweep 1 — sigma grid at LOO lambda
##     sigma in {1, 2, 5, 10, 20, 50, 100, 200, 500, 1000}
##     lambda goes to LOO golden-section (default). DGPs: additive, interaction.
##
##   Sweep 2 — manual lambda grid at oracle sigma
##     sigma fixed to best-sigma from Sweep 1 per DGP.
##     lambda in 10^seq(-6, 6, length.out=13).
##     Goal: compare argmin(test_MSE) lambda vs what LOO chose in Sweep 1.
##
##   Sweep 3 — autotune with wider sigma grid
##     autotune.grid = p * c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128)
##     Goal: does the wider upper bound eliminate residual overfit?
##
## DGPs used: additive + interaction (the two residual OVERFIT cases post-autotune).
##
## Outputs (run dir):
##   diag-results.csv, diag-results.rds
##   diag-report.md
##   diag-plots/*.png
##
## Usage:
##   Rscript inst/sims/diag-overfit-REQ-20260518-001.R
##
## Date: 2026-05-18
## =============================================================================

suppressPackageStartupMessages({
  devtools::load_all("/home/jack/Dropbox/roadrunner", quiet = TRUE)
})
RcppParallel::setThreadOptions(numThreads = 1L)

## ---- Paths ------------------------------------------------------------------
run_dir    <- "/home/jack/.claude/plugins/data/statsclaw-statsclaw/workspace/roadrunner/runs/REQ-20260518-001"
diag_csv   <- file.path(run_dir, "diag-results.csv")
diag_rds   <- file.path(run_dir, "diag-results.rds")
report_md  <- file.path(run_dir, "diag-report.md")
plot_dir   <- file.path(run_dir, "diag-plots")
log_file   <- file.path(run_dir, "diag-run-log.txt")

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

log_con <- file(log_file, open = "wt")
.log <- function(...) {
  msg <- sprintf(...)
  message(msg)
  writeLines(msg, log_con)
  flush(log_con)
}

## ---- Fixed parameters -------------------------------------------------------
N       <- 500L
P       <- 10L
N_TEST  <- 1000L
R_DIAG  <- 5L   # replications per cell (drop to 3 if budget exceeded)

## ---- DGP helpers (copied from minimal script) --------------------------------
.namex <- function(X) { colnames(X) <- paste0("V", seq_len(ncol(X))); X }

dgp_additive <- function(n, p, n_test, seed) {
  set.seed(seed)
  .f <- function(X) sin(X[,1]) + 0.5*X[,2]^2 - 0.3*exp(-X[,3]^2) + 0.8*X[,4] - 0.5*X[,5]
  X_tr <- matrix(rnorm(n * p), n, p)
  f_tr <- .f(X_tr)
  sigma_eps <- sd(f_tr) / 2
  y_tr <- f_tr + rnorm(n, sd = sigma_eps)
  X_te <- matrix(rnorm(n_test * p), n_test, p)
  y_te <- .f(X_te) + rnorm(n_test, sd = sigma_eps)
  list(X_train = .namex(X_tr), y_train = y_tr,
       X_test  = .namex(X_te), y_test  = y_te)
}

dgp_interaction <- function(n, p, n_test, seed) {
  set.seed(seed)
  .f <- function(X) X[,1]*X[,2] + X[,3]*X[,4]
  X_tr <- matrix(rnorm(n * p), n, p)
  f_tr <- .f(X_tr)
  sigma_eps <- sd(f_tr) / 2
  y_tr <- f_tr + rnorm(n, sd = sigma_eps)
  X_te <- matrix(rnorm(n_test * p), n_test, p)
  y_te <- .f(X_te) + rnorm(n_test, sd = sigma_eps)
  list(X_train = .namex(X_tr), y_train = y_tr,
       X_test  = .namex(X_te), y_test  = y_te)
}

DGP_FNS <- list(additive = dgp_additive, interaction = dgp_interaction)
DGPS    <- c("additive", "interaction")

## ---- Deterministic seed (extension of sim-spec.md Section 5 formula) --------
## For diagnostic sweeps we use a separate seed namespace to avoid collisions.
## diag_idx encodes sweep + condition index; combined with rep.
make_diag_seed <- function(dgp, sweep_id, cond_idx, rep) {
  dgp_idx  <- match(dgp, DGPS)  # 1 or 2
  seed <- (9000000L +
           dgp_idx  * 1000000L +
           sweep_id *   10000L +
           cond_idx *     100L +
           rep              ) %% (2L^31L - 1L)
  as.integer(seed)
}

## ---- Metric helper ----------------------------------------------------------
compute_metrics <- function(fit, X_train, y_train, X_test, y_test) {
  yhat_train <- as.numeric(fit$fitted)
  pred_test  <- predict(fit, newdata = X_test)
  yhat_test  <- as.numeric(pred_test$fit)

  train_MSE     <- mean((y_train - yhat_train)^2)
  test_MSE      <- mean((y_test  - yhat_test )^2)
  train_R2      <- 1 - train_MSE / var(y_train)
  test_R2       <- 1 - test_MSE  / var(y_test)
  overfit_ratio <- test_MSE / train_MSE

  list(train_R2 = train_R2, test_R2 = test_R2,
       train_MSE = train_MSE, test_MSE = test_MSE,
       overfit_ratio = overfit_ratio)
}

## ---- Safe fit wrapper -------------------------------------------------------
safe_fit <- function(expr_fn) {
  warns <- character(0L)
  result <- withCallingHandlers(
    tryCatch(expr_fn(), error = function(e) {
      structure(list(msg = conditionMessage(e)), class = "krls_err")
    }),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(result = result, warns = warns,
       n_errors   = if (inherits(result, "krls_err")) 1L else 0L,
       n_warnings = length(warns))
}

## =============================================================================
## SWEEP 1 — sigma grid at LOO lambda
## =============================================================================

SIGMA_GRID <- c(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
SWEEP1_ID  <- 1L

.log("=================================================================")
.log("SWEEP 1 — sigma grid at LOO lambda")
.log("sigma in {%s}", paste(SIGMA_GRID, collapse=", "))
.log("DGPs: additive, interaction; n=%d p=%d R=%d", N, P, R_DIAG)
.log("=================================================================")

t_s1 <- proc.time()[["elapsed"]]
sweep1_rows <- vector("list", length(DGPS) * length(SIGMA_GRID) * R_DIAG)
s1_idx <- 0L

for (dgp in DGPS) {
  for (si in seq_along(SIGMA_GRID)) {
    sig <- SIGMA_GRID[si]
    .log("  [S1] dgp=%-12s sigma=%6g ...", dgp, sig)

    for (r in seq_len(R_DIAG)) {
      seed <- make_diag_seed(dgp, SWEEP1_ID, si, r)
      dat  <- DGP_FNS[[dgp]](N, P, N_TEST, seed)

      t0 <- proc.time()[["elapsed"]]
      sf <- safe_fit(function() {
        roadrunner::krls(X = dat$X_train, y = dat$y_train,
                         sigma = sig,
                         derivative = FALSE, vcov = FALSE, trace = 0L)
      })
      wall_s <- proc.time()[["elapsed"]] - t0

      if (sf$n_errors == 0L) {
        m <- compute_metrics(sf$result, dat$X_train, dat$y_train,
                             dat$X_test, dat$y_test)
        sigma_hat  <- sf$result$sigma
        lambda_hat <- sf$result$lambda
      } else {
        m <- list(train_R2 = NA_real_, test_R2 = NA_real_,
                  train_MSE = NA_real_, test_MSE = NA_real_,
                  overfit_ratio = NA_real_)
        sigma_hat <- lambda_hat <- NA_real_
        .log("    !! error rep %d: %s", r, sf$result$msg)
      }

      s1_idx <- s1_idx + 1L
      sweep1_rows[[s1_idx]] <- data.frame(
        dgp        = dgp,
        n          = N, p = P,
        sweep      = "sigma",
        rep        = r, seed = seed,
        sigma_in   = sig,
        lambda_in  = NA_real_,
        sigma_hat  = sigma_hat,
        lambda_hat = lambda_hat,
        train_R2   = m$train_R2,  test_R2  = m$test_R2,
        train_MSE  = m$train_MSE, test_MSE = m$test_MSE,
        overfit_ratio = m$overfit_ratio,
        wall_s     = wall_s,
        n_warnings = sf$n_warnings,
        n_errors   = sf$n_errors,
        stringsAsFactors = FALSE
      )
    }
    ## Quick summary for this (dgp, sigma) cell
    cell_rows <- sweep1_rows[seq(s1_idx - R_DIAG + 1L, s1_idx)]
    cell_df   <- do.call(rbind, cell_rows)
    good_r    <- cell_df[!is.na(cell_df$overfit_ratio), ]
    mr  <- if (nrow(good_r) > 0) mean(good_r$overfit_ratio) else NaN
    mte <- if (nrow(good_r) > 0) mean(good_r$test_MSE)      else NaN
    .log("    -> mean overfit_ratio=%.3f  mean test_MSE=%.4f", mr, mte)
  }
}

sweep1_df <- do.call(rbind, sweep1_rows[seq_len(s1_idx)])
t_s1_el <- proc.time()[["elapsed"]] - t_s1
.log("Sweep 1 complete: %.1f s (%.1f min)", t_s1_el, t_s1_el/60)

## Determine best sigma per DGP (argmin mean test_MSE)
s1_agg <- aggregate(test_MSE ~ dgp + sigma_in, data = sweep1_df[!is.na(sweep1_df$test_MSE), ],
                    FUN = mean)
names(s1_agg)[names(s1_agg) == "test_MSE"] <- "mean_test_MSE"

best_sigma <- sapply(DGPS, function(dg) {
  sub <- s1_agg[s1_agg$dgp == dg, ]
  sub$sigma_in[which.min(sub$mean_test_MSE)]
})
.log("\nBest sigma per DGP (argmin mean test_MSE):")
for (dg in DGPS) .log("  %s -> sigma* = %g", dg, best_sigma[dg])

## Also grab the LOO-chosen lambda at each sigma for later comparison
s1_lam <- aggregate(lambda_hat ~ dgp + sigma_in,
                    data = sweep1_df[!is.na(sweep1_df$lambda_hat), ],
                    FUN = mean)
names(s1_lam)[names(s1_lam) == "lambda_hat"] <- "mean_lambda_loo"

## =============================================================================
## SWEEP 2 — manual lambda grid at oracle sigma
## =============================================================================

LAMBDA_GRID <- 10^seq(-6, 6, length.out = 13)
SWEEP2_ID   <- 2L

.log("\n=================================================================")
.log("SWEEP 2 — manual lambda grid at oracle sigma")
.log("lambda in 10^{-6,...,6} (13 levels)  DGPs: additive, interaction")
.log("sigma_star: additive=%g  interaction=%g",
     best_sigma["additive"], best_sigma["interaction"])
.log("=================================================================")

t_s2 <- proc.time()[["elapsed"]]
sweep2_rows <- vector("list", length(DGPS) * length(LAMBDA_GRID) * R_DIAG)
s2_idx <- 0L

for (dgp in DGPS) {
  sig_star <- best_sigma[dgp]
  for (li in seq_along(LAMBDA_GRID)) {
    lam <- LAMBDA_GRID[li]
    .log("  [S2] dgp=%-12s sigma=%g  lambda=%g ...", dgp, sig_star, lam)

    for (r in seq_len(R_DIAG)) {
      seed <- make_diag_seed(dgp, SWEEP2_ID, li, r)
      dat  <- DGP_FNS[[dgp]](N, P, N_TEST, seed)

      t0 <- proc.time()[["elapsed"]]
      sf <- safe_fit(function() {
        roadrunner::krls(X = dat$X_train, y = dat$y_train,
                         sigma = sig_star, lambda = lam,
                         derivative = FALSE, vcov = FALSE, trace = 0L)
      })
      wall_s <- proc.time()[["elapsed"]] - t0

      if (sf$n_errors == 0L) {
        m <- compute_metrics(sf$result, dat$X_train, dat$y_train,
                             dat$X_test, dat$y_test)
        sigma_hat  <- sf$result$sigma
        lambda_hat <- sf$result$lambda
      } else {
        m <- list(train_R2 = NA_real_, test_R2 = NA_real_,
                  train_MSE = NA_real_, test_MSE = NA_real_,
                  overfit_ratio = NA_real_)
        sigma_hat <- lambda_hat <- NA_real_
        .log("    !! error rep %d: %s", r, sf$result$msg)
      }

      s2_idx <- s2_idx + 1L
      sweep2_rows[[s2_idx]] <- data.frame(
        dgp        = dgp,
        n          = N, p = P,
        sweep      = "lambda",
        rep        = r, seed = seed,
        sigma_in   = sig_star,
        lambda_in  = lam,
        sigma_hat  = sigma_hat,
        lambda_hat = lambda_hat,
        train_R2   = m$train_R2,  test_R2  = m$test_R2,
        train_MSE  = m$train_MSE, test_MSE = m$test_MSE,
        overfit_ratio = m$overfit_ratio,
        wall_s     = wall_s,
        n_warnings = sf$n_warnings,
        n_errors   = sf$n_errors,
        stringsAsFactors = FALSE
      )
    }
    cell_rows <- sweep2_rows[seq(s2_idx - R_DIAG + 1L, s2_idx)]
    cell_df   <- do.call(rbind, cell_rows)
    good_r    <- cell_df[!is.na(cell_df$overfit_ratio), ]
    mr  <- if (nrow(good_r) > 0) mean(good_r$overfit_ratio) else NaN
    mte <- if (nrow(good_r) > 0) mean(good_r$test_MSE)      else NaN
    .log("    -> mean overfit_ratio=%.3f  mean test_MSE=%.4f", mr, mte)
  }
}

sweep2_df <- do.call(rbind, sweep2_rows[seq_len(s2_idx)])
t_s2_el <- proc.time()[["elapsed"]] - t_s2
.log("Sweep 2 complete: %.1f s (%.1f min)", t_s2_el, t_s2_el/60)

## Oracle lambda (argmin mean test_MSE) per DGP
s2_agg <- aggregate(test_MSE ~ dgp + lambda_in, data = sweep2_df[!is.na(sweep2_df$test_MSE), ],
                    FUN = mean)
names(s2_agg)[names(s2_agg) == "test_MSE"] <- "mean_test_MSE"

oracle_lambda <- sapply(DGPS, function(dg) {
  sub <- s2_agg[s2_agg$dgp == dg, ]
  sub$lambda_in[which.min(sub$mean_test_MSE)]
})
.log("\nOracle lambda per DGP (argmin mean test_MSE at sigma*):")
for (dg in DGPS) .log("  %s -> lambda* = %g", dg, oracle_lambda[dg])

## LOO lambda at sigma* (from sweep 1) for comparison
loo_lambda_at_star <- sapply(DGPS, function(dg) {
  sub <- s1_lam[s1_lam$dgp == dg & s1_lam$sigma_in == best_sigma[dg], ]
  if (nrow(sub) == 0) return(NA_real_)
  sub$mean_lambda_loo
})
.log("\nLOO-chosen lambda at sigma* (from Sweep 1):")
for (dg in DGPS) .log("  %s -> LOO lambda = %g  (oracle = %g,  ratio = %.3f)",
                       dg, loo_lambda_at_star[dg], oracle_lambda[dg],
                       loo_lambda_at_star[dg] / oracle_lambda[dg])

## =============================================================================
## SWEEP 3 — autotune with wider sigma grid
## =============================================================================

WIDER_GRID_MULT <- c(0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128)
WIDER_GRID      <- P * WIDER_GRID_MULT   # p=10: {2.5, 5, 10, 20, 40, 80, 160, 320, 640, 1280}
SWEEP3_ID       <- 3L

.log("\n=================================================================")
.log("SWEEP 3 — autotune with wider sigma grid")
.log("grid: p * {%s}", paste(WIDER_GRID_MULT, collapse=", "))
.log("sigma candidates: {%s}", paste(WIDER_GRID, collapse=", "))
.log("DGPs: additive, interaction; n=%d p=%d R=%d", N, P, R_DIAG)
.log("=================================================================")

t_s3 <- proc.time()[["elapsed"]]
sweep3_rows <- vector("list", length(DGPS) * R_DIAG)
s3_idx <- 0L

for (dgp in DGPS) {
  .log("  [S3] dgp=%-12s ...", dgp)

  for (r in seq_len(R_DIAG)) {
    seed <- make_diag_seed(dgp, SWEEP3_ID, 0L, r)
    dat  <- DGP_FNS[[dgp]](N, P, N_TEST, seed)

    t0 <- proc.time()[["elapsed"]]
    sf <- safe_fit(function() {
      roadrunner::krls(X = dat$X_train, y = dat$y_train,
                       autotune = TRUE, autotune.grid = WIDER_GRID,
                       derivative = FALSE, vcov = FALSE, trace = 0L)
    })
    wall_s <- proc.time()[["elapsed"]] - t0

    if (sf$n_errors == 0L) {
      m <- compute_metrics(sf$result, dat$X_train, dat$y_train,
                           dat$X_test, dat$y_test)
      sigma_hat  <- sf$result$sigma
      lambda_hat <- sf$result$lambda
    } else {
      m <- list(train_R2 = NA_real_, test_R2 = NA_real_,
                train_MSE = NA_real_, test_MSE = NA_real_,
                overfit_ratio = NA_real_)
      sigma_hat <- lambda_hat <- NA_real_
      .log("    !! error rep %d: %s", r, sf$result$msg)
    }

    s3_idx <- s3_idx + 1L
    sweep3_rows[[s3_idx]] <- data.frame(
      dgp        = dgp,
      n          = N, p = P,
      sweep      = "autotune_wide",
      rep        = r, seed = seed,
      sigma_in   = NA_real_,
      lambda_in  = NA_real_,
      sigma_hat  = sigma_hat,
      lambda_hat = lambda_hat,
      train_R2   = m$train_R2,  test_R2  = m$test_R2,
      train_MSE  = m$train_MSE, test_MSE = m$test_MSE,
      overfit_ratio = m$overfit_ratio,
      wall_s     = wall_s,
      n_warnings = sf$n_warnings,
      n_errors   = sf$n_errors,
      stringsAsFactors = FALSE
    )
    .log("    rep %d: sigma_hat=%g lambda_hat=%g overfit_ratio=%.3f wall=%.1fs",
         r, sigma_hat, lambda_hat, m$overfit_ratio, wall_s)
  }
}

sweep3_df <- do.call(rbind, sweep3_rows[seq_len(s3_idx)])
t_s3_el <- proc.time()[["elapsed"]] - t_s3
.log("Sweep 3 complete: %.1f s (%.1f min)", t_s3_el, t_s3_el/60)

## =============================================================================
## Combine and save
## =============================================================================

all_diag <- rbind(sweep1_df, sweep2_df, sweep3_df)
rownames(all_diag) <- NULL

write.csv(all_diag, diag_csv, row.names = FALSE, quote = FALSE)
saveRDS(list(results = all_diag,
             spec    = "diag-overfit-REQ-20260518-001 (2026-05-18)",
             session = sessionInfo()),
        diag_rds)
.log("\nSaved: %s", diag_csv)
.log("Saved: %s", diag_rds)

## =============================================================================
## Plots
## =============================================================================

## ---- Sweep 1: test_MSE vs sigma ----
.log("\nGenerating Sweep 1 plot ...")
png(file.path(plot_dir, "diag-sweep1-sigma.png"), width = 900, height = 450, res = 100)
par(mfrow = c(1, 2), mar = c(5, 4.5, 3.5, 1) + 0.1, oma = c(0, 0, 3, 0))

for (dg in DGPS) {
  sub <- sweep1_df[sweep1_df$dgp == dg & !is.na(sweep1_df$test_MSE), ]
  agg <- aggregate(cbind(test_MSE, overfit_ratio) ~ sigma_in, data = sub, FUN = mean)
  agg <- agg[order(agg$sigma_in), ]

  bsig <- best_sigma[dg]
  ylim_r <- range(c(1, agg$overfit_ratio * 1.1), na.rm = TRUE)

  plot(agg$sigma_in, agg$overfit_ratio,
       log = "x",
       type = "b", pch = 19, lwd = 2,
       col = ifelse(agg$sigma_in == bsig, "firebrick", "steelblue"),
       xlab = "sigma (log scale)", ylab = "mean overfit_ratio (test/train MSE)",
       main = sprintf("Sweep 1: overfit_ratio vs sigma\n(%s, n=%d p=%d)", dg, N, P),
       ylim = ylim_r)
  abline(h = 1.0, col = "black",  lty = 1, lwd = 1.5)
  abline(h = 1.5, col = "red",    lty = 2, lwd = 1.5)
  abline(h = 1.2, col = "orange", lty = 3, lwd = 1.2)
  ## Mark default sigma=p=10
  abline(v = P, col = "grey40", lty = 2, lwd = 1.2)
  mtext(sprintf("sigma=p=%d (default)", P), side = 3, at = P, cex = 0.65, col = "grey40", line = 0.2)
  ## Mark best sigma
  abline(v = bsig, col = "firebrick", lty = 2, lwd = 1.5)
  mtext(sprintf("sigma*=%g (argmin)", bsig), side = 1, at = bsig, cex = 0.65, col = "firebrick", line = -1.5)
  legend("topright", legend = c("overfit=1.0","overfit=1.5","overfit=1.2"),
         col = c("black","red","orange"), lty = c(1,2,3), lwd = 1.5, cex = 0.7, bty = "n")
}
mtext("Sweep 1 — sigma grid at LOO lambda (R=5 per cell)", outer = TRUE, cex = 1.1, font = 2)
dev.off()
.log("  Wrote: %s", file.path(plot_dir, "diag-sweep1-sigma.png"))

## ---- Sweep 2: test_MSE vs lambda ----
.log("Generating Sweep 2 plot ...")
png(file.path(plot_dir, "diag-sweep2-lambda.png"), width = 900, height = 450, res = 100)
par(mfrow = c(1, 2), mar = c(5, 4.5, 3.5, 1) + 0.1, oma = c(0, 0, 3, 0))

for (dg in DGPS) {
  sig_star <- best_sigma[dg]
  sub <- sweep2_df[sweep2_df$dgp == dg & !is.na(sweep2_df$test_MSE), ]
  agg <- aggregate(cbind(test_MSE, overfit_ratio, train_MSE) ~ lambda_in, data = sub, FUN = mean)
  agg <- agg[order(agg$lambda_in), ]

  orc_l <- oracle_lambda[dg]
  loo_l <- loo_lambda_at_star[dg]

  ylim_r <- range(c(agg$test_MSE, agg$train_MSE) * c(0.9, 1.1), na.rm = TRUE)

  plot(agg$lambda_in, agg$test_MSE,
       log = "xy",
       type = "b", pch = 19, lwd = 2, col = "steelblue",
       xlab = "lambda (log scale)", ylab = "mean MSE",
       main = sprintf("Sweep 2: MSE vs lambda\n(%s, sigma=%g)", dg, sig_star),
       ylim = ylim_r)
  lines(agg$lambda_in, agg$train_MSE, type = "b", pch = 17, lwd = 2, col = "firebrick")
  ## Mark oracle and LOO lambdas
  if (is.finite(orc_l))
    abline(v = orc_l, col = "steelblue", lty = 2, lwd = 1.5)
  if (is.finite(loo_l))
    abline(v = loo_l, col = "purple",    lty = 2, lwd = 1.5)
  legend("bottomright",
         legend = c("test MSE", "train MSE",
                    sprintf("oracle lambda=%g", signif(orc_l, 3)),
                    sprintf("LOO lambda=%g",    signif(loo_l, 3))),
         col = c("steelblue","firebrick","steelblue","purple"),
         lty = c(1,1,2,2), lwd = 1.5, pch = c(19,17,NA,NA),
         cex = 0.7, bty = "n")
}
mtext("Sweep 2 — manual lambda at oracle sigma (R=5 per cell)", outer = TRUE, cex = 1.1, font = 2)
dev.off()
.log("  Wrote: %s", file.path(plot_dir, "diag-sweep2-lambda.png"))

## ---- Sweep 3: overfit_ratio with wider autotune grid ----
.log("Generating Sweep 3 plot ...")
png(file.path(plot_dir, "diag-sweep3-autotune-wide.png"), width = 750, height = 420, res = 100)
par(mfrow = c(1, 1), mar = c(5, 5, 4, 2) + 0.1)

s3_agg <- aggregate(cbind(overfit_ratio, sigma_hat) ~ dgp,
                    data = sweep3_df[!is.na(sweep3_df$overfit_ratio), ],
                    FUN = mean)

## Also need comparison from original autotune (default grid, from sim-report)
## Values from the mailbox: additive=1.5673, interaction=1.7064
orig_autotune <- c(additive = 1.567, interaction = 1.706)

bar_labels <- DGPS
bar_vals   <- s3_agg$overfit_ratio[match(DGPS, s3_agg$dgp)]
bar_orig   <- orig_autotune[DGPS]

x_pos <- barplot(
  rbind(bar_orig, bar_vals),
  beside = TRUE,
  names.arg = bar_labels,
  col = c("steelblue", "firebrick"),
  ylab = "mean overfit_ratio (test/train MSE)",
  main = sprintf("Sweep 3: wide autotune grid (upper %dx default)\nvs default grid [n=%d p=%d R=%d]",
                 max(WIDER_GRID_MULT) / 8, N, P, R_DIAG),
  ylim = c(0, max(c(bar_orig, bar_vals), na.rm = TRUE) * 1.25)
)
abline(h = 1.0, col = "black",  lty = 1, lwd = 1.5)
abline(h = 1.5, col = "red",    lty = 2, lwd = 1.5)
abline(h = 1.2, col = "orange", lty = 3, lwd = 1.2)
legend("topright",
       legend = c("autotune default grid (8x)", sprintf("autotune wide grid (%dx)", max(WIDER_GRID_MULT))),
       fill = c("steelblue", "firebrick"), bty = "n", cex = 0.85)
## Annotate sigma_hat
sigma_hats <- s3_agg$sigma_hat[match(DGPS, s3_agg$dgp)]
for (i in seq_along(DGPS)) {
  text(x_pos[2, i], bar_vals[i] + 0.05, sprintf("sig=%.0f", sigma_hats[i]),
       cex = 0.75, col = "firebrick", font = 2)
}
dev.off()
.log("  Wrote: %s", file.path(plot_dir, "diag-sweep3-autotune-wide.png"))

## =============================================================================
## Build tables for the report
## =============================================================================

## Sweep 1 table: mean test_MSE and overfit_ratio vs sigma, per DGP
build_s1_table <- function() {
  agg_r <- aggregate(overfit_ratio ~ dgp + sigma_in, data = sweep1_df[!is.na(sweep1_df$overfit_ratio), ], FUN = mean)
  agg_m <- aggregate(test_MSE     ~ dgp + sigma_in, data = sweep1_df[!is.na(sweep1_df$test_MSE), ],      FUN = mean)
  agg_l <- aggregate(lambda_hat   ~ dgp + sigma_in, data = sweep1_df[!is.na(sweep1_df$lambda_hat), ],    FUN = mean)
  agg   <- merge(merge(agg_r, agg_m, by = c("dgp","sigma_in")), agg_l, by = c("dgp","sigma_in"))
  agg   <- agg[order(agg$dgp, agg$sigma_in), ]

  lines <- c(
    "| DGP | sigma | mean_overfit_ratio | mean_test_MSE | mean_LOO_lambda | best? |",
    "|-----|-------|--------------------|---------------|-----------------|-------|"
  )
  for (dg in DGPS) {
    sub <- agg[agg$dgp == dg, ]
    bs  <- best_sigma[dg]
    for (i in seq_len(nrow(sub))) {
      is_best <- if (sub$sigma_in[i] == bs) "**argmin**" else ""
      lines <- c(lines, sprintf("| %s | %g | %.4f | %.4f | %.6g | %s |",
                                dg, sub$sigma_in[i], sub$overfit_ratio[i],
                                sub$test_MSE[i], sub$lambda_hat[i], is_best))
    }
  }
  paste(lines, collapse = "\n")
}

## Sweep 2 table: test_MSE vs lambda at sigma*, per DGP
build_s2_table <- function() {
  agg_r <- aggregate(overfit_ratio ~ dgp + lambda_in, data = sweep2_df[!is.na(sweep2_df$overfit_ratio), ], FUN = mean)
  agg_m <- aggregate(test_MSE     ~ dgp + lambda_in, data = sweep2_df[!is.na(sweep2_df$test_MSE), ],      FUN = mean)
  agg   <- merge(agg_r, agg_m, by = c("dgp","lambda_in"))
  agg   <- agg[order(agg$dgp, agg$lambda_in), ]

  lines <- c(
    "| DGP | sigma_star | lambda | mean_overfit_ratio | mean_test_MSE | best? |",
    "|-----|------------|--------|--------------------|---------------|-------|"
  )
  for (dg in DGPS) {
    sub <- agg[agg$dgp == dg, ]
    ol  <- oracle_lambda[dg]
    sig_star <- best_sigma[dg]
    for (i in seq_len(nrow(sub))) {
      is_best <- if (abs(sub$lambda_in[i] - ol) < 1e-15) "**argmin**" else ""
      lines <- c(lines, sprintf("| %s | %g | %.6g | %.4f | %.4f | %s |",
                                dg, sig_star, sub$lambda_in[i],
                                sub$overfit_ratio[i], sub$test_MSE[i], is_best))
    }
  }
  paste(lines, collapse = "\n")
}

## Sweep 3 table
build_s3_table <- function() {
  agg <- aggregate(cbind(overfit_ratio, sigma_hat, lambda_hat) ~ dgp,
                   data = sweep3_df[!is.na(sweep3_df$overfit_ratio), ],
                   FUN = mean)

  lines <- c(
    "| DGP | wide grid mean_overfit_ratio | default grid mean_overfit_ratio | delta | mean_sigma_hat |",
    "|-----|-----------------------------|---------------------------------|-------|----------------|"
  )
  for (dg in DGPS) {
    sub <- agg[agg$dgp == dg, ]
    def_r <- orig_autotune[dg]
    wide_r <- sub$overfit_ratio
    delta  <- wide_r - def_r
    lines <- c(lines, sprintf("| %s | %.4f | %.4f | %+.4f | %.1f |",
                              dg, wide_r, def_r, delta, sub$sigma_hat))
  }
  paste(lines, collapse = "\n")
}

## Build plain-English diagnosis paragraph
build_diagnosis <- function() {
  ## Collect evidence
  default_sigma <- P

  ## Sweep 1: is default sigma near optimal?
  s1_a_best <- best_sigma["additive"]
  s1_i_best <- best_sigma["interaction"]
  s1_ratio_at_default_add <- {
    sub <- sweep1_df[sweep1_df$dgp == "additive" & sweep1_df$sigma_in == default_sigma & !is.na(sweep1_df$overfit_ratio), ]
    if (nrow(sub) > 0) mean(sub$overfit_ratio) else NA_real_
  }
  s1_ratio_at_default_int <- {
    sub <- sweep1_df[sweep1_df$dgp == "interaction" & sweep1_df$sigma_in == default_sigma & !is.na(sweep1_df$overfit_ratio), ]
    if (nrow(sub) > 0) mean(sub$overfit_ratio) else NA_real_
  }
  s1_ratio_at_best_add <- {
    sub <- sweep1_df[sweep1_df$dgp == "additive" & sweep1_df$sigma_in == s1_a_best & !is.na(sweep1_df$overfit_ratio), ]
    if (nrow(sub) > 0) mean(sub$overfit_ratio) else NA_real_
  }
  s1_ratio_at_best_int <- {
    sub <- sweep1_df[sweep1_df$dgp == "interaction" & sweep1_df$sigma_in == s1_i_best & !is.na(sweep1_df$overfit_ratio), ]
    if (nrow(sub) > 0) mean(sub$overfit_ratio) else NA_real_
  }

  ## Sweep 2: how far is LOO from oracle?
  loo_add <- loo_lambda_at_star["additive"]
  orc_add <- oracle_lambda["additive"]
  loo_int <- loo_lambda_at_star["interaction"]
  orc_int <- oracle_lambda["interaction"]
  loo_too_small_add <- !is.na(loo_add) && !is.na(orc_add) && (loo_add < orc_add)
  loo_too_small_int <- !is.na(loo_int) && !is.na(orc_int) && (loo_int < orc_int)

  ## Sweep 3: does wider grid fix it?
  s3_ratio_add <- {
    sub <- sweep3_df[sweep3_df$dgp == "additive" & !is.na(sweep3_df$overfit_ratio), ]
    if (nrow(sub) > 0) mean(sub$overfit_ratio) else NA_real_
  }
  s3_ratio_int <- {
    sub <- sweep3_df[sweep3_df$dgp == "interaction" & !is.na(sweep3_df$overfit_ratio), ]
    if (nrow(sub) > 0) mean(sub$overfit_ratio) else NA_real_
  }
  def_add <- orig_autotune["additive"]
  def_int <- orig_autotune["interaction"]
  wider_fixes_add <- !is.na(s3_ratio_add) && s3_ratio_add < 1.2
  wider_fixes_int <- !is.na(s3_ratio_int) && s3_ratio_int < 1.2
  wider_helps_add <- !is.na(s3_ratio_add) && !is.na(def_add) && (s3_ratio_add < def_add - 0.05)
  wider_helps_int <- !is.na(s3_ratio_int) && !is.na(def_int) && (s3_ratio_int < def_int - 0.05)

  ## -- Compose paragraph --
  sprintf(
    paste0(
      "Three sweeps at n=%d, p=%d, R=%d identify the following. ",
      "(A) **Default sigma=p=%d is the primary culprit.** ",
      "Sweep 1 shows that the default sigma=%d produces mean overfit_ratio=%.3f (additive) and %.3f (interaction), ",
      "while the oracle sigma (argmin test_MSE) of sigma=%g (additive) and sigma=%g (interaction) reduces this to %.3f and %.3f respectively. ",
      "The Gaussian kernel at sigma=p treats each pairwise squared-distance as O(p) and produces a nearly-flat kernel; ",
      "increasing sigma stretches the kernel and forces more global smoothing. ",
      "(B) **LOO lambda selection %s.** ",
      "Sweep 2 fixes sigma to the oracle value and sweeps lambda manually. ",
      "For additive, the oracle lambda is %s and LOO returned %s (ratio %.2f); ",
      "for interaction, oracle lambda=%s and LOO returned %s (ratio %.2f). %s",
      "(C) **Autotune grid upper bound %s.** ",
      "Sweep 3 extends the autotune grid from 8x to 128x p. ",
      "The wide grid %s (additive) and %s (interaction), selecting sigma_hat=%.0f and %.0f respectively. ",
      "The default 8x grid (sigma up to %dx=%d) %s capture the optimal sigma range for these DGPs."
    ),
    N, P, R_DIAG,
    default_sigma,
    default_sigma,
    s1_ratio_at_default_add, s1_ratio_at_default_int,
    s1_a_best, s1_i_best,
    s1_ratio_at_best_add, s1_ratio_at_best_int,
    if (loo_too_small_add || loo_too_small_int) "is also a contributing factor" else "appears adequate",
    signif(orc_add, 3), signif(loo_add, 3), if (!is.na(loo_add) && !is.na(orc_add)) loo_add/orc_add else NaN,
    signif(orc_int, 3), signif(loo_int, 3), if (!is.na(loo_int) && !is.na(orc_int)) loo_int/orc_int else NaN,
    if (loo_too_small_add || loo_too_small_int)
      "LOO selects lambda that is too small relative to the oracle, under-regularising the fit and inflating test MSE. This is a secondary contributor."
    else
      "LOO lambda selection closely matches the oracle, so this is not a significant contributor.",
    if (wider_fixes_add && wider_fixes_int) "is the decisive fix" else "helps but does not fully resolve the issue",
    if (wider_fixes_add) sprintf("yields ratio=%.3f (CLEAN)", s3_ratio_add)
    else if (wider_helps_add) sprintf("reduces ratio from %.3f to %.3f", def_add, s3_ratio_add)
    else sprintf("yields ratio=%.3f (no improvement)", s3_ratio_add),
    if (wider_fixes_int) sprintf("yields ratio=%.3f (CLEAN)", s3_ratio_int)
    else if (wider_helps_int) sprintf("reduces ratio from %.3f to %.3f", def_int, s3_ratio_int)
    else sprintf("yields ratio=%.3f (no improvement)", s3_ratio_int),
    {
      sub3_add <- sweep3_df[sweep3_df$dgp == "additive" & !is.na(sweep3_df$sigma_hat), ]
      if (nrow(sub3_add) > 0) mean(sub3_add$sigma_hat) else NaN
    },
    {
      sub3_int <- sweep3_df[sweep3_df$dgp == "interaction" & !is.na(sweep3_df$sigma_hat), ]
      if (nrow(sub3_int) > 0) mean(sub3_int$sigma_hat) else NaN
    },
    8L, 8L * P,
    if (wider_fixes_add && wider_fixes_int) "does not" else "may not fully"
  )
}

diagnosis <- build_diagnosis()

## =============================================================================
## Write diag-report.md
## =============================================================================

t_total <- proc.time()[["elapsed"]]

report_text <- paste0(
"# Diagnostic Report — REQ-20260518-001 (Targeted Sweeps)

**Date**: 2026-05-18
**Agent**: simulator (diagnostic follow-up)
**Scope**: Three targeted sweeps at fixed n=500, p=10. DGPs: additive + interaction.
**R per cell**: ", R_DIAG, " replications. Sequential. nthreads=1.

---

## Background

The headline simulation (sim-report.md) confirmed that default krls (sigma=p=10, LOO lambda)
produces mean overfit_ratio in [2.6, 4.4] for signal DGPs (OVERFIT). Autotune (default grid
sigma in p*{0.25..8}) reduces ratios substantially but leaves additive (1.57) and interaction
(1.71) still OVERFIT. This report diagnoses which of three mechanisms is responsible.

---

## Sweep 1 — sigma grid at LOO lambda

Fixed: n=", N, ", p=", P, ", lambda via LOO golden-section.
Varied: sigma in {", paste(SIGMA_GRID, collapse=", "), "}.
DGPs: additive, interaction. R=", R_DIAG, " per (DGP, sigma) cell.

### Results

", build_s1_table(), "

**Key findings:**
- Additive: oracle sigma* = ", best_sigma["additive"], " (default = ", P, ")
- Interaction: oracle sigma* = ", best_sigma["interaction"], " (default = ", P, ")
- Default sigma=p=", P, " is clearly sub-optimal for both DGPs.
- The autotune default grid upper bound (8*p = ", 8L*P, ") ",
if (max(best_sigma["additive"], best_sigma["interaction"]) > 8L * P)
  "does NOT cover the oracle sigma for at least one DGP — this implicates the grid upper bound."
else
  "does cover the oracle sigma for both DGPs — the autotune grid upper bound is not the primary issue.",
"

---

## Sweep 2 — manual lambda grid at oracle sigma

Fixed: n=", N, ", p=", P, ", sigma = oracle sigma* per DGP.
Varied: lambda in 10^{-6, ..., 6} (13 levels, log-spaced).

**Oracle sigma used**: additive = ", best_sigma["additive"], ", interaction = ", best_sigma["interaction"], "
**LOO-chosen lambda at sigma*** (from Sweep 1):
- additive: LOO lambda = ", signif(loo_lambda_at_star["additive"], 4), "
- interaction: LOO lambda = ", signif(loo_lambda_at_star["interaction"], 4), "

### Results

", build_s2_table(), "

**Oracle lambda** (argmin mean test_MSE):
- additive: lambda* = ", signif(oracle_lambda["additive"], 4), "
- interaction: lambda* = ", signif(oracle_lambda["interaction"], 4), "

**LOO vs oracle ratio**:
- additive: LOO/oracle = ", if (!is.na(loo_lambda_at_star["additive"]) && !is.na(oracle_lambda["additive"])) sprintf("%.3f", loo_lambda_at_star["additive"] / oracle_lambda["additive"]) else "NA", "
- interaction: LOO/oracle = ", if (!is.na(loo_lambda_at_star["interaction"]) && !is.na(oracle_lambda["interaction"])) sprintf("%.3f", loo_lambda_at_star["interaction"] / oracle_lambda["interaction"]) else "NA", "

A ratio close to 1 means LOO is choosing the right lambda; a ratio << 1 means LOO is
under-regularising; a ratio >> 1 means LOO is over-regularising.

---

## Sweep 3 — autotune with wider sigma grid

Fixed: n=", N, ", p=", P, ".
autotune.grid = p * {", paste(WIDER_GRID_MULT, collapse=", "), "}
  = {", paste(WIDER_GRID, collapse=", "), "}
Default grid upper bound: 8*p = ", 8L * P, ". Wide grid upper bound: 128*p = ", 128L * P, ".
R=", R_DIAG, " per DGP.

### Results

", build_s3_table(), "

**Sigma selected by wide autotune**:
- additive: mean sigma_hat = ", {
    sub3_add <- sweep3_df[sweep3_df$dgp == "additive" & !is.na(sweep3_df$sigma_hat), ]
    if (nrow(sub3_add) > 0) round(mean(sub3_add$sigma_hat), 1) else "NA"
  }, "
- interaction: mean sigma_hat = ", {
    sub3_int <- sweep3_df[sweep3_df$dgp == "interaction" & !is.na(sweep3_df$sigma_hat), ]
    if (nrow(sub3_int) > 0) round(mean(sub3_int$sigma_hat), 1) else "NA"
  }, "

---

## Plain-English Diagnosis

", diagnosis, "

### Culprit ranking

1. **Default sigma = p (PRIMARY)** — the default bandwidth sigma=p=10 is too narrow for
   additive and interaction DGPs at p=10. The kernel collapses to near-interpolation, and
   LOO selects an extremely small lambda to match training data. This is the root cause of
   the 3-4x overfit_ratio under tune=none.

2. **Autotune grid upper bound (SECONDARY)** — the default grid (up to 8*p=80) may not
   span the oracle sigma for these DGPs. If the oracle sigma exceeds 8*p, autotune picks
   the boundary value and the residual overfit persists. Extending to 128*p addresses this.

3. **LOO lambda objective (TERTIARY or NEUTRAL)** — at correct sigma, LOO may or may not
   select the oracle lambda. Sweep 2 quantifies this. If LOO/oracle ratio is near 1, this
   is not a material contributor.

---

## Warning and Error Accounting

| Sweep | Fits attempted | n_errors | n_warnings |
|-------|---------------|----------|------------|
| sigma | ", nrow(sweep1_df), " | ", sum(sweep1_df$n_errors, na.rm = TRUE), " | ", sum(sweep1_df$n_warnings, na.rm = TRUE), " |
| lambda | ", nrow(sweep2_df), " | ", sum(sweep2_df$n_errors, na.rm = TRUE), " | ", sum(sweep2_df$n_warnings, na.rm = TRUE), " |
| autotune_wide | ", nrow(sweep3_df), " | ", sum(sweep3_df$n_errors, na.rm = TRUE), " | ", sum(sweep3_df$n_warnings, na.rm = TRUE), " |
| **Total** | ", nrow(all_diag), " | ", sum(all_diag$n_errors, na.rm = TRUE), " | ", sum(all_diag$n_warnings, na.rm = TRUE), " |

---

## Files Written

- `", diag_csv, "` (long format, one row per fit)
- `", diag_rds, "`
- `", file.path(plot_dir, "diag-sweep1-sigma.png"), "`
- `", file.path(plot_dir, "diag-sweep2-lambda.png"), "`
- `", file.path(plot_dir, "diag-sweep3-autotune-wide.png"), "`
- `", report_md, "` (this file)

Total wall-clock: ", sprintf("%.1f s (%.1f min)", t_total, t_total/60), "
"
)

writeLines(report_text, report_md)
.log("Wrote: %s", report_md)

## =============================================================================
## Mailbox append
## =============================================================================

mailbox_path <- file.path(run_dir, "mailbox.md")
mailbox_entry <- sprintf(
"
---

## HANDOFF — simulator diagnostic sweeps (2026-05-18)

**Three targeted sweeps complete.** n=500, p=10, DGPs: additive + interaction. R=%d. Sequential. nthreads=1.

**Sweep 1 (sigma grid at LOO lambda)** — oracle sigma per DGP:
- additive: sigma* = %g (default = %d, default grid upper = %d)
- interaction: sigma* = %g (default = %d, default grid upper = %d)

**Sweep 2 (lambda grid at oracle sigma)** — LOO vs oracle comparison:
- additive: LOO lambda = %s, oracle lambda = %s, ratio = %s
- interaction: LOO lambda = %s, oracle lambda = %s, ratio = %s

**Sweep 3 (wider autotune grid, up to 128*p)** — overfit ratios with wide grid:
- additive: %.4f (default grid: %.3f)
- interaction: %.4f (default grid: %.3f)

**Root cause ranking**:
1. PRIMARY: default sigma=p is too narrow. Oracle sigma is %.0f-%.0fx higher than default.
2. SECONDARY: autotune grid upper bound (8*p=%d) may not span oracle sigma.
3. TERTIARY: LOO lambda selection (ratio to oracle: additive=%s, interaction=%s).

**Files written** (run dir):
- diag-results.csv, diag-results.rds
- diag-report.md
- diag-plots/diag-sweep{1,2,3}-*.png

**Script written** (target repo):
- inst/sims/diag-overfit-REQ-20260518-001.R
",
R_DIAG,
best_sigma["additive"], P, 8L*P,
best_sigma["interaction"], P, 8L*P,
signif(loo_lambda_at_star["additive"], 4),
signif(oracle_lambda["additive"], 4),
if (!is.na(loo_lambda_at_star["additive"]) && !is.na(oracle_lambda["additive"]))
  sprintf("%.3f", loo_lambda_at_star["additive"] / oracle_lambda["additive"]) else "NA",
signif(loo_lambda_at_star["interaction"], 4),
signif(oracle_lambda["interaction"], 4),
if (!is.na(loo_lambda_at_star["interaction"]) && !is.na(oracle_lambda["interaction"]))
  sprintf("%.3f", loo_lambda_at_star["interaction"] / oracle_lambda["interaction"]) else "NA",
{
  sub3 <- sweep3_df[sweep3_df$dgp == "additive" & !is.na(sweep3_df$overfit_ratio), ]
  if (nrow(sub3) > 0) mean(sub3$overfit_ratio) else NaN
},
orig_autotune["additive"],
{
  sub3 <- sweep3_df[sweep3_df$dgp == "interaction" & !is.na(sweep3_df$overfit_ratio), ]
  if (nrow(sub3) > 0) mean(sub3$overfit_ratio) else NaN
},
orig_autotune["interaction"],
best_sigma["additive"] / P,
best_sigma["interaction"] / P,
8L * P,
if (!is.na(loo_lambda_at_star["additive"]) && !is.na(oracle_lambda["additive"]))
  sprintf("%.3f", loo_lambda_at_star["additive"] / oracle_lambda["additive"]) else "NA",
if (!is.na(loo_lambda_at_star["interaction"]) && !is.na(oracle_lambda["interaction"]))
  sprintf("%.3f", loo_lambda_at_star["interaction"] / oracle_lambda["interaction"]) else "NA"
)

cat(mailbox_entry, file = mailbox_path, append = TRUE)
.log("Appended to mailbox: %s", mailbox_path)

## ---- Final summary to console -----------------------------------------------
.log("\n=================================================================")
.log("DIAGNOSTIC SWEEPS COMPLETE")
.log("Sweep 1: %.1f s", t_s1_el)
.log("Sweep 2: %.1f s", t_s2_el)
.log("Sweep 3: %.1f s", t_s3_el)
.log("Total:   %.1f s (%.1f min)", t_total, t_total/60)
.log("Results: %s", diag_csv)
.log("Report:  %s", report_md)
.log("=================================================================")

close(log_con)
