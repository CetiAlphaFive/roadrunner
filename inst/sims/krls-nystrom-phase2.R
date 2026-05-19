#!/usr/bin/env Rscript
##
## EMP-PHASE2 - Nystrom speedup verification for v0.0.0.9043.
##
## Mandatory sim rules (this session):
##   1. Smoke test first (1 cell, R=2)
##   2. Intermediate CSV per cell (append after every fit)
##   3. Wall-clock cap 5 min total (abort if projected over)

suppressPackageStartupMessages({
  devtools::load_all(".", quiet = TRUE)
})

OUT_CSV     <- "inst/sims/results/krls-nystrom-phase2.csv"
WALL_CAP_S  <- 300
MASTER_SEED <- 20260519L
R_REPS      <- 5L

dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
if (file.exists(OUT_CSV)) file.remove(OUT_CSV)

dgp_additive <- function(n, p, seed) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  y <- as.numeric(sin(X[, 1]) + X[, 2]^2 + X[, 3] + 0.5 * rnorm(n))
  list(X = X, y = y)
}

run_one <- function(n, p, approx, rep_seed) {
  d_tr <- dgp_additive(n, p, rep_seed)
  d_te <- dgp_additive(1000L, p, rep_seed + 1L)
  t0 <- proc.time()[["elapsed"]]
  if (approx == "exact") {
    fit <- roadrunner::krls(X = d_tr$X, y = d_tr$y,
                            derivative = FALSE, vcov = FALSE)
  } else {
    fit <- roadrunner::krls(X = d_tr$X, y = d_tr$y,
                            approx = "nystrom", landmark_seed = 42L,
                            derivative = FALSE, vcov = FALSE)
  }
  wall <- proc.time()[["elapsed"]] - t0
  yhat <- as.numeric(predict(fit, newdata = d_te$X)$fit)
  test_RMSE <- sqrt(mean((d_te$y - yhat)^2))
  list(wall = wall, test_RMSE = test_RMSE,
       sigma = fit$sigma, lambda = fit$lambda,
       nystrom_m = if (approx == "nystrom") fit$nystrom_m else NA_integer_)
}

append_row <- function(row) {
  write.table(row, OUT_CSV, sep = ",", row.names = FALSE,
              col.names = !file.exists(OUT_CSV),
              append = file.exists(OUT_CSV), quote = FALSE)
}

## ---------- SMOKE ----------
cat("[smoke] n=500 p=5 approx=nystrom R=2 ...\n")
sm1 <- run_one(500L, 5L, "nystrom", MASTER_SEED + 9999L)
sm2 <- run_one(500L, 5L, "nystrom", MASTER_SEED + 10000L)
stopifnot(is.finite(sm1$wall), is.finite(sm2$wall),
          is.finite(sm1$test_RMSE), is.finite(sm2$test_RMSE))
cat(sprintf("[smoke OK] wall %.2fs/%.2fs  RMSE %.4f/%.4f\n",
            sm1$wall, sm2$wall, sm1$test_RMSE, sm2$test_RMSE))

## ---------- FULL ----------
cells <- expand.grid(
  n = c(2000L, 5000L),
  p = c(10L),
  approx = c("exact", "nystrom"),
  stringsAsFactors = FALSE
)

wall_start <- proc.time()[["elapsed"]]
for (ci in seq_len(nrow(cells))) {
  n <- cells$n[ci]; p <- cells$p[ci]; approx <- cells$approx[ci]
  for (r in seq_len(R_REPS)) {
    elapsed <- proc.time()[["elapsed"]] - wall_start
    if (elapsed > WALL_CAP_S - 30) {
      cat(sprintf("[abort] wall-cap approaching (%.0fs).\n", elapsed))
      stop("EMP-PHASE2: 5-min wall cap hit.")
    }
    # Pair exact vs nystrom on identical DGP draws at the same (n, p, rep).
    # Seed key excludes approx so both modes see the same (X, y) realization.
    n_idx <- match(n, c(2000L, 5000L))  # 1 or 2
    rep_seed <- MASTER_SEED + (n_idx - 1L) * R_REPS * 2L + (r - 1L) * 2L
    res <- run_one(n, p, approx, rep_seed)
    row <- data.frame(n = n, p = p, approx = approx, rep = r,
                      seed = rep_seed, wall_s = res$wall,
                      test_RMSE = res$test_RMSE,
                      sigma = res$sigma, lambda = res$lambda,
                      nystrom_m = res$nystrom_m)
    append_row(row)
    cat(sprintf("[%d/%d %s n=%d rep=%d] wall=%.2fs RMSE=%.4f\n",
                ci, nrow(cells), approx, n, r, res$wall, res$test_RMSE))
  }
}

## ---------- SUMMARY ----------
df <- read.csv(OUT_CSV)
agg <- aggregate(cbind(wall_s, test_RMSE) ~ n + p + approx, data = df,
                 FUN = median)
print(agg)

for (n_val in c(2000L, 5000L)) {
  w_e <- agg$wall_s[agg$n == n_val & agg$approx == "exact"]
  w_n <- agg$wall_s[agg$n == n_val & agg$approx == "nystrom"]
  r_e <- agg$test_RMSE[agg$n == n_val & agg$approx == "exact"]
  r_n <- agg$test_RMSE[agg$n == n_val & agg$approx == "nystrom"]
  cat(sprintf("\nn=%d: speedup=%.2fx  RMSE exact=%.4f  RMSE nystrom=%.4f  RMSE delta=%.1f%%\n",
              n_val, w_e/w_n, r_e, r_n, 100 * (r_n - r_e) / r_e))
}

## ---------- PAIRED DELTA ----------
chk <- merge(
  df[df$approx == "exact",   c("n", "rep", "test_RMSE")],
  df[df$approx == "nystrom", c("n", "rep", "test_RMSE")],
  by = c("n", "rep"), suffixes = c("_exact", "_nystrom")
)
chk$rmse_delta_pct <- 100 * (chk$test_RMSE_nystrom - chk$test_RMSE_exact) /
                      chk$test_RMSE_exact
agg2 <- aggregate(rmse_delta_pct ~ n, chk,
                  function(x) c(mean = mean(x), median = median(x),
                                min = min(x), max = max(x)))
print(agg2)

cat("\nResults at:", OUT_CSV, "\n")
