## v0.26 KRLS speed baseline -- fixed benchmark grid for KRLS speedup work.
## Patterned on inst/sims/v0.26-speed-baseline.R (ares).
##
## Grid:
##   * Gaussian + ls (byte-locked: must produce identical digest pre/post)
##   * Several (n, p) cells covering typical / large fits
##   * Setups: default, autotune (fast), cv5 (lambda.method=cv)
##
## Outputs:
##   inst/sims/results/v0.26-krls-speed-baseline.csv
##   inst/sims/results/v0.26-krls-speed-baseline.rds
##   inst/sims/results/v0.26-krls-iter-NN.csv         (iter > 0)
##   inst/sims/results/v0.26-krls-iter-NN.rds         (iter > 0)
##
## Usage:
##   Rscript inst/sims/v0.26-krls-speed-baseline.R          # baseline
##   Rscript inst/sims/v0.26-krls-speed-baseline.R 1        # iter=1
##   Rscript inst/sims/v0.26-krls-speed-baseline.R 7 quick  # iter=7, 2 reps

suppressPackageStartupMessages({
  library(bench)
})
devtools::load_all("/home/jack/Dropbox/roadrunner", quiet = TRUE)

args   <- commandArgs(TRUE)
iter   <- if (length(args) >= 1L) as.integer(args[1]) else 0L
quick  <- length(args) >= 2L && tolower(args[2]) == "quick"
n_reps <- if (quick) 2L else 3L

out_dir <- "/home/jack/Dropbox/roadrunner/inst/sims/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
tag <- if (iter == 0L) {
  "v0.26-krls-speed-baseline"
} else {
  sprintf("v0.26-krls-iter-%02d", iter)
}
csv_path <- file.path(out_dir, paste0(tag, ".csv"))
rds_path <- file.path(out_dir, paste0(tag, ".rds"))

# ---- DGP ----------------------------------------------------------------
.namex <- function(x) { colnames(x) <- paste0("V", seq_len(ncol(x))); x }

dgp_friedman1_g <- function(n, p, seed) {
  set.seed(seed)
  x <- matrix(stats::runif(n * p), n, p)
  f <- 10 * sin(pi * x[, 1] * x[, 2]) + 20 * (x[, 3] - 0.5)^2 +
       10 * x[, 4] + 5 * x[, 5]
  sigma <- sqrt(stats::var(f) / 3)
  y <- f + stats::rnorm(n, sd = sigma)
  list(x = .namex(x), y = y)
}

# Train + held-out test split for OOS RMSE (a "test set performance"
# guard alongside the byte-id digest).
dgp_train_test <- function(n, p, seed) {
  tr <- dgp_friedman1_g(n,        p, seed)
  te <- dgp_friedman1_g(max(200L, n %/% 2L), p, seed + 7919L)
  list(x = tr$x, y = tr$y, x_te = te$x, y_te = te$y)
}

# ---- Grid (cells = setup x (n, p)) --------------------------------------
cells_def <- list(
  list(n = 400L,  p = 5L,  seed = 200001L),
  list(n = 800L,  p = 10L, seed = 200002L),
  list(n = 1500L, p = 10L, seed = 200003L)
)

setups <- list(
  default = list(
    label = "default",
    args_fn = function(cell) list()  # lambda.method default = "loo"
  ),
  gcv = list(
    label = "gcv",
    args_fn = function(cell) list(lambda.method = "gcv")
  ),
  autotune_fast = list(
    label = "autotune_fast",
    args_fn = function(cell) list(autotune = TRUE,
                                  autotune.speed = "fast",
                                  autotune.nthreads = 12L)
  )
)

# ---- Helpers -------------------------------------------------------------
fit_one <- function(cell, setup, nthreads) {
  args <- c(list(X = cell$x, y = cell$y, nthreads = as.integer(nthreads)),
            setup$args_fn(cell))
  do.call(krls, args)
}

# Canonical numeric digest helper.
round_d <- function(z) {
  if (is.null(z)) return("NULL")
  z <- as.numeric(z)
  z[!is.finite(z)] <- NA_real_
  z <- signif(z, 12)
  paste(z, collapse = ",")
}

fit_digest <- function(fit) {
  sigma_str <- if (!is.null(fit$sigma_vec)) {
    round_d(fit$sigma_vec)
  } else {
    round_d(fit$sigma)
  }
  digest_str <- paste(
    round_d(fit$coeffs),
    round_d(fit$lambda),
    sigma_str,
    round_d(fit$Looe),
    round_d(fit$fitted),
    sep = "|"
  )
  list(
    digest = digest_str,
    coef = as.numeric(fit$coeffs),
    lambda = as.numeric(fit$lambda),
    sigma = as.numeric(fit$sigma),
    Looe = as.numeric(fit$Looe),
    R2 = as.numeric(fit$R2)
  )
}

oos_rmse <- function(fit, x_te, y_te) {
  yhat <- predict(fit, newdata = x_te)
  if (is.list(yhat)) yhat <- yhat$fit
  sqrt(mean((as.numeric(yhat) - as.numeric(y_te))^2))
}

run_cell <- function(cell_data, cell_meta, setup) {
  cell_id <- sprintf("n%d_p%d__%s",
                     cell_meta$n, cell_meta$p, setup$label)
  cat(sprintf("[%s] timing (nthreads=4, reps=%d)... ", cell_id, n_reps))
  cell_for_fit <- list(x = cell_data$x, y = cell_data$y)
  t0 <- Sys.time()
  bm <- bench::mark(
    fit = fit_one(cell_for_fit, setup, nthreads = 12L),
    iterations = n_reps,
    min_iterations = n_reps,
    max_iterations = n_reps,
    check = FALSE,
    filter_gc = FALSE
  )
  wall <- as.numeric(Sys.time() - t0, units = "secs")
  cat(sprintf("median=%.3fs (wall=%.1fs)\n",
              as.numeric(bm$median), wall))

  # Determinism reference fit at nthreads = 1
  cat(sprintf("[%s] nthreads=1 digest + OOS... ", cell_id))
  fit1 <- fit_one(cell_for_fit, setup, nthreads = 1L)
  dig1 <- fit_digest(fit1)
  rmse_te <- oos_rmse(fit1, cell_data$x_te, cell_data$y_te)
  cat(sprintf("rmse_te=%.4f\n", rmse_te))

  # Determinism check: nthreads = 4
  fit4 <- fit_one(cell_for_fit, setup, nthreads = 12L)
  dig4 <- fit_digest(fit4)
  det_ok <- identical(dig1$digest, dig4$digest)
  if (!det_ok) {
    warning(sprintf("DETERMINISM BREAK for %s: nthreads=1 != nthreads=4",
                    cell_id))
  }

  list(
    row = data.frame(
      cell_id = cell_id,
      n = cell_meta$n,
      p = cell_meta$p,
      setup = setup$label,
      n_reps = n_reps,
      median_s = as.numeric(bm$median),
      min_s = as.numeric(bm$min),
      mem_alloc_b = as.numeric(bm$mem_alloc),
      det_ok = det_ok,
      rmse_te = rmse_te,
      lambda = dig1$lambda,
      sigma = dig1$sigma[1],
      R2 = dig1$R2,
      stringsAsFactors = FALSE
    ),
    digest = dig1
  )
}

# ---- Run -----------------------------------------------------------------
all_cells <- list()
for (ce in cells_def) {
  for (st in setups) {
    all_cells[[length(all_cells) + 1L]] <- list(cell_meta = ce, setup = st)
  }
}

cat(sprintf("=== %s : %d cells, %d reps each ===\n",
            tag, length(all_cells), n_reps))

rows <- vector("list", length(all_cells))
digs <- vector("list", length(all_cells))
names(digs) <- vapply(all_cells, function(s)
  sprintf("n%d_p%d__%s", s$cell_meta$n, s$cell_meta$p, s$setup$label),
  character(1))

for (i in seq_along(all_cells)) {
  s <- all_cells[[i]]
  cd <- dgp_train_test(s$cell_meta$n, s$cell_meta$p, s$cell_meta$seed)
  res <- run_cell(cd, s$cell_meta, s$setup)
  rows[[i]] <- res$row
  digs[[i]] <- res$digest
}

df <- do.call(rbind, rows)
df$tag <- tag
df$iter <- iter
df$timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

write.csv(df, csv_path, row.names = FALSE)
saveRDS(list(meta = list(tag = tag, iter = iter,
                         timestamp = df$timestamp[1],
                         n_reps = n_reps),
             digests = digs,
             timings = df),
        rds_path)

cat(sprintf("\nWrote: %s\nWrote: %s\n", csv_path, rds_path))

# ---- Comparison (if not baseline) ----------------------------------------
if (iter > 0L) {
  base_rds <- file.path(out_dir, "v0.26-krls-speed-baseline.rds")
  if (file.exists(base_rds)) {
    base <- readRDS(base_rds)
    cmp <- merge(base$timings[, c("cell_id", "median_s", "rmse_te")],
                 df[, c("cell_id", "median_s", "rmse_te")],
                 by = "cell_id", suffixes = c("_base", "_iter"))
    cmp$speedup <- cmp$median_s_base / cmp$median_s_iter
    cmp$pct <- 100 * (cmp$speedup - 1)
    cmp$rmse_delta <- cmp$rmse_te_iter - cmp$rmse_te_base
    cmp$digest_ok <- vapply(cmp$cell_id, function(k) {
      isTRUE(identical(base$digests[[k]]$digest, digs[[k]]$digest))
    }, logical(1))
    cmp$coef_linf <- vapply(cmp$cell_id, function(k) {
      a <- base$digests[[k]]$coef
      b <- digs[[k]]$coef
      if (is.null(a) || is.null(b) || length(a) != length(b)) return(Inf)
      max(abs(a - b))
    }, numeric(1))
    cmp$lambda_delta <- vapply(cmp$cell_id, function(k) {
      abs(base$digests[[k]]$lambda - digs[[k]]$lambda)
    }, numeric(1))
    cat("\n=== Comparison vs baseline ===\n")
    print(cmp[order(-cmp$speedup), ], row.names = FALSE)
    gmean <- exp(mean(log(cmp$speedup)))
    cat(sprintf("\nGeometric-mean speedup: %.4fx (%.2f%%)\n",
                gmean, 100 * (gmean - 1)))
    write.csv(cmp, file.path(out_dir, paste0(tag, "-cmp.csv")),
              row.names = FALSE)
  } else {
    cat("(no krls baseline.rds found; skipping comparison)\n")
  }
}

cat("\nDONE.\n")
