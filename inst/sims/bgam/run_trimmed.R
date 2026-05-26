# Trimmed bgam vs baselines sim — leader-supervised
# 4 cells, 50 reps each, gaussian DGP1 only
# Target wall: < 10 minutes serial

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
})

set.seed(20260526)

# DGP1: sparse smooth additive
make_data <- function(n, p, sigma) {
  X <- matrix(runif(n * p, -1, 1), n, p)
  y <- sin(2 * pi * X[, 1]) + 0.5 * X[, 2]^2 + sigma * rnorm(n)
  list(X = X, y = y)
}

# Cells
cells <- list(
  list(n = 200,  p = 10, sigma = 0.5),
  list(n = 200,  p = 50, sigma = 0.5),
  list(n = 1000, p = 10, sigma = 0.5),
  list(n = 1000, p = 50, sigma = 0.5)
)

reps <- 50
results <- list()

cat(sprintf("[sim] %d cells x %d reps x 4 learners\n", length(cells), reps))
t0_total <- Sys.time()

for (ci in seq_along(cells)) {
  cell <- cells[[ci]]
  cat(sprintf("[cell %d/%d] n=%d p=%d sigma=%.2f\n",
              ci, length(cells), cell$n, cell$p, cell$sigma))
  t0_cell <- Sys.time()

  cell_results <- vector("list", reps)
  for (r in seq_len(reps)) {
    set.seed(20260526 + 1000 * ci + r)
    tr <- make_data(cell$n, cell$p, cell$sigma)
    te <- make_data(cell$n, cell$p, cell$sigma)

    out <- list()

    # bgam
    t1 <- Sys.time()
    fit_bgam <- tryCatch(
      bgam(tr$X, tr$y, mstop = 100, nu = 0.1, autotune = FALSE),
      error = function(e) NULL
    )
    if (!is.null(fit_bgam)) {
      yhat <- predict(fit_bgam, newdata = te$X)
      out$bgam <- list(rmse = sqrt(mean((te$y - yhat)^2)),
                       time = as.numeric(Sys.time() - t1, units = "secs"))
    }

    # ares
    t1 <- Sys.time()
    fit_ares <- tryCatch(ares(tr$X, tr$y), error = function(e) NULL)
    if (!is.null(fit_ares)) {
      yhat <- predict(fit_ares, newdata = te$X)
      out$ares <- list(rmse = sqrt(mean((te$y - yhat)^2)),
                       time = as.numeric(Sys.time() - t1, units = "secs"))
    }

    # ols (with degree-3 poly only on first 5 features to keep deg controlled)
    t1 <- Sys.time()
    df_ols <- data.frame(y = tr$y, tr$X)
    names(df_ols) <- c("y", sprintf("x%d", seq_len(cell$p)))
    df_te_ols <- data.frame(te$X)
    names(df_te_ols) <- sprintf("x%d", seq_len(cell$p))
    fit_ols <- tryCatch(ols(y ~ ., data = df_ols), error = function(e) NULL)
    if (!is.null(fit_ols)) {
      yhat <- predict(fit_ols, newdata = df_te_ols)
      out$ols <- list(rmse = sqrt(mean((te$y - yhat)^2)),
                      time = as.numeric(Sys.time() - t1, units = "secs"))
    }

    cell_results[[r]] <- out
  }

  results[[ci]] <- list(cell = cell, reps = cell_results,
                        wall = as.numeric(Sys.time() - t0_cell, units = "secs"))
  cat(sprintf("  done in %.1fs\n", results[[ci]]$wall))
}

cat(sprintf("[sim] total wall: %.1fs\n", as.numeric(Sys.time() - t0_total, units = "secs")))

# Aggregate
summary_rows <- list()
for (ci in seq_along(results)) {
  cell <- results[[ci]]$cell
  rr <- results[[ci]]$reps
  for (lr in c("bgam", "ares", "ols")) {
    rmses <- sapply(rr, function(x) if (is.null(x[[lr]])) NA else x[[lr]]$rmse)
    times <- sapply(rr, function(x) if (is.null(x[[lr]])) NA else x[[lr]]$time)
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      n = cell$n, p = cell$p, sigma = cell$sigma, learner = lr,
      rmse_mean = mean(rmses, na.rm = TRUE),
      rmse_mcse = sd(rmses, na.rm = TRUE) / sqrt(sum(!is.na(rmses))),
      time_mean = mean(times, na.rm = TRUE),
      n_ok = sum(!is.na(rmses))
    )
  }
}
sum_df <- do.call(rbind, summary_rows)
print(sum_df, row.names = FALSE)

saveRDS(list(cells = cells, results = results, summary = sum_df,
             pkg_version = utils::packageVersion("roadrunner")),
        "inst/sims/results/bgam-trimmed-0.0.0.9059.rds")

cat("[sim] wrote inst/sims/results/bgam-trimmed-0.0.0.9059.rds\n")
