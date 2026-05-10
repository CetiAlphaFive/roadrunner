# Monte Carlo benchmark for ares vs earth.
# Smoke-grid by default. Run with ARES_SIM_FULL=1 for the full grid.
#
# Reproducible: master seed below; per-cell seed = master + cell_index.

suppressPackageStartupMessages({
  library(ares)
  if (!requireNamespace("earth", quietly = TRUE))
    stop("'earth' is required for the simulation comparison.")
  library(earth)
})

source(file.path("inst/sims", "dgp_friedman.R"))
source(file.path("inst/sims", "dgp_additive.R"))
source(file.path("inst/sims", "dgp_interaction.R"))

MASTER_SEED <- 20260509L
FULL <- nzchar(Sys.getenv("ARES_SIM_FULL"))

dgps <- list(
  friedman1   = dgp_friedman1,
  friedman2   = dgp_friedman2,
  friedman3   = dgp_friedman3,
  additive    = dgp_additive,
  interaction = dgp_interaction
)

if (FULL) {
  grid <- expand.grid(
    dgp     = names(dgps),
    n       = c(50L, 200L, 1000L, 5000L),
    p       = c(5L, 10L, 25L),
    snr     = c(1, 3, 10),
    degree  = c(1L, 2L),
    nrep    = 20L,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
} else {
  grid <- expand.grid(
    dgp     = c("friedman1", "additive", "interaction"),
    n       = c(200L, 1000L),
    p       = c(5L, 10L),
    snr     = 3,
    degree  = c(1L, 2L),
    nrep    = 5L,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
}

run_one <- function(row_idx, dgp, n, p, snr, degree, nrep) {
  # Skip cells where p is incompatible with DGP
  if (dgp == "friedman2" && p < 5L) return(NULL)
  if (dgp == "friedman3" && p < 5L) return(NULL)

  reps <- replicate(nrep, simplify = FALSE, expr = {
    NULL
  })
  cell_seed <- MASTER_SEED + row_idx
  out <- vector("list", nrep)
  for (rep in seq_len(nrep)) {
    rep_seed <- as.integer(cell_seed * 17 + rep) %% .Machine$integer.max
    d <- dgps[[dgp]](n = n, p = p, snr = snr, seed = rep_seed)
    # ares (single-thread) and ares (multi-thread, capped at 2 for CRAN safety)
    t_a1 <- system.time(af1 <- ares::ares(d$x, d$y, degree = degree, nthreads = 1))[["elapsed"]]
    t_a2 <- system.time(af2 <- ares::ares(d$x, d$y, degree = degree, nthreads = 2))[["elapsed"]]
    t_e  <- system.time(ef  <- earth::earth(d$x, d$y, degree = degree))[["elapsed"]]

    nt_a <- length(af1$selected.terms)
    nt_e <- length(ef$selected.terms)
    pred_a <- predict(af1, d$x); pred_e <- predict(ef, d$x)
    rmse_diff <- sqrt(mean((pred_a - pred_e)^2)) / max(stats::sd(d$y), 1e-9)

    out[[rep]] <- data.frame(
      cell_index = row_idx, rep = rep, dgp = dgp, n = n, p = p, snr = snr,
      degree = degree,
      t_ares_st = t_a1, t_ares_mt = t_a2, t_earth = t_e,
      rss_ares = af1$rss, rss_earth = ef$rss,
      gcv_ares = af1$gcv, gcv_earth = ef$gcv,
      nterms_ares = nt_a, nterms_earth = nt_e,
      pred_rmse_diff = rmse_diff,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

cat("Running", nrow(grid), "cells with smoke=", !FULL, "\n")
results <- vector("list", nrow(grid))
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  cat(sprintf("[%d/%d] %s n=%d p=%d snr=%g deg=%d nrep=%d\n",
              i, nrow(grid), g$dgp, g$n, g$p, g$snr, g$degree, g$nrep))
  r <- run_one(i, g$dgp, g$n, g$p, g$snr, g$degree, g$nrep)
  if (!is.null(r)) results[[i]] <- r
}
sim_raw <- do.call(rbind, results)

# Per-cell summary
cells <- unique(sim_raw[, c("cell_index", "dgp", "n", "p", "snr", "degree")])
sim_summary <- do.call(rbind, lapply(seq_len(nrow(cells)), function(k) {
  c0 <- cells[k, ]
  rows <- sim_raw[sim_raw$cell_index == c0$cell_index, ]
  data.frame(
    c0,
    nrep = nrow(rows),
    rss_match    = stats::median(abs(rows$rss_ares - rows$rss_earth) / pmax(rows$rss_earth, 1e-12)),
    gcv_ratio    = stats::median(rows$gcv_ares / rows$gcv_earth),
    nterms_diff  = stats::median(rows$nterms_ares - rows$nterms_earth),
    pred_rmse    = stats::median(rows$pred_rmse_diff),
    speedup_st   = stats::median(rows$t_earth / rows$t_ares_st),
    speedup_mt   = stats::median(rows$t_earth / rows$t_ares_mt),
    stringsAsFactors = FALSE
  )
}))

dir.create("inst/sims/results", showWarnings = FALSE, recursive = TRUE)
utils::write.csv(sim_raw, "inst/sims/results/sim_raw.csv", row.names = FALSE)
utils::write.csv(sim_summary, "inst/sims/results/sim_summary.csv", row.names = FALSE)
cat("Wrote inst/sims/results/sim_summary.csv (", nrow(sim_summary), " cells)\n", sep = "")
print(sim_summary)
