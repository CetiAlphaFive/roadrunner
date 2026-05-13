## v0.26 speed baseline -- fixed benchmark grid for the speed-improvement
## auto-research loop. Every iteration compares against the artifacts
## produced by THIS script.
##
## Grid: (family x dgp x setup) -- see `cells` below.
## Reps: 5 per cell via bench::mark.
## Threads: nthreads = 4 for timing; nthreads = 1 also recorded once per
##   cell for the determinism / digest reference (used by the iteration
##   gate to verify byte-identical fits).
##
## Outputs:
##   inst/sims/results/v0.26-speed-baseline.csv   (long; one row per cell)
##   inst/sims/results/v0.26-speed-baseline.rds   (list keyed by cell_id
##     with the nthreads=1 fit digest + the raw nthreads=4 timing)
##
## Re-run with `iter` >= 1 to produce iteration artifacts; the script
## writes:
##   inst/sims/results/v0.26-iter-NN.csv
##   inst/sims/results/v0.26-iter-NN.rds
##
## Usage:
##   Rscript inst/sims/v0.26-speed-baseline.R            # baseline (iter=0)
##   Rscript inst/sims/v0.26-speed-baseline.R 1          # iter=1
##   Rscript inst/sims/v0.26-speed-baseline.R 7 quick    # iter=7, 2 reps

suppressPackageStartupMessages({
  library(bench)
  library(mlbench)
})
devtools::load_all("/home/jack/Dropbox/roadrunner", quiet = TRUE)

args   <- commandArgs(TRUE)
iter   <- if (length(args) >= 1L) as.integer(args[1]) else 0L
quick  <- length(args) >= 2L && tolower(args[2]) == "quick"
n_reps <- if (quick) 2L else 5L

out_dir <- "/home/jack/Dropbox/roadrunner/inst/sims/results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
tag <- if (iter == 0L) {
  "v0.26-speed-baseline"
} else {
  sprintf("v0.26-iter-%02d", iter)
}
csv_path <- file.path(out_dir, paste0(tag, ".csv"))
rds_path <- file.path(out_dir, paste0(tag, ".rds"))

# ---- DGPs ----------------------------------------------------------------
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

dgp_twonorm_b <- function(n, p, seed) {
  set.seed(seed)
  # mlbench.twonorm is 20-D; if p != 20 we adapt by adding nuisance cols
  d <- mlbench::mlbench.twonorm(n)
  x <- d$x
  if (p > ncol(x)) {
    extra <- matrix(stats::rnorm(n * (p - ncol(x))), n, p - ncol(x))
    x <- cbind(x, extra)
  } else if (p < ncol(x)) {
    x <- x[, seq_len(p), drop = FALSE]
  }
  list(x = .namex(x), y = as.integer(d$classes) - 1L)
}

dgp_friedman1_p <- function(n, p, seed) {
  set.seed(seed)
  x <- matrix(stats::runif(n * p), n, p)
  eta <- 0.5 * sin(pi * x[, 1] * x[, 2]) + 0.3 * (x[, 3] - 0.5) +
         0.4 * x[, 4] - 0.2 * x[, 5]
  mu <- exp(eta + 1.5)  # mean ~ exp(1.5) ~ 4.5
  y <- stats::rpois(n, lambda = mu)
  list(x = .namex(x), y = y)
}

dgp_friedman1_gamma <- function(n, p, seed) {
  set.seed(seed)
  x <- matrix(stats::runif(n * p), n, p)
  eta <- 0.4 * sin(pi * x[, 1] * x[, 2]) + 0.3 * (x[, 3] - 0.5) +
         0.3 * x[, 4] - 0.2 * x[, 5]
  mu <- exp(eta + 1)
  shape <- 4
  y <- stats::rgamma(n, shape = shape, rate = shape / mu)
  list(x = .namex(x), y = y)
}

# ---- Grid (cells = family x dgp x setup) ---------------------------------
cells_def <- list(
  list(family = "gaussian", dgp_fn = dgp_friedman1_g,
       dgp = "friedman1", n = 500L,  p = 10L, seed = 100001L),
  list(family = "gaussian", dgp_fn = dgp_friedman1_g,
       dgp = "friedman1", n = 2000L, p = 10L, seed = 100002L),
  list(family = "gaussian", dgp_fn = dgp_friedman1_g,
       dgp = "friedman1_highdim", n = 1000L, p = 20L, seed = 100003L),
  list(family = "binomial", dgp_fn = dgp_twonorm_b,
       dgp = "twonorm", n = 1500L, p = 20L, seed = 100004L),
  list(family = "poisson", dgp_fn = dgp_friedman1_p,
       dgp = "friedman1_poisson", n = 1000L, p = 10L, seed = 100005L),
  list(family = "gamma", dgp_fn = dgp_friedman1_gamma,
       dgp = "friedman1_gamma", n = 1000L, p = 10L, seed = 100006L)
)

# Setups applied to every cell. `args_fn(cell)` returns a list of named
# args to pass to ares() in addition to (x = , y = , family = , nthreads = ).
setups <- list(
  default = list(
    label = "default",
    args_fn = function(cell) list()
  ),
  cv5 = list(
    label = "cv5",
    args_fn = function(cell) list(pmethod = "cv", nfold = 5L, ncross = 1L)
  ),
  autotune = list(
    label = "autotune",
    args_fn = function(cell) list(autotune = TRUE,
                                  autotune.speed = "fast",
                                  autotune.warmstart = TRUE)
  ),
  bag20 = list(
    label = "bag20",
    args_fn = function(cell) list(n.boot = 20L)
  )
)

# ---- Helpers -------------------------------------------------------------
fit_one <- function(cell, setup, nthreads) {
  d <- cell$dgp_fn(cell$n, cell$p, cell$seed)
  args <- c(list(x = d$x, y = d$y, family = cell$family,
                 nthreads = as.integer(nthreads)),
            setup$args_fn(cell))
  do.call(ares, args)
}

# Digest of a fit -- bag-aware. For non-bagged, capture forward+backward
# selected basis, coef, gcv. For bagged, capture aggregate digests across
# replicates (each replicate is itself a small ares fit).
fit_digest <- function(fit) {
  # canonical numeric digest helper: round to 1e-12 then digest::digest
  round_d <- function(z) {
    if (is.null(z)) return("NULL")
    z <- as.numeric(z)
    z[!is.finite(z)] <- NA_real_
    z <- signif(z, 12)
    paste(z, collapse = ",")
  }
  if (!is.null(fit$boot$replicates)) {
    # bagged: digest each replicate's coef + selected
    pieces <- vapply(fit$boot$replicates, function(r) {
      paste(round_d(r$coefficients),
            paste(r$selected, collapse = ","),
            round_d(r$gcv),
            sep = "|")
    }, character(1))
    digest_str <- paste(pieces, collapse = ";;")
  } else {
    digest_str <- paste(round_d(fit$coefficients),
                        paste(fit$selected, collapse = ","),
                        round_d(fit$gcv),
                        sep = "|")
  }
  list(
    digest = digest_str,
    coef = fit$coefficients,
    selected = fit$selected,
    gcv = fit$gcv,
    rss = fit$rss
  )
}

run_cell <- function(cell, setup) {
  cell_id <- sprintf("%s__%s__n%d_p%d__%s",
                     cell$family, cell$dgp, cell$n, cell$p, setup$label)
  cat(sprintf("[%s] timing (nthreads=4, reps=%d)... ", cell_id, n_reps))
  t0 <- Sys.time()
  bm <- bench::mark(
    fit = fit_one(cell, setup, nthreads = 4L),
    iterations = n_reps,
    min_iterations = n_reps,
    max_iterations = n_reps,
    check = FALSE,
    filter_gc = FALSE
  )
  wall <- as.numeric(Sys.time() - t0, units = "secs")
  cat(sprintf("median=%.3fs (wall=%.1fs)\n",
              as.numeric(bm$median), wall))

  # Determinism reference fit at nthreads = 1 (one rep)
  cat(sprintf("[%s] nthreads=1 digest... ", cell_id))
  fit1 <- fit_one(cell, setup, nthreads = 1L)
  dig1 <- fit_digest(fit1)
  cat("done\n")

  # Determinism check: nthreads = 4 single fit, must equal nthreads = 1
  fit4 <- fit_one(cell, setup, nthreads = 4L)
  dig4 <- fit_digest(fit4)
  det_ok <- identical(dig1$digest, dig4$digest)
  if (!det_ok) {
    warning(sprintf("DETERMINISM BREAK for %s: nthreads=1 != nthreads=4",
                    cell_id))
  }

  list(
    row = data.frame(
      cell_id = cell_id,
      family = cell$family,
      dgp = cell$dgp,
      n = cell$n,
      p = cell$p,
      setup = setup$label,
      n_reps = n_reps,
      median_s = as.numeric(bm$median),
      min_s = as.numeric(bm$min),
      mem_alloc_b = as.numeric(bm$mem_alloc),
      det_ok = det_ok,
      stringsAsFactors = FALSE
    ),
    digest = dig1
  )
}

# ---- Run -----------------------------------------------------------------
all_cells <- list()
for (ce in cells_def) {
  for (st in setups) {
    all_cells[[length(all_cells) + 1L]] <- list(cell = ce, setup = st)
  }
}

cat(sprintf("=== %s : %d cells, %d reps each ===\n",
            tag, length(all_cells), n_reps))

rows <- vector("list", length(all_cells))
digs <- vector("list", length(all_cells))
names(digs) <- vapply(all_cells, function(s)
  sprintf("%s__%s__n%d_p%d__%s",
          s$cell$family, s$cell$dgp, s$cell$n, s$cell$p, s$setup$label),
  character(1))

for (i in seq_along(all_cells)) {
  s <- all_cells[[i]]
  res <- run_cell(s$cell, s$setup)
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
  base_rds <- file.path(out_dir, "v0.26-speed-baseline.rds")
  if (file.exists(base_rds)) {
    base <- readRDS(base_rds)
    cmp <- merge(base$timings[, c("cell_id", "median_s")],
                 df[, c("cell_id", "median_s")],
                 by = "cell_id", suffixes = c("_base", "_iter"))
    cmp$speedup <- cmp$median_s_base / cmp$median_s_iter
    cmp$pct <- 100 * (cmp$speedup - 1)
    # Digest equivalence
    cmp$digest_ok <- vapply(cmp$cell_id, function(k) {
      isTRUE(identical(base$digests[[k]]$digest, digs[[k]]$digest))
    }, logical(1))
    cmp$coef_linf <- vapply(cmp$cell_id, function(k) {
      a <- base$digests[[k]]$coef
      b <- digs[[k]]$coef
      if (is.null(a) || is.null(b) || length(a) != length(b)) return(Inf)
      max(abs(a - b))
    }, numeric(1))
    cmp$gcv_diff <- vapply(cmp$cell_id, function(k) {
      abs(base$digests[[k]]$gcv - digs[[k]]$gcv)
    }, numeric(1))
    cat("\n=== Comparison vs baseline ===\n")
    print(cmp[order(-cmp$speedup), ], row.names = FALSE)
    gmean <- exp(mean(log(cmp$speedup)))
    cat(sprintf("\nGeometric-mean speedup: %.4fx (%.2f%%)\n",
                gmean, 100 * (gmean - 1)))
    write.csv(cmp, file.path(out_dir, paste0(tag, "-cmp.csv")),
              row.names = FALSE)
  } else {
    cat("(no baseline.rds found; skipping comparison)\n")
  }
}

cat("\nDONE.\n")
