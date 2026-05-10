## Baseline benchmark: ares (v0.0.0.9000) vs earth — pre-v0.1.
## Output: inst/sims/results/baseline-pre-v0.1.csv

suppressPackageStartupMessages({
  library(earth)
})
devtools::load_all(quiet = TRUE)

dgp_friedman <- function(n, seed) {
  set.seed(seed)
  p <- 10
  x <- matrix(runif(n * p), n, p)
  y <- 10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
       10 * x[,4] + 5 * x[,5] + rnorm(n, sd = 1)
  list(x = x, y = y)
}

dgp_additive <- function(n, seed) {
  set.seed(seed)
  p <- 8
  x <- matrix(rnorm(n * p), n, p)
  y <- x[,1] + sin(x[,2]) + abs(x[,3]) + 0.5 * x[,4]^2 + rnorm(n, sd = 0.5)
  list(x = x, y = y)
}

dgp_interaction <- function(n, seed) {
  set.seed(seed)
  p <- 6
  x <- matrix(rnorm(n * p), n, p)
  y <- x[,1] * x[,2] + x[,3] * x[,4] + 0.5 * x[,5] + rnorm(n, sd = 0.7)
  list(x = x, y = y)
}

dgps <- list(friedman = dgp_friedman, additive = dgp_additive, interaction = dgp_interaction)

cells <- expand.grid(
  n = c(500, 1500),
  degree = c(1, 2),
  dgp = c("friedman", "additive", "interaction"),
  threads = c(1L, 2L, 4L),
  stringsAsFactors = FALSE
)

rows <- list()
for (i in seq_len(nrow(cells))) {
  cell <- cells[i, ]
  d <- dgps[[cell$dgp]](cell$n, seed = 100 + i)

  ## ares
  t0 <- Sys.time()
  fit_a <- ares(d$x, d$y, nk = 21, degree = cell$degree, nthreads = cell$threads)
  t_ares <- as.numeric(Sys.time() - t0, units = "secs")
  rss_a <- fit_a$rss
  gcv_a <- fit_a$gcv

  ## earth (single-threaded by definition)
  if (cell$threads == 1L) {
    t0 <- Sys.time()
    fit_e <- earth::earth(d$x, d$y, nk = 21, degree = cell$degree)
    t_earth <- as.numeric(Sys.time() - t0, units = "secs")
    rss_e <- sum((d$y - predict(fit_e))^2)
    gcv_e <- fit_e$gcv
  } else {
    t_earth <- NA_real_; rss_e <- NA_real_; gcv_e <- NA_real_
  }

  rows[[i]] <- data.frame(
    n = cell$n, degree = cell$degree, dgp = cell$dgp, threads = cell$threads,
    t_ares = t_ares, t_earth = t_earth,
    rss_ares = rss_a, rss_earth = rss_e,
    gcv_ares = gcv_a, gcv_earth = gcv_e,
    rss_rel_err = if (!is.na(rss_e) && rss_e > 0) abs(rss_a - rss_e) / rss_e else NA_real_,
    speed_ratio = if (!is.na(t_earth)) t_ares / t_earth else NA_real_
  )
  cat(sprintf("[%2d/%d] %s n=%d deg=%d t=%d  ares=%.3fs  earth=%s  rss_relerr=%s\n",
              i, nrow(cells), cell$dgp, cell$n, cell$degree, cell$threads,
              t_ares,
              if (is.na(t_earth)) "NA" else sprintf("%.3fs", t_earth),
              if (is.na(rows[[i]]$rss_rel_err)) "NA" else sprintf("%.3f%%", 100 * rows[[i]]$rss_rel_err)))
}

out <- do.call(rbind, rows)
dir.create("inst/sims/results", showWarnings = FALSE, recursive = TRUE)
write.csv(out, "inst/sims/results/baseline-pre-v0.1.csv", row.names = FALSE)
cat("\n--- summary at threads=1 ---\n")
sub <- subset(out, threads == 1)
print(summary(sub[, c("t_ares", "t_earth", "rss_rel_err", "speed_ratio")]))
cat("\nWritten: inst/sims/results/baseline-pre-v0.1.csv\n")
