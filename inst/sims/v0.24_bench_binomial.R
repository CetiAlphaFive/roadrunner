## v0.24 binomial benchmark: ares (family="binomial") vs earth(glm=binomial())
## vs ranger(classification = TRUE) on four mlbench classification DGPs.
##
## DGPs: twonorm, threenorm, ringnorm, waveform (waveform collapsed to
## class 1 vs not-1 so it's binary for Phase A).
## Sample sizes: n = 500 and n = 1500 (train + 30% holdout split).
## Repetitions: nrep = 10 per (dgp, n) cell.
## Metrics: holdout AUC, log-loss, classification accuracy, fit time.
##
## Run with:  Rscript inst/sims/v0.24_bench_binomial.R [nrep]
## Output:    inst/sims/results/v0.24-binomial.csv (long, one row per
##            (dgp, n, rep, model)).

suppressPackageStartupMessages({
  library(earth)
  library(ranger)
  library(mlbench)
})
devtools::load_all(quiet = TRUE)

# ---- DGP wrappers ---------------------------------------------------------
# Each wrapper takes (n, seed) and returns list(x, y) with y in {0, 1}.
.namex <- function(x) { colnames(x) <- paste0("V", seq_len(ncol(x))); x }
dgp_twonorm <- function(n, seed) {
  set.seed(seed)
  d <- mlbench.twonorm(n)
  list(x = .namex(d$x), y = as.integer(d$classes) - 1L)
}
dgp_threenorm <- function(n, seed) {
  set.seed(seed)
  d <- mlbench.threenorm(n)
  list(x = .namex(d$x), y = as.integer(d$classes) - 1L)
}
dgp_ringnorm <- function(n, seed) {
  set.seed(seed)
  d <- mlbench.ringnorm(n)
  list(x = .namex(d$x), y = as.integer(d$classes) - 1L)
}
dgp_waveform <- function(n, seed) {
  set.seed(seed)
  d <- mlbench.waveform(n)
  # 3-class -> 2-class (class 1 vs not-1).
  list(x = .namex(d$x), y = as.integer(as.integer(d$classes) == 1L))
}
dgps <- list(twonorm = dgp_twonorm,
             threenorm = dgp_threenorm,
             ringnorm = dgp_ringnorm,
             waveform = dgp_waveform)

# ---- Metrics --------------------------------------------------------------
auc_simple <- function(prob, y) {
  # Mann-Whitney based AUC. y in {0,1}.
  if (length(unique(y)) < 2L) return(NA_real_)
  r <- rank(prob, ties.method = "average")
  n1 <- sum(y == 1L); n0 <- sum(y == 0L)
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n0 * n1)
}
logloss <- function(prob, y) {
  eps <- 1e-15
  p <- pmin(pmax(prob, eps), 1 - eps)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}
acc <- function(prob, y) mean((prob >= 0.5) == (y == 1L))

# ---- Cell runner ----------------------------------------------------------
nrep <- if (length(commandArgs(TRUE)) >= 1L) as.integer(commandArgs(TRUE)[1]) else 10L
cells <- expand.grid(dgp = names(dgps), n = c(500L, 1500L),
                     rep = seq_len(nrep),
                     stringsAsFactors = FALSE)

rows <- list()
for (i in seq_len(nrow(cells))) {
  ce <- cells[i, ]
  d  <- dgps[[ce$dgp]](ce$n, seed = 1000 + i)
  # 70/30 split
  set.seed(ce$rep * 31L + 7L)
  test_idx <- sample.int(ce$n, max(1L, as.integer(0.30 * ce$n)))
  train_idx <- setdiff(seq_len(ce$n), test_idx)
  xtr <- d$x[train_idx, , drop = FALSE]; ytr <- d$y[train_idx]
  xte <- d$x[test_idx,  , drop = FALSE]; yte <- d$y[test_idx]

  # --- ares (family = binomial) ---
  t0 <- Sys.time()
  fit_a <- ares(xtr, ytr, family = "binomial", nthreads = 2)
  t_a <- as.numeric(Sys.time() - t0, units = "secs")
  p_a <- predict(fit_a, xte)        # response scale
  rows[[length(rows) + 1L]] <- data.frame(
    dgp = ce$dgp, n = ce$n, rep = ce$rep, model = "ares",
    auc = auc_simple(p_a, yte),
    logloss = logloss(p_a, yte),
    acc = acc(p_a, yte),
    t_fit = t_a, stringsAsFactors = FALSE
  )

  # --- earth (glm = binomial) ---
  t0 <- Sys.time()
  fit_e <- tryCatch(
    earth::earth(xtr, ytr, glm = list(family = stats::binomial())),
    error = function(e) NULL
  )
  t_e <- as.numeric(Sys.time() - t0, units = "secs")
  if (!is.null(fit_e)) {
    p_e <- as.numeric(predict(fit_e, xte, type = "response"))
    rows[[length(rows) + 1L]] <- data.frame(
      dgp = ce$dgp, n = ce$n, rep = ce$rep, model = "earth",
      auc = auc_simple(p_e, yte),
      logloss = logloss(p_e, yte),
      acc = acc(p_e, yte),
      t_fit = t_e, stringsAsFactors = FALSE
    )
  } else {
    rows[[length(rows) + 1L]] <- data.frame(
      dgp = ce$dgp, n = ce$n, rep = ce$rep, model = "earth",
      auc = NA_real_, logloss = NA_real_, acc = NA_real_,
      t_fit = NA_real_, stringsAsFactors = FALSE
    )
  }

  # --- ranger (classification) ---
  ytr_f <- factor(ytr, levels = c(0L, 1L))
  t0 <- Sys.time()
  fit_r <- ranger::ranger(x = xtr, y = ytr_f, num.trees = 500,
                          probability = TRUE, num.threads = 2,
                          seed = ce$rep + 1L)
  t_r <- as.numeric(Sys.time() - t0, units = "secs")
  p_r <- predict(fit_r, data = xte)$predictions[, "1"]
  rows[[length(rows) + 1L]] <- data.frame(
    dgp = ce$dgp, n = ce$n, rep = ce$rep, model = "ranger",
    auc = auc_simple(p_r, yte),
    logloss = logloss(p_r, yte),
    acc = acc(p_r, yte),
    t_fit = t_r, stringsAsFactors = FALSE
  )

  cat(sprintf("[%3d/%3d] %-10s n=%4d rep=%2d  ares=%.3fs/AUC=%.3f  earth=%.3fs/AUC=%.3f  ranger=%.3fs/AUC=%.3f\n",
              i, nrow(cells), ce$dgp, ce$n, ce$rep,
              t_a, auc_simple(p_a, yte),
              t_e, if (!is.null(fit_e)) auc_simple(p_e, yte) else NA_real_,
              t_r, auc_simple(p_r, yte)))
}

out <- do.call(rbind, rows)
dir.create("inst/sims/results", showWarnings = FALSE, recursive = TRUE)
write.csv(out, "inst/sims/results/v0.24-binomial.csv", row.names = FALSE)

# ---- Summary ----
cat("\n--- mean +/- sd by (dgp, n, model) ---\n")
agg <- aggregate(cbind(auc, logloss, acc, t_fit) ~ dgp + n + model,
                 data = out, FUN = function(z) c(mean = mean(z, na.rm = TRUE),
                                                  sd = sd(z, na.rm = TRUE)))
print(agg)
cat("\nWritten: inst/sims/results/v0.24-binomial.csv\n")
