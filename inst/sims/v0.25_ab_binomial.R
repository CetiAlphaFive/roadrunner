## v0.25 A/B test: should the binomial default flip to
## auto.linpreds = TRUE, adjust.endspan = 2 ?
##
## Compares ares-default (auto.linpreds=FALSE, adjust.endspan=1) against
## ares-flipped (auto.linpreds=TRUE, adjust.endspan=2) on the four mlbench
## binary DGPs at two sample sizes, with a few repetitions.
##
## Decision rule per the roadmap:
##   - if FLIPPED holdout AUC tie / improve on ALL cells, flip the default;
##   - if any cell regresses, keep the conservative default.
##
## Run: Rscript inst/sims/v0.25_ab_binomial.R [nrep]
## Output: inst/sims/results/v0.25-ab-binomial.csv

suppressPackageStartupMessages({
  library(mlbench)
})
devtools::load_all(quiet = TRUE)

.namex <- function(x) { colnames(x) <- paste0("V", seq_len(ncol(x))); x }
dgps <- list(
  twonorm   = function(n, seed) { set.seed(seed); d <- mlbench.twonorm(n);
                                  list(x = .namex(d$x),
                                       y = as.integer(d$classes) - 1L) },
  threenorm = function(n, seed) { set.seed(seed); d <- mlbench.threenorm(n);
                                  list(x = .namex(d$x),
                                       y = as.integer(d$classes) - 1L) },
  ringnorm  = function(n, seed) { set.seed(seed); d <- mlbench.ringnorm(n);
                                  list(x = .namex(d$x),
                                       y = as.integer(d$classes) - 1L) },
  waveform  = function(n, seed) { set.seed(seed); d <- mlbench.waveform(n);
                                  list(x = .namex(d$x),
                                       y = as.integer(as.integer(d$classes) == 1L)) }
)

auc_simple <- function(prob, y) {
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

nrep <- if (length(commandArgs(TRUE)) >= 1L) as.integer(commandArgs(TRUE)[1]) else 6L
cells <- expand.grid(dgp = names(dgps), n = c(500L, 1500L),
                     rep = seq_len(nrep),
                     stringsAsFactors = FALSE)

rows <- list()
for (i in seq_len(nrow(cells))) {
  ce <- cells[i, ]
  d  <- dgps[[ce$dgp]](ce$n, seed = 2000 + i)
  set.seed(ce$rep * 31L + 11L)
  test_idx <- sample.int(ce$n, max(1L, as.integer(0.30 * ce$n)))
  tr_idx <- setdiff(seq_len(ce$n), test_idx)
  xtr <- d$x[tr_idx,, drop=FALSE]; ytr <- d$y[tr_idx]
  xte <- d$x[test_idx,, drop=FALSE]; yte <- d$y[test_idx]

  # --- A: conservative defaults ---
  fa <- ares(xtr, ytr, family = "binomial", nthreads = 2,
             auto.linpreds = FALSE, adjust.endspan = 1L)
  pa <- predict(fa, xte)
  rows[[length(rows)+1L]] <- data.frame(
    dgp = ce$dgp, n = ce$n, rep = ce$rep, variant = "default",
    auc = auc_simple(pa, yte), logloss = logloss(pa, yte),
    stringsAsFactors = FALSE
  )

  # --- B: flipped (earth-like) ---
  fb <- ares(xtr, ytr, family = "binomial", nthreads = 2,
             auto.linpreds = TRUE, adjust.endspan = 2L)
  pb <- predict(fb, xte)
  rows[[length(rows)+1L]] <- data.frame(
    dgp = ce$dgp, n = ce$n, rep = ce$rep, variant = "flipped",
    auc = auc_simple(pb, yte), logloss = logloss(pb, yte),
    stringsAsFactors = FALSE
  )
  cat(sprintf("[%3d/%3d] %-10s n=%4d rep=%2d  default AUC=%.4f  flipped AUC=%.4f  (delta=%+.4f)\n",
              i, nrow(cells), ce$dgp, ce$n, ce$rep,
              auc_simple(pa, yte), auc_simple(pb, yte),
              auc_simple(pb, yte) - auc_simple(pa, yte)))
}

out <- do.call(rbind, rows)
dir.create("inst/sims/results", showWarnings = FALSE, recursive = TRUE)
write.csv(out, "inst/sims/results/v0.25-ab-binomial.csv", row.names = FALSE)

cat("\n--- mean AUC by (dgp, n, variant) ---\n")
agg <- aggregate(cbind(auc, logloss) ~ dgp + n + variant,
                 data = out,
                 FUN = function(z) round(mean(z, na.rm = TRUE), 4))
print(agg)

cat("\n--- per-cell verdict (flipped - default, AUC) ---\n")
wide <- reshape(out[, c("dgp","n","rep","variant","auc")],
                idvar = c("dgp","n","rep"), timevar = "variant",
                direction = "wide")
wide$delta_auc <- wide$auc.flipped - wide$auc.default
agg2 <- aggregate(delta_auc ~ dgp + n, data = wide,
                  FUN = function(z) round(c(mean = mean(z),
                                             min = min(z),
                                             max = max(z)), 4))
print(agg2)

cat("\nWritten: inst/sims/results/v0.25-ab-binomial.csv\n")
