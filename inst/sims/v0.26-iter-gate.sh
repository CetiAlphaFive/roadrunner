#!/usr/bin/env bash
# v0.26 iteration gate runner.
# Usage: ./inst/sims/v0.26-iter-gate.sh <iter_num>
# Runs the per-iteration gate: rebuild + tests + R CMD check + baseline rerun.
# Exits non-zero on any gate failure.
set -euo pipefail
ITER="${1:?usage: v0.26-iter-gate.sh <iter_num>}"
cd /home/jack/Dropbox/roadrunner

GATE_LOG="/tmp/v0.26-iter-${ITER}-gate.log"
echo "=== v0.26 iter-${ITER} gate ===" | tee "$GATE_LOG"

echo "[1/4] Rebuild + ARES_FULL_TESTS=1 testthat..." | tee -a "$GATE_LOG"
ARES_FULL_TESTS=1 Rscript -e '
  devtools::clean_dll(quiet=TRUE)
  devtools::load_all(quiet=TRUE)
  res <- testthat::test_local(reporter = "summary", stop_on_failure = FALSE)
  df <- as.data.frame(res)
  cat(sprintf("totals: pass=%d fail=%d warn=%d skip=%d\n",
              sum(df$passed), sum(df$failed), sum(df$warning), sum(df$skipped)))
  if (sum(df$failed) > 0L) {
    cat("FAILED tests:\n")
    print(df[df$failed > 0L, c("file", "test", "failed")])
    quit(status = 1L)
  }
' 2>&1 | tee -a "$GATE_LOG"

echo "[2/4] R CMD check (no manual, no vignettes)..." | tee -a "$GATE_LOG"
Rscript -e '
  res <- rcmdcheck::rcmdcheck(args = c("--no-manual", "--no-build-vignettes"),
                              quiet = TRUE,
                              error_on = "warning")
  if (length(res$errors) > 0 || length(res$warnings) > 0) {
    cat("R CMD check failures:\n")
    print(res$errors); print(res$warnings)
    quit(status = 1L)
  }
  cat(sprintf("R CMD check: 0 errors, 0 warnings, %d notes\n",
              length(res$notes)))
' 2>&1 | tee -a "$GATE_LOG"

echo "[3/4] Baseline timing rerun (iter=${ITER})..." | tee -a "$GATE_LOG"
Rscript inst/sims/v0.26-speed-baseline.R "$ITER" 2>&1 | tee -a "$GATE_LOG"

echo "[4/4] Comparison vs baseline.rds..." | tee -a "$GATE_LOG"
Rscript -e "
  out_dir <- 'inst/sims/results'
  base <- readRDS(file.path(out_dir, 'v0.26-speed-baseline.rds'))
  iter <- readRDS(file.path(out_dir, sprintf('v0.26-iter-%02d.rds', ${ITER})))
  cmp <- merge(base\$timings[, c('cell_id', 'median_s')],
               iter\$timings[, c('cell_id', 'median_s')],
               by = 'cell_id', suffixes = c('_base', '_iter'))
  cmp\$speedup <- cmp\$median_s_base / cmp\$median_s_iter
  cmp\$digest_ok <- vapply(cmp\$cell_id, function(k)
    isTRUE(identical(base\$digests[[k]]\$digest, iter\$digests[[k]]\$digest)),
    logical(1))
  cmp\$coef_linf <- vapply(cmp\$cell_id, function(k) {
    a <- base\$digests[[k]]\$coef; b <- iter\$digests[[k]]\$coef
    if (is.null(a) || is.null(b) || length(a) != length(b)) return(Inf)
    max(abs(a - b))
  }, numeric(1))
  cmp\$gcv_diff <- vapply(cmp\$cell_id, function(k)
    abs(base\$digests[[k]]\$gcv - iter\$digests[[k]]\$gcv),
    numeric(1))
  print(cmp[order(-cmp\$speedup), ], row.names = FALSE)
  gmean <- exp(mean(log(cmp\$speedup)))
  cat(sprintf('\nGeometric-mean speedup: %.4fx (%.2f%%)\n', gmean, 100*(gmean-1)))
  cat(sprintf('Digests identical: %d / %d cells\n', sum(cmp\$digest_ok), nrow(cmp)))
  cat(sprintf('Max coef L_inf diff: %.3e\n', max(cmp\$coef_linf)))
  cat(sprintf('Max GCV diff: %.3e\n', max(cmp\$gcv_diff)))
  # Gate
  ok_det <- all(cmp\$digest_ok | (cmp\$coef_linf <= 1e-10 & cmp\$gcv_diff <= 1e-12))
  ok_speed <- gmean >= 1.0
  if (!ok_det) { cat('GATE FAIL: fit equivalence broken\n'); quit(status=1L) }
  if (!ok_speed) { cat('GATE FAIL: net regression\n'); quit(status=1L) }
  cat('GATE PASS\n')
" 2>&1 | tee -a "$GATE_LOG"

echo "=== iter-${ITER} gate PASSED ==="
