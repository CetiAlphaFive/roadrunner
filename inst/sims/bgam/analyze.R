# Analysis script: produce summary tables and SIM-VAL checklist
# REQ-20260526-013349-mboost-gam
#
# Input:  inst/sims/results/bgam-vs-baselines-0.0.0.9059.rds
# Output: inst/sims/results/bgam-vs-baselines-0.0.0.9059.md

analyze_bgam_sim <- function(rds_path = NULL,
                              md_path  = NULL,
                              pkg_ver  = NULL) {

  if (is.null(rds_path)) {
    rds_path <- "inst/sims/results/bgam-vs-baselines-0.0.0.9059.rds"
  }
  if (is.null(md_path)) {
    md_path <- "inst/sims/results/bgam-vs-baselines-0.0.0.9059.md"
  }

  res <- readRDS(rds_path)

  # ---- helpers -------------------------------------------------------------
  fmt_mean_mcse <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) return("—")
    m  <- mean(x)
    se <- stats::sd(x) / sqrt(length(x))
    sprintf("%.4f (±%.4f)", m, se)
  }

  fmt_prop <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0L) return("—")
    m  <- mean(x)
    se <- sqrt(m * (1 - m) / length(x))
    sprintf("%.3f (±%.3f)", m, se)
  }

  # ---- aggregate by (dgp, cell_id, n, p, sigma, learner) ------------------
  summarize_cell <- function(sub) {
    list(
      n_reps   = nrow(sub),
      n_fail   = sum(sub$failed),
      fail_pct = mean(sub$failed) * 100,
      rmse_fmt     = fmt_mean_mcse(sub$rmse[!sub$failed]),
      bias_fmt     = fmt_mean_mcse(sub$bias[!sub$failed]),
      tp_fmt       = fmt_mean_mcse(sub$tp_rate[!sub$failed]),
      cov_fmt      = fmt_prop(sub$coverage_95[!sub$failed]),
      piw_fmt      = fmt_mean_mcse(sub$pi_width[!sub$failed]),
      brier_fmt    = fmt_mean_mcse(sub$brier[!sub$failed]),
      auc_fmt      = fmt_mean_mcse(sub$auc[!sub$failed]),
      cal_fmt      = fmt_mean_mcse(sub$cal_slope[!sub$failed]),
      # raw means for AC checks
      rmse_mean    = mean(sub$rmse[!sub$failed], na.rm = TRUE),
      tp_mean      = mean(sub$tp_rate[!sub$failed & !is.na(sub$tp_rate)], na.rm = TRUE),
      auc_mean     = mean(sub$auc[!sub$failed], na.rm = TRUE),
      cal_mean     = mean(sub$cal_slope[!sub$failed], na.rm = TRUE),
      cov_mean     = mean(sub$coverage_95[!sub$failed], na.rm = TRUE),
      fail_rate    = mean(sub$failed)
    )
  }

  keys <- c("dgp", "cell_id", "n", "p", "sigma", "learner")
  cell_learner <- unique(res[, keys])

  slist <- lapply(seq_len(nrow(cell_learner)), function(i) {
    k   <- cell_learner[i, ]
    sub <- res[res$dgp     == k$dgp     &
               res$cell_id == k$cell_id &
               res$learner == k$learner, ]
    c(as.list(k), summarize_cell(sub))
  })

  # ---- SIM-VAL checks ------------------------------------------------------
  simval <- list()

  # SIM-VAL-01: failure_rate < 0.05 for all (learner, dgp)
  agg_fail <- stats::aggregate(failed ~ learner + dgp, res, mean)
  sv01 <- all(agg_fail$failed < 0.05)
  simval[["SIM-VAL-01"]] <- list(
    pass = sv01,
    note = if (sv01) {
      "All failure rates < 5%."
    } else {
      bad <- agg_fail[agg_fail$failed >= 0.05, ]
      paste("FAIL:", paste(sprintf("%s/%s=%.1f%%", bad$learner, bad$dgp,
                                  bad$failed * 100), collapse = "; "))
    }
  )

  # SIM-VAL-02: bgam RMSE < 1.20 * best at n=1000, DGP1
  dgp1_1000 <- res[res$dgp == "DGP1" & res$n == 1000L & !res$failed, ]
  if (nrow(dgp1_1000) > 0L && "bgam" %in% dgp1_1000$learner) {
    bgam_r  <- mean(dgp1_1000$rmse[dgp1_1000$learner == "bgam"],  na.rm = TRUE)
    ares_r  <- mean(dgp1_1000$rmse[dgp1_1000$learner == "ares"],  na.rm = TRUE)
    krls_r  <- mean(dgp1_1000$rmse[dgp1_1000$learner == "krls"],  na.rm = TRUE)
    ols_r   <- mean(dgp1_1000$rmse[dgp1_1000$learner == "ols"],   na.rm = TRUE)
    best_r  <- min(ares_r, krls_r, ols_r, na.rm = TRUE)
    ratio   <- bgam_r / best_r
    sv02    <- ratio < 1.20
    simval[["SIM-VAL-02"]] <- list(
      pass = sv02,
      note = sprintf("bgam/best RMSE ratio at n=1000 DGP1 = %.3f (threshold 1.20)", ratio)
    )
  } else {
    simval[["SIM-VAL-02"]] <- list(pass = NA, note = "bgam results not available (blocked)")
  }

  # SIM-VAL-03: bgam TP > 0.40 at n=200, DGP1
  dgp1_200_bgam <- res[res$dgp == "DGP1" & res$n == 200L &
                         res$learner == "bgam" & !res$failed, ]
  if (nrow(dgp1_200_bgam) > 0L) {
    tp_mean <- mean(dgp1_200_bgam$tp_rate, na.rm = TRUE)
    sv03    <- isTRUE(tp_mean > 0.40)
    simval[["SIM-VAL-03"]] <- list(
      pass = sv03,
      note = sprintf("bgam TP rate at n=200 DGP1 = %.3f (threshold 0.40)", tp_mean)
    )
  } else {
    simval[["SIM-VAL-03"]] <- list(pass = NA, note = "bgam results not available (blocked)")
  }

  # SIM-VAL-04: bgam AUC > 0.75 at n=1000, DGP3
  dgp3_1000_bgam <- res[res$dgp == "DGP3" & res$n == 1000L &
                           res$learner == "bgam" & !res$failed, ]
  if (nrow(dgp3_1000_bgam) > 0L) {
    auc_mean <- mean(dgp3_1000_bgam$auc, na.rm = TRUE)
    sv04     <- isTRUE(auc_mean > 0.75)
    simval[["SIM-VAL-04"]] <- list(
      pass = sv04,
      note = sprintf("bgam AUC at n=1000 DGP3 = %.3f (threshold 0.75)", auc_mean)
    )
  } else {
    simval[["SIM-VAL-04"]] <- list(pass = NA, note = "bgam results not available (blocked)")
  }

  # SIM-VAL-05: bgam calibration slope in (0.50, 1.50), DGP3
  dgp3_bgam <- res[res$dgp == "DGP3" & res$learner == "bgam" & !res$failed, ]
  if (nrow(dgp3_bgam) > 0L) {
    cal_mean <- mean(dgp3_bgam$cal_slope, na.rm = TRUE)
    sv05     <- isTRUE(cal_mean > 0.50 && cal_mean < 1.50)
    simval[["SIM-VAL-05"]] <- list(
      pass = sv05,
      note = sprintf("bgam cal slope DGP3 = %.3f (range 0.50–1.50)", cal_mean)
    )
  } else {
    simval[["SIM-VAL-05"]] <- list(pass = NA, note = "bgam results not available (blocked)")
  }

  # SIM-VAL-06: total row count == 16 * 4 * 200 = 12800
  expected <- 16L * 4L * 200L
  sv06     <- nrow(res) == expected
  simval[["SIM-VAL-06"]] <- list(
    pass = sv06,
    note = sprintf("nrow(res) = %d (expected %d)", nrow(res), expected)
  )

  # ---- build markdown ------------------------------------------------------
  lines <- character(0)

  lines <- c(lines,
    sprintf("# bgam vs Baselines — Simulation Results"),
    sprintf("# Package version: %s", if (!is.null(pkg_ver)) pkg_ver else "unknown"),
    sprintf("# Date: %s", Sys.time()),
    sprintf("# RDS: %s", rds_path),
    "",
    "---",
    "",
    "## Simulation Validation Checklist",
    "",
    "| ID | Description | Pass/Fail | Note |",
    "|---|---|---|---|",
    sprintf("| SIM-VAL-01 | Failure rate < 5%% all cells | %s | %s |",
            sv_icon(simval[["SIM-VAL-01"]]$pass), simval[["SIM-VAL-01"]]$note),
    sprintf("| SIM-VAL-02 | bgam RMSE competitive at n=1000 DGP1 | %s | %s |",
            sv_icon(simval[["SIM-VAL-02"]]$pass), simval[["SIM-VAL-02"]]$note),
    sprintf("| SIM-VAL-03 | bgam TP rate > 0.40 at n=200 DGP1 | %s | %s |",
            sv_icon(simval[["SIM-VAL-03"]]$pass), simval[["SIM-VAL-03"]]$note),
    sprintf("| SIM-VAL-04 | bgam AUC > 0.75 at n=1000 DGP3 | %s | %s |",
            sv_icon(simval[["SIM-VAL-04"]]$pass), simval[["SIM-VAL-04"]]$note),
    sprintf("| SIM-VAL-05 | bgam cal slope in (0.50, 1.50) DGP3 | %s | %s |",
            sv_icon(simval[["SIM-VAL-05"]]$pass), simval[["SIM-VAL-05"]]$note),
    sprintf("| SIM-VAL-06 | Row count = 12,800 | %s | %s |",
            sv_icon(simval[["SIM-VAL-06"]]$pass), simval[["SIM-VAL-06"]]$note),
    ""
  )

  # ---- DGP1 table ----------------------------------------------------------
  lines <- c(lines,
    "---",
    "",
    "## DGP1: Sparse Smooth Additive (y = sin(2pi*x1) + 0.5*x2^2 + epsilon)",
    "",
    "| Learner | n | p | sigma | RMSE (+-MCSE) | Bias (+-MCSE) | TP Rate (+-MCSE) | Fail% |",
    "|---------|---|---|-------|--------------|--------------|-----------------|-------|"
  )
  dgp1_sl <- slist[sapply(slist, function(x) x$dgp == "DGP1")]
  dgp1_sl <- dgp1_sl[order(sapply(dgp1_sl, function(x) x$cell_id),
                            sapply(dgp1_sl, function(x) x$learner))]
  for (s in dgp1_sl) {
    lines <- c(lines, sprintf(
      "| %s | %d | %d | %.1f | %s | %s | %s | %.1f%% |",
      s$learner, s$n, s$p, ifelse(is.na(s$sigma), 0, s$sigma),
      s$rmse_fmt, s$bias_fmt,
      ifelse(s$learner == "bgam", s$tp_fmt, "—"),
      s$fail_pct
    ))
  }
  lines <- c(lines, "")

  # ---- DGP2 table ----------------------------------------------------------
  lines <- c(lines,
    "---",
    "",
    "## DGP2: Heteroskedastic (y = sin(2pi*x1) + 0.5*x2^2 + N(0,(1+|x1|)^2*sigma^2))",
    "",
    "| Learner | n | p | sigma | RMSE (+-MCSE) | Coverage 95% (+-MCSE) | PI Width (+-MCSE) | Fail% |",
    "|---------|---|---|-------|--------------|----------------------|------------------|-------|"
  )
  dgp2_sl <- slist[sapply(slist, function(x) x$dgp == "DGP2")]
  dgp2_sl <- dgp2_sl[order(sapply(dgp2_sl, function(x) x$cell_id),
                            sapply(dgp2_sl, function(x) x$learner))]
  for (s in dgp2_sl) {
    cov_col <- if (s$learner %in% c("ares", "ols")) s$cov_fmt else "—"
    piw_col <- if (s$learner %in% c("ares", "ols")) s$piw_fmt else "—"
    lines <- c(lines, sprintf(
      "| %s | %d | %d | %.1f | %s | %s | %s | %.1f%% |",
      s$learner, s$n, s$p, ifelse(is.na(s$sigma), 0, s$sigma),
      s$rmse_fmt, cov_col, piw_col, s$fail_pct
    ))
  }
  lines <- c(lines, "")

  # ---- DGP3 table ----------------------------------------------------------
  lines <- c(lines,
    "---",
    "",
    "## DGP3: Binary Outcome (logit(p) = sin(pi*x1) + x2 - 1)",
    "",
    "| Learner | n | p | Brier (+-MCSE) | AUC (+-MCSE) | Cal Slope (+-MCSE) | Fail% |",
    "|---------|---|---|---------------|-------------|-------------------|-------|"
  )
  dgp3_sl <- slist[sapply(slist, function(x) x$dgp == "DGP3")]
  dgp3_sl <- dgp3_sl[order(sapply(dgp3_sl, function(x) x$cell_id),
                            sapply(dgp3_sl, function(x) x$learner))]
  for (s in dgp3_sl) {
    lines <- c(lines, sprintf(
      "| %s | %d | %d | %s | %s | %s | %.1f%% |",
      s$learner, s$n, s$p,
      s$brier_fmt, s$auc_fmt, s$cal_fmt, s$fail_pct
    ))
  }
  lines <- c(lines, "")

  # ---- AC notes ------------------------------------------------------------
  lines <- c(lines,
    "---",
    "",
    "## Acceptance Criteria Notes",
    "",
    "| AC | Criterion | Evaluation |",
    "|---|---|---|",
    sprintf("| AC-1 | bgam RMSE < 1.20 * best at n=1000 | %s |",
            simval[["SIM-VAL-02"]]$note),
    sprintf("| AC-2 | bgam TP > 0.40 at n=200 DGP1 | %s |",
            simval[["SIM-VAL-03"]]$note),
    sprintf("| AC-3 | Failure rate < 5%% | %s |",
            simval[["SIM-VAL-01"]]$note),
    "| AC-4 | ares + ols PI coverage near 95% | see DGP2 table above |",
    sprintf("| AC-5 | bgam AUC > 0.75 at n=1000 DGP3 | %s |",
            simval[["SIM-VAL-04"]]$note),
    sprintf("| AC-6 | bgam cal slope in (0.50, 1.50) | %s |",
            simval[["SIM-VAL-05"]]$note),
    ""
  )

  writeLines(lines, md_path)
  cat("Wrote", md_path, "\n")
  invisible(simval)
}

# Helper: pass/fail icon (text only — no emojis per project conventions)
sv_icon <- function(x) {
  if (is.na(x)) return("BLOCKED")
  if (isTRUE(x)) return("PASS") else return("FAIL")
}
