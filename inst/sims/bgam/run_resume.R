# Resume script: runs remaining cells after crash at DGP1 cell 8
# REQ-20260526-013349-mboost-gam
#
# Context: mclapply + TBB (inside krls) is fork-unsafe at n=1000, p=50.
# This script runs ALL cells serially (cores=1) to avoid the fork/TBB crash.
# Because we have no partial RDS from the crashed run, we re-run all 16 cells.
# (DGP1 cells 1-7 results were held in memory and lost when the process crashed.)
#
# Usage (from repo root):
#   Rscript inst/sims/bgam/run_resume.R

SIM_DIR  <- "inst/sims/bgam"
RDS_PATH <- "inst/sims/results/bgam-vs-baselines-0.0.0.9059.rds"
MD_PATH  <- "inst/sims/results/bgam-vs-baselines-0.0.0.9059.md"

cat("=== bgam vs baselines simulation (serial resume) ===\n")
cat("Date:    ", format(Sys.time()), "\n")
cat("Cores:    1 (serial — avoids TBB/fork crash at n=1000, p=50)\n")
cat("Reps:     200\n\n")

if (!requireNamespace("roadrunner", quietly = TRUE)) {
  stop("roadrunner not loadable.")
}
suppressPackageStartupMessages(library(roadrunner))
pkg_ver <- as.character(utils::packageVersion("roadrunner"))
cat("roadrunner version:", pkg_ver, "\n\n")

source(file.path(SIM_DIR, "dgps.R"))
source(file.path(SIM_DIR, "harness.R"))
source(file.path(SIM_DIR, "analyze.R"))

t0 <- proc.time()[["elapsed"]]

results <- run_simulation(
  n_reps   = 200L,
  learners = c("bgam", "ares", "krls", "ols"),
  cores    = 1L,
  verbose  = TRUE
)

elapsed <- proc.time()[["elapsed"]] - t0
cat(sprintf("\nTotal wall clock: %.1f minutes (%.0f seconds)\n",
            elapsed / 60, elapsed))

dir.create(dirname(RDS_PATH), showWarnings = FALSE, recursive = TRUE)
saveRDS(results, RDS_PATH)
cat("Saved RDS:", RDS_PATH, "\n")
cat("Rows:     ", nrow(results), "\n")

simval <- analyze_bgam_sim(rds_path = RDS_PATH, md_path = MD_PATH, pkg_ver = pkg_ver)

cat("\n=== SIM-VAL CHECKLIST ===\n")
for (nm in names(simval)) {
  status <- if (is.na(simval[[nm]]$pass)) "BLOCKED" else
            if (simval[[nm]]$pass) "PASS" else "FAIL"
  cat(sprintf("  %s: [%s] %s\n", nm, status, simval[[nm]]$note))
}
cat("\nDone. Results in", MD_PATH, "\n")
