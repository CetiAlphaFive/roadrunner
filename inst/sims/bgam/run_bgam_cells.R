# Retry script for bgam cells after builder completes feat/bgam
# REQ-20260526-013349-mboost-gam
#
# Usage: run AFTER 'devtools::load_all()' on feat/bgam branch confirms bgam exists.
#
#   Rscript inst/sims/bgam/run_bgam_cells.R
#
# This script:
#   1. Checks bgam availability
#   2. Loads existing partial results (baseline cells already completed)
#   3. Runs ONLY the bgam learner across all 16 cells x 200 reps
#   4. Merges with baseline results and re-writes the RDS + MD

RDS_PATH  <- "inst/sims/results/bgam-vs-baselines-0.0.0.9059.rds"
RDS_BGAM  <- "inst/sims/results/bgam-cells-0.0.0.9059.rds"
MD_PATH   <- "inst/sims/results/bgam-vs-baselines-0.0.0.9059.md"
SIM_DIR   <- "inst/sims/bgam"

cat("=== bgam cells retry script ===\n")
cat("Date:", format(Sys.time()), "\n")

# ---- check bgam ----
if (!requireNamespace("roadrunner", quietly = TRUE)) {
  stop("roadrunner package not loadable. Run devtools::load_all() first.")
}
if (!exists("bgam", envir = asNamespace("roadrunner"), inherits = FALSE)) {
  stop(paste0(
    "bgam not found in roadrunner namespace. ",
    "Ensure feat/bgam builder has completed and devtools::load_all() was run."
  ))
}
cat("bgam found in roadrunner namespace.\n")

suppressPackageStartupMessages(library(roadrunner))

source(file.path(SIM_DIR, "dgps.R"))
source(file.path(SIM_DIR, "harness.R"))

pkg_ver <- as.character(utils::packageVersion("roadrunner"))
cat("roadrunner version:", pkg_ver, "\n")

N_CORES <- min(2L, parallel::detectCores(logical = FALSE))

cat(sprintf("Running bgam cells: 16 cells x 200 reps on %d core(s)\n", N_CORES))
t0 <- proc.time()[["elapsed"]]

bgam_results <- run_simulation(
  n_reps   = 200L,
  learners = "bgam",
  cores    = N_CORES,
  verbose  = TRUE
)

elapsed <- proc.time()[["elapsed"]] - t0
cat(sprintf("bgam cells done in %.1f minutes.\n", elapsed / 60))

# ---- save bgam-only results ----
saveRDS(bgam_results, RDS_BGAM)
cat("Wrote bgam-only results to", RDS_BGAM, "\n")

# ---- merge with existing baseline results ----
if (file.exists(RDS_PATH)) {
  existing <- readRDS(RDS_PATH)
  # Remove any old bgam rows (placeholder failures) from the existing file
  existing_no_bgam <- existing[existing$learner != "bgam", ]
  combined <- rbind(existing_no_bgam, bgam_results)
  saveRDS(combined, RDS_PATH)
  cat("Merged and updated", RDS_PATH, "\n")
  cat("Total rows:", nrow(combined), "\n")
} else {
  saveRDS(bgam_results, RDS_PATH)
  cat("Wrote combined results to", RDS_PATH, "\n")
}

# ---- regenerate markdown ----
source(file.path(SIM_DIR, "analyze.R"))
analyze_bgam_sim(rds_path = RDS_PATH, md_path = MD_PATH, pkg_ver = pkg_ver)
cat("Done.\n")
