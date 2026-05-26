# Top-level simulation driver: bgam vs ares vs krls vs ols
# REQ-20260526-013349-mboost-gam
#
# Runs all 16 cells x 200 reps x 4 learners = 12,800 fits.
# bgam cells require feat/bgam to be compiled; if bgam is unavailable they are
# recorded as failures and can be rerun with run_bgam_cells.R later.
#
# Usage (from roadrunner repo root):
#   Rscript inst/sims/bgam/run_all.R [cores]
#
# cores: optional integer argument (default = 1). Default is serial because
#        mclapply + TBB (inside krls) is fork-unsafe on Linux at large n/p
#        (n=1000, p=50 will crash the forked worker). Set cores=2 only if
#        you have verified your TBB+fork combination is safe for your krls cells.

args <- commandArgs(trailingOnly = TRUE)
CORES <- if (length(args) >= 1L) as.integer(args[1L]) else 1L
CORES <- max(1L, CORES)

SIM_DIR  <- "inst/sims/bgam"
RDS_PATH <- "inst/sims/results/bgam-vs-baselines-0.0.0.9059.rds"
MD_PATH  <- "inst/sims/results/bgam-vs-baselines-0.0.0.9059.md"

cat("=== bgam vs baselines simulation ===\n")
cat("Date:    ", format(Sys.time()), "\n")
cat("Cores:   ", CORES, "\n")
cat("Reps:     200\n")
cat("Cells:    16 (DGP1: 8, DGP2: 4, DGP3: 4)\n")
cat("Learners: bgam, ares, krls, ols\n\n")

# ---- load package ----
if (!requireNamespace("roadrunner", quietly = TRUE)) {
  stop("roadrunner not loadable. Run devtools::load_all() first.")
}
suppressPackageStartupMessages(library(roadrunner))
pkg_ver <- as.character(utils::packageVersion("roadrunner"))
cat("roadrunner version:", pkg_ver, "\n")

# ---- check bgam availability (non-fatal) ----
bgam_loaded <- tryCatch(
  exists("bgam", envir = asNamespace("roadrunner"), inherits = FALSE),
  error = function(e) FALSE
)
if (!bgam_loaded) {
  cat(paste0(
    "\nNOTE: bgam not found in roadrunner namespace.\n",
    "  bgam cells will be recorded as failures.\n",
    "  After builder completes, run inst/sims/bgam/run_bgam_cells.R\n",
    "  to fill in bgam results and regenerate the summary.\n\n"
  ))
} else {
  cat("bgam found. All learners will be attempted.\n\n")
}

# ---- source simulation files ----
source(file.path(SIM_DIR, "dgps.R"))
source(file.path(SIM_DIR, "harness.R"))
source(file.path(SIM_DIR, "analyze.R"))

# ---- run ----
t0 <- proc.time()[["elapsed"]]

results <- run_simulation(
  n_reps   = 200L,
  learners = c("bgam", "ares", "krls", "ols"),
  cores    = CORES,
  verbose  = TRUE
)

elapsed <- proc.time()[["elapsed"]] - t0
cat(sprintf("\nTotal wall clock: %.1f minutes (%.0f seconds)\n",
            elapsed / 60, elapsed))

# ---- save RDS ----
dir.create(dirname(RDS_PATH), showWarnings = FALSE, recursive = TRUE)
saveRDS(results, RDS_PATH)
cat("Saved RDS:", RDS_PATH, "\n")
cat("Rows:     ", nrow(results), "\n")

# ---- analyze and write markdown ----
simval <- analyze_bgam_sim(rds_path = RDS_PATH, md_path = MD_PATH,
                            pkg_ver = pkg_ver)

# ---- print SIM-VAL summary to console ----
cat("\n=== SIM-VAL CHECKLIST ===\n")
for (nm in names(simval)) {
  status <- if (is.na(simval[[nm]]$pass)) "BLOCKED" else
            if (simval[[nm]]$pass) "PASS" else "FAIL"
  cat(sprintf("  %s: [%s] %s\n", nm, status, simval[[nm]]$note))
}

cat("\nDone. Results in", MD_PATH, "\n")
