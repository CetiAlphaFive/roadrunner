# Post-process MC results into headline tables.
sim <- utils::read.csv("inst/sims/results/sim_summary.csv", stringsAsFactors = FALSE)
cat("Median rss_match:    ", stats::median(sim$rss_match),  "\n")
cat("Median gcv_ratio:    ", stats::median(sim$gcv_ratio),  "\n")
cat("Median speedup ST:   ", stats::median(sim$speedup_st), "\n")
cat("Median speedup MT:   ", stats::median(sim$speedup_mt), "\n")
cat("Median pred_rmse_diff:", stats::median(sim$pred_rmse), "\n")
