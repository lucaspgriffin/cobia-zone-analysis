# =============================================================================
# run_all.R — Master script for Cobia M&R NOAA Zone Analysis Pipeline
# Run from repo root: source("run_all.R")
# =============================================================================

cat("==============================================================\n")
cat("  Cobia Mark-and-Recapture NOAA Zone Analysis Pipeline\n")
cat("==============================================================\n\n")

t_start <- Sys.time()

source("00_setup.R")
source("01_load_data.R")
source("02_assign_zones.R")
source("03_temporal_analysis.R")
source("04_visualize.R")
source("05_sst_analysis.R")

t_end <- Sys.time()
elapsed <- round(difftime(t_end, t_start, units = "secs"), 1)

cat("==============================================================\n")
cat("  Pipeline complete in", elapsed, "seconds\n")
cat("  Outputs in: output/\n")
cat("==============================================================\n")
