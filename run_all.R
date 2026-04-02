# =============================================================================
# run_all.R — Cobia M&R NOAA Zone Analysis Pipeline v2
# Run from repo root: source("run_all.R")
# =============================================================================

cat("==============================================================\n")
cat("  Cobia M&R NOAA Zone Analysis Pipeline v2\n")
cat("  Year-matched SST | Circular stats | GAM framework\n")
cat("==============================================================\n\n")

t_start <- Sys.time()

source("00_setup.R")
source("01_load_data.R")
source("02_assign_zones.R")
source("03_sst_extract.R")
source("04_yearly_peaks.R")
source("05_gam_analysis.R")
source("06_visualize.R")

t_end <- Sys.time()
elapsed <- round(difftime(t_end, t_start, units = "secs"), 1)

cat("==============================================================\n")
cat("  Pipeline complete in", elapsed, "seconds\n")
cat("  Outputs in:", OUTPUT_DIR, "\n")
cat("==============================================================\n")
