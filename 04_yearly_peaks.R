# =============================================================================
# 04_yearly_peaks.R — Year-specific peak detection with circular statistics
# Per zone per year: find peak week, SST at peak, daylength at peak
# Detect bimodal zones using Hartigan's dip test
# =============================================================================

cat("Computing yearly peaks with circular statistics...\n")

if (!requireNamespace("diptest", quietly = TRUE)) {
  cat("  Installing diptest package...\n")
  install.packages("diptest", repos = "https://cloud.r-project.org", quiet = TRUE)
}
library(diptest)

# Use only OISST-era data with assigned zones and adequate sample
zoned <- mr_oisst %>% filter(!is.na(tag_zone) & !is.na(tag_week))

MIN_TAGS_YEAR <- 5   # minimum tags per zone-year to estimate a peak
MIN_TAGS_ZONE <- 30  # minimum total tags per zone to include in analysis

# --- Identify zones with enough data -----------------------------------------
zone_n <- zoned %>% count(tag_zone, name = "n_total")
valid_zones <- zone_n %>% filter(n_total >= MIN_TAGS_ZONE) %>% pull(tag_zone)
cat("  Zones with >=", MIN_TAGS_ZONE, "tags:", paste(sort(valid_zones), collapse = ", "), "\n")

# --- Circular statistics per zone (pooled across years for diagnostics) -------
circ_summary <- zoned %>%
  filter(tag_zone %in% valid_zones) %>%
  group_by(tag_zone) %>%
  summarise(
    n = n(),
    # Convert Julian day to circular (radians)
    circ_mean_rad = {
      theta <- circular(tag_julian * 2 * pi / 365, type = "angles",
                        units = "radians", rotation = "clock")
      as.numeric(mean.circular(theta))
    },
    # Mean resultant length (concentration: 0=uniform, 1=perfectly concentrated)
    rho_bar = {
      theta <- circular(tag_julian * 2 * pi / 365, type = "angles",
                        units = "radians", rotation = "clock")
      rho.circular(theta)
    },
    # Hartigan dip test for multimodality (on Julian day, not circular)
    dip_stat = dip.test(tag_julian)$statistic,
    dip_p    = dip.test(tag_julian)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    circ_mean_julian = (circ_mean_rad * 365 / (2 * pi)) %% 365,
    is_bimodal = dip_p < 0.05
  )

cat("\n  CIRCULAR STATISTICS (pooled across years):\n")
cat(sprintf("  %4s %5s %8s %6s %8s %s\n",
            "Zone", "n", "CircMean", "Rho", "DipP", "Bimodal?"))
for (i in seq_len(nrow(circ_summary))) {
  cs <- circ_summary[i, ]
  cat(sprintf("  %4d %5d %8.0f %6.3f %8.4f %s\n",
              cs$tag_zone, cs$n, cs$circ_mean_julian, cs$rho_bar,
              cs$dip_p, ifelse(cs$is_bimodal, "YES", "")))
}

bimodal_zones <- circ_summary %>% filter(is_bimodal) %>% pull(tag_zone)
if (length(bimodal_zones) > 0) {
  cat("  Bimodal zones (dip test p<0.05):", paste(bimodal_zones, collapse = ", "), "\n")
}

# --- Yearly peaks per zone ----------------------------------------------------
# For each zone-year with enough data, find the peak week and associated SST
yearly_peaks <- zoned %>%
  filter(tag_zone %in% valid_zones) %>%
  group_by(tag_zone, tag_year) %>%
  filter(n() >= MIN_TAGS_YEAR) %>%
  summarise(
    n_tags     = n(),
    # Peak week = mode of tag_week
    peak_week  = {
      wk_counts <- table(tag_week)
      as.integer(names(which.max(wk_counts)))
    },
    # Circular mean Julian day
    circ_mean_julian = {
      theta <- circular(tag_julian * 2 * pi / 365, type = "angles",
                        units = "radians", rotation = "clock")
      (as.numeric(mean.circular(theta)) * 365 / (2 * pi)) %% 365
    },
    # Median Julian day (for comparison)
    median_julian = median(tag_julian),
    # Mean SST during the peak week (year-matched!)
    peak_sst = {
      pk_wk <- as.integer(names(which.max(table(tag_week))))
      mean(sst_mean[tag_week == pk_wk], na.rm = TRUE)
    },
    # Mean SST across ALL tags that year-zone (observed temp)
    mean_sst    = mean(sst_mean, na.rm = TRUE),
    # Mean daylength at tagging
    mean_daylength = mean(daylength_hrs, na.rm = TRUE),
    # Daylength at peak week
    peak_daylength = {
      pk_wk <- as.integer(names(which.max(table(tag_week))))
      mean(daylength_hrs[tag_week == pk_wk], na.rm = TRUE)
    },
    .groups = "drop"
  )

cat("\n  Yearly peaks:", nrow(yearly_peaks), "zone-year combinations\n")
cat("  Years covered:", min(yearly_peaks$tag_year), "-", max(yearly_peaks$tag_year), "\n")
cat("  Zones covered:", length(unique(yearly_peaks$tag_zone)), "\n")

# --- Summary: Mean peak SST and peak week per zone (across years) -------------
zone_peak_summary <- yearly_peaks %>%
  group_by(tag_zone) %>%
  summarise(
    n_years    = n(),
    mean_peak_week = mean(peak_week, na.rm = TRUE),
    sd_peak_week   = sd(peak_week, na.rm = TRUE),
    mean_peak_sst  = mean(peak_sst, na.rm = TRUE),
    sd_peak_sst    = sd(peak_sst, na.rm = TRUE),
    mean_peak_daylength = mean(peak_daylength, na.rm = TRUE),
    sd_peak_daylength   = sd(peak_daylength, na.rm = TRUE),
    mean_obs_sst   = mean(mean_sst, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n  ZONE PEAK SUMMARY (averaged across years):\n")
cat(sprintf("  %4s %4s %8s %7s %8s %7s %8s\n",
            "Zone", "nYr", "PkWk", "sdWk", "PkSST", "sdSST", "PkDayL"))
for (i in seq_len(nrow(zone_peak_summary))) {
  z <- zone_peak_summary[i, ]
  cat(sprintf("  %4d %4d %8.1f %7.1f %7.1fC %6.1fC %7.1fh\n",
              z$tag_zone, z$n_years, z$mean_peak_week, z$sd_peak_week,
              z$mean_peak_sst, z$sd_peak_sst, z$mean_peak_daylength))
}

# --- CV analysis: SST vs Julian day vs daylength (across zone-years) ----------
cat("\n  CV ANALYSIS (what's more consistent at peaks across zone-years?):\n")

valid_yp <- yearly_peaks %>% filter(!is.na(peak_sst) & !is.na(peak_daylength))

cv_sst <- sd(valid_yp$peak_sst, na.rm = TRUE) / mean(valid_yp$peak_sst, na.rm = TRUE)
cv_jul <- sd(valid_yp$circ_mean_julian, na.rm = TRUE) / mean(valid_yp$circ_mean_julian, na.rm = TRUE)
cv_day <- sd(valid_yp$peak_daylength, na.rm = TRUE) / mean(valid_yp$peak_daylength, na.rm = TRUE)

cat("    Peak SST:       CV =", round(cv_sst, 4), "\n")
cat("    Julian day:     CV =", round(cv_jul, 4), "\n")
cat("    Daylength:      CV =", round(cv_day, 4), "\n")
cat("    Lowest CV = most consistent across zone-years = likely driver\n")

# Save outputs
write.csv(yearly_peaks, file.path(OUTPUT_DIR, "yearly_peaks.csv"), row.names = FALSE)
write.csv(zone_peak_summary, file.path(OUTPUT_DIR, "zone_peak_summary.csv"), row.names = FALSE)
write.csv(circ_summary, file.path(OUTPUT_DIR, "circular_summary.csv"), row.names = FALSE)

cat("\n04_yearly_peaks.R complete.\n\n")
