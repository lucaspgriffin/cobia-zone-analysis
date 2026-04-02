# =============================================================================
# 03_temporal_analysis.R — Weekly aggregation, peak detection, movement matrix
# Cobia M&R NOAA Zone Analysis Pipeline
# =============================================================================

cat("Running temporal analysis...\n")

# Only use records with assigned zones
zoned_data <- mr_data %>% filter(!is.na(tag_zone))

# --- Weekly catch by zone (pooled across years) ------------------------------
weekly_catch <- zoned_data %>%
  group_by(tag_zone, tag_week) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(tag_zone = 1:21, tag_week = 1:52, fill = list(n = 0))

cat("  Weekly catch table:", nrow(weekly_catch), "rows (zone x week)\n")

# --- Weekly catch by zone and year (for diagnostics) -------------------------
weekly_catch_yearly <- zoned_data %>%
  group_by(tag_zone, tag_year, tag_week) %>%
  summarise(n = n(), .groups = "drop")

# --- Peak weeks per zone -----------------------------------------------------
# Top-3 weeks and weeks above 75th percentile
peak_weeks <- weekly_catch %>%
  filter(n > 0) %>%
  group_by(tag_zone) %>%
  mutate(
    q75 = quantile(n, 0.75),
    above_q75 = n >= q75
  ) %>%
  arrange(tag_zone, desc(n)) %>%
  mutate(rank = row_number()) %>%
  filter(rank <= 3 | above_q75) %>%
  ungroup()

cat("  Peak weeks identified:", nrow(peak_weeks), "zone-week combos\n")

# --- Zone summary ------------------------------------------------------------
zone_summary <- zoned_data %>%
  group_by(tag_zone) %>%
  summarise(
    n_tags         = n(),
    n_recaptured   = sum(recaptured),
    recap_rate     = round(mean(recaptured) * 100, 1),
    median_julian  = median(tag_julian, na.rm = TRUE),
    iqr_julian     = IQR(tag_julian, na.rm = TRUE),
    peak_week      = {
      wk_counts <- table(tag_week)
      as.integer(names(which.max(wk_counts)))
    },
    peak_month     = {
      mo_counts <- table(tag_month)
      as.integer(names(which.max(mo_counts)))
    },
    mean_length_cm = round(mean(length_cm, na.rm = TRUE), 1),
    min_year       = min(tag_year, na.rm = TRUE),
    max_year       = max(tag_year, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(tag_zone)

cat("\n  ZONE SUMMARY:\n")
cat("  ", paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("  %4s %6s %5s %6s %7s %6s %6s %8s\n",
            "Zone", "nTags", "nRec", "Rec%", "MedJul", "PkWk", "PkMo", "MeanLen"))
cat("  ", paste(rep("-", 80), collapse = ""), "\n")
for (i in seq_len(nrow(zone_summary))) {
  z <- zone_summary[i, ]
  cat(sprintf("  %4d %6d %5d %5.1f%% %7.0f %6d %6d %7.1f cm\n",
              z$tag_zone, z$n_tags, z$n_recaptured, z$recap_rate,
              z$median_julian, z$peak_week, z$peak_month, z$mean_length_cm))
}

# Flag low-sample zones
low_n_zones <- zone_summary %>% filter(n_tags < 10) %>% pull(tag_zone)
if (length(low_n_zones) > 0) {
  cat("\n  WARNING: Zones with < 10 tags (excluded from peak analysis):",
      paste(low_n_zones, collapse = ", "), "\n")
}

# --- Movement matrix (tag zone x recap zone) ---------------------------------
recap_zoned <- mr_data %>%
  filter(recaptured & !is.na(tag_zone) & !is.na(recap_zone))

if (nrow(recap_zoned) > 0) {
  movement_matrix <- table(
    Tag_Zone = recap_zoned$tag_zone,
    Recap_Zone = recap_zoned$recap_zone
  )

  # Convert to long format for plotting
  movement_long <- as.data.frame(movement_matrix) %>%
    mutate(
      Tag_Zone = as.integer(as.character(Tag_Zone)),
      Recap_Zone = as.integer(as.character(Recap_Zone))
    )

  cat("\n  Movement matrix:", nrow(recap_zoned), "recaptured fish with both zones\n")
  cat("  Same-zone recaptures:", sum(recap_zoned$tag_zone == recap_zoned$recap_zone), "\n")
  cat("  Cross-zone recaptures:", sum(recap_zoned$tag_zone != recap_zoned$recap_zone), "\n")
} else {
  movement_long <- data.frame(Tag_Zone = integer(), Recap_Zone = integer(), Freq = integer())
  cat("\n  No recaptured fish with both tag and recap zones assigned.\n")
}

# --- Peak Julian day by zone (for map coloring) ------------------------------
peak_by_zone <- zone_summary %>%
  select(tag_zone, median_julian, peak_week, peak_month, n_tags)

# Save summary table
write.csv(zone_summary, file.path(OUTPUT_DIR, "zone_summary.csv"), row.names = FALSE)
cat("  Zone summary saved to", file.path(OUTPUT_DIR, "zone_summary.csv"), "\n")

cat("\n03_temporal_analysis.R complete.\n\n")
