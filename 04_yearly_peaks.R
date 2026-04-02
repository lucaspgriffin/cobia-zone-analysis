# =============================================================================
# 04_yearly_peaks.R — Season-split peak detection with circular statistics
#
# Approach: Split at summer solstice (Julian 172 / ~week 25)
# - Bimodal zones (dip test p < 0.05): compute spring AND fall peaks
# - Unimodal zones: assign all records to their dominant season only
# - Each zone-year-season = one observation for GAM analysis
# =============================================================================

cat("Computing seasonal peaks (spring/fall split)...\n")

if (!requireNamespace("diptest", quietly = TRUE)) {
  install.packages("diptest", repos = "https://cloud.r-project.org", quiet = TRUE)
}
library(diptest)

SOLSTICE_JULIAN <- 172  # June 21 = summer solstice boundary
MIN_TAGS_SEASON <- 5    # minimum tags per zone-year-season to estimate a peak
MIN_TAGS_ZONE   <- 30   # minimum total tags per zone

zoned <- mr_oisst %>% filter(!is.na(tag_zone) & !is.na(tag_week))

# --- Identify valid zones ----------------------------------------------------
zone_n <- zoned %>% count(tag_zone, name = "n_total")
valid_zones <- zone_n %>% filter(n_total >= MIN_TAGS_ZONE) %>% pull(tag_zone)
cat("  Zones with >=", MIN_TAGS_ZONE, "tags:", paste(sort(valid_zones), collapse = ", "), "\n")

# --- Bimodality classification per zone (pooled) -----------------------------
circ_summary <- zoned %>%
  filter(tag_zone %in% valid_zones) %>%
  group_by(tag_zone) %>%
  summarise(
    n = n(),
    circ_mean_rad = {
      theta <- circular(tag_julian * 2 * pi / 365, type = "angles",
                        units = "radians", rotation = "clock")
      as.numeric(mean.circular(theta))
    },
    rho_bar = {
      theta <- circular(tag_julian * 2 * pi / 365, type = "angles",
                        units = "radians", rotation = "clock")
      rho.circular(theta)
    },
    dip_stat = dip.test(tag_julian)$statistic,
    dip_p    = dip.test(tag_julian)$p.value,
    # Proportion in spring vs fall
    pct_spring = mean(tag_julian < SOLSTICE_JULIAN) * 100,
    pct_fall   = mean(tag_julian >= SOLSTICE_JULIAN) * 100,
    .groups = "drop"
  ) %>%
  mutate(
    circ_mean_julian = (circ_mean_rad * 365 / (2 * pi)) %% 365,
    is_bimodal = dip_p < 0.05,
    # For unimodal zones: which season dominates?
    dominant_season = ifelse(circ_mean_julian < SOLSTICE_JULIAN, "spring", "fall"),
    # Classification
    modality = case_when(
      is_bimodal ~ "bimodal",
      TRUE ~ paste0("unimodal_", dominant_season)
    )
  )

cat("\n  ZONE MODALITY CLASSIFICATION:\n")
cat(sprintf("  %4s %5s %8s %6s %6s %6s %8s %s\n",
            "Zone", "n", "CircMn", "Rho", "Spr%", "Fal%", "DipP", "Class"))
for (i in seq_len(nrow(circ_summary))) {
  cs <- circ_summary[i, ]
  cat(sprintf("  %4d %5d %8.0f %6.3f %5.1f%% %5.1f%% %8.4f %s\n",
              cs$tag_zone, cs$n, cs$circ_mean_julian, cs$rho_bar,
              cs$pct_spring, cs$pct_fall, cs$dip_p, cs$modality))
}

bimodal_zones <- circ_summary %>% filter(is_bimodal) %>% pull(tag_zone)
unimodal_spring <- circ_summary %>% filter(modality == "unimodal_spring") %>% pull(tag_zone)
unimodal_fall <- circ_summary %>% filter(modality == "unimodal_fall") %>% pull(tag_zone)

cat("\n  Bimodal:", paste(bimodal_zones, collapse = ", "), "\n")
if (length(unimodal_spring) > 0) cat("  Unimodal spring:", paste(unimodal_spring, collapse = ", "), "\n")
if (length(unimodal_fall) > 0) cat("  Unimodal fall:", paste(unimodal_fall, collapse = ", "), "\n")

# --- Assign season to each tag record ----------------------------------------
zone_modality <- circ_summary %>% select(tag_zone, modality, dominant_season)

zoned <- zoned %>%
  left_join(zone_modality, by = "tag_zone") %>%
  mutate(
    season = case_when(
      # Bimodal zones: split at solstice
      modality == "bimodal" & tag_julian < SOLSTICE_JULIAN ~ "spring",
      modality == "bimodal" & tag_julian >= SOLSTICE_JULIAN ~ "fall",
      # Unimodal zones: all records get the dominant season
      TRUE ~ dominant_season
    )
  )

cat("  Spring records:", sum(zoned$season == "spring", na.rm = TRUE), "\n")
cat("  Fall records:", sum(zoned$season == "fall", na.rm = TRUE), "\n")

# --- Compute seasonal peaks per zone-year ------------------------------------
seasonal_peaks <- zoned %>%
  filter(tag_zone %in% valid_zones) %>%
  group_by(tag_zone, tag_year, season) %>%
  filter(n() >= MIN_TAGS_SEASON) %>%
  summarise(
    n_tags = n(),
    modality = first(modality),
    peak_week = {
      wk_counts <- table(tag_week)
      as.integer(names(which.max(wk_counts)))
    },
    circ_mean_julian = {
      theta <- circular(tag_julian * 2 * pi / 365, type = "angles",
                        units = "radians", rotation = "clock")
      (as.numeric(mean.circular(theta)) * 365 / (2 * pi)) %% 365
    },
    median_julian = median(tag_julian),
    peak_sst = {
      pk_wk <- as.integer(names(which.max(table(tag_week))))
      mean(sst_mean[tag_week == pk_wk], na.rm = TRUE)
    },
    mean_sst = mean(sst_mean, na.rm = TRUE),
    peak_daylength = {
      pk_wk <- as.integer(names(which.max(table(tag_week))))
      mean(daylength_hrs[tag_week == pk_wk], na.rm = TRUE)
    },
    mean_daylength = mean(daylength_hrs, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n  Seasonal peaks:", nrow(seasonal_peaks), "zone-year-season observations\n")
cat("  Spring peaks:", sum(seasonal_peaks$season == "spring"), "\n")
cat("  Fall peaks:", sum(seasonal_peaks$season == "fall"), "\n")

# --- Season-specific summaries ------------------------------------------------
for (s in c("spring", "fall")) {
  sp <- seasonal_peaks %>% filter(season == s)
  if (nrow(sp) == 0) next

  cat(sprintf("\n  %s PEAK SUMMARY (across zone-years):\n", toupper(s)))
  sp_summary <- sp %>%
    group_by(tag_zone) %>%
    summarise(
      n_years = n(), mean_pk_wk = mean(peak_week), sd_pk_wk = sd(peak_week),
      mean_sst = mean(peak_sst, na.rm = TRUE), sd_sst = sd(peak_sst, na.rm = TRUE),
      mean_dl = mean(peak_daylength, na.rm = TRUE),
      .groups = "drop"
    )
  cat(sprintf("  %4s %4s %7s %6s %7s %6s %7s\n",
              "Zone", "nYr", "PkWk", "sdWk", "SST", "sdSST", "DayL"))
  for (i in seq_len(nrow(sp_summary))) {
    z <- sp_summary[i, ]
    cat(sprintf("  %4d %4d %7.1f %6.1f %6.1fC %5.1fC %6.1fh\n",
                z$tag_zone, z$n_years, z$mean_pk_wk, z$sd_pk_wk,
                z$mean_sst, z$sd_sst, z$mean_dl))
  }
}

# --- CV analysis by season ----------------------------------------------------
cat("\n  CV ANALYSIS BY SEASON:\n")
for (s in c("spring", "fall")) {
  sp <- seasonal_peaks %>% filter(season == s & !is.na(peak_sst) & !is.na(peak_daylength))
  if (nrow(sp) < 5) { cat("    ", toupper(s), ": insufficient data\n"); next }

  cv_sst <- sd(sp$peak_sst, na.rm = TRUE) / mean(sp$peak_sst, na.rm = TRUE)
  cv_jul <- sd(sp$circ_mean_julian, na.rm = TRUE) / mean(sp$circ_mean_julian, na.rm = TRUE)
  cv_day <- sd(sp$peak_daylength, na.rm = TRUE) / mean(sp$peak_daylength, na.rm = TRUE)

  cat(sprintf("  %6s (n=%d): SST CV=%.4f | Julian CV=%.4f | Daylength CV=%.4f\n",
              toupper(s), nrow(sp), cv_sst, cv_jul, cv_day))
}

# Save
write.csv(seasonal_peaks, file.path(OUTPUT_DIR, "seasonal_peaks.csv"), row.names = FALSE)
write.csv(circ_summary, file.path(OUTPUT_DIR, "circular_summary.csv"), row.names = FALSE)

# Keep yearly_peaks for backward compat (combined, for overall analysis)
yearly_peaks <- seasonal_peaks  # GAM script uses this name

cat("\n04_yearly_peaks.R complete.\n\n")
