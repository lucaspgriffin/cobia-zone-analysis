# =============================================================================
# 05_sst_analysis.R — Extract SST per zone, correlate with peak timing
# Cobia M&R NOAA Zone Analysis Pipeline
#
# SST source: NOAA OISST v2.1 (daily, 0.25-deg) via ERDDAP
# Uses 5-year climatology (2019-2023) at zone centroids
#
# Questions:
#   1) Is peak tagging timing better predicted by SST or Julian day?
#   2) Does the order of peaks across zones follow a temperature wave
#      (migration tracking isotherms) or a fixed calendar signal?
# =============================================================================

cat("Running SST analysis...\n")

sf_use_s2(FALSE)

# =============================================================================
# STEP 1: Compute zone centroids for SST extraction
# =============================================================================

zones <- st_read(ZONE_CACHE, quiet = TRUE) %>% st_make_valid()
zone_centroids <- st_centroid(zones)
centroid_coords <- st_coordinates(zone_centroids)

zone_coords <- data.frame(
  StatZone = zones$StatZone,
  lon = round(centroid_coords[, 1], 2),
  lat = round(centroid_coords[, 2], 2)
)

cat("  Zone centroids computed for SST extraction\n")

# =============================================================================
# STEP 2: Extract SST climatology per zone from NOAA OISST via ERDDAP
# =============================================================================

SST_CACHE <- file.path(SHAPE_DIR, "zone_sst_climatology.csv")

# ERDDAP OISST v2.1 requires:
#   - zlev dimension = 0.0 (surface)
#   - longitude in 0-360 range (not -180 to 180)
#   - time with T12:00:00Z

if (!file.exists(SST_CACHE)) {
  cat("  Downloading SST from NOAA OISST v2.1 via ERDDAP...\n")
  cat("  (5-year daily climatology 2019-2023 at zone centroids)\n")

  erddap_base <- "https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOisst21Agg.csv"

  all_sst <- list()

  for (i in seq_len(nrow(zone_coords))) {
    z <- zone_coords[i, ]
    lon360 <- z$lon + 360  # Convert -180:180 to 0:360

    cat("    Zone", z$StatZone, "(", z$lat, ",", z$lon, " -> lon360=", lon360, ")...")

    # Query 5 years of daily SST at this point
    url <- paste0(erddap_base,
                  "?sst[(2019-01-01T12:00:00Z):1:(2023-12-31T12:00:00Z)]",
                  "[(0.0):1:(0.0)]",
                  "[(", z$lat, "):1:(", z$lat, ")]",
                  "[(", lon360, "):1:(", lon360, ")]")

    result <- tryCatch({
      tmp <- tempfile(fileext = ".csv")
      download.file(url, tmp, quiet = TRUE, method = "libcurl")
      d <- read.csv(tmp, skip = 1, stringsAsFactors = FALSE)
      names(d) <- c("time", "zlev", "latitude", "longitude", "sst")
      d$sst <- as.numeric(d$sst)
      d$time <- as.Date(substr(d$time, 1, 10))
      d$StatZone <- z$StatZone
      d$julian <- yday(d$time)
      d$week <- isoweek(d$time)
      d$month <- month(d$time)
      d$year <- year(d$time)
      d <- d %>% select(time, StatZone, sst, julian, week, month, year)
      unlink(tmp)
      cat(" OK (", nrow(d), "days)\n")
      d
    }, error = function(e) {
      cat(" FAILED:", e$message, "\n")
      NULL
    })

    if (!is.null(result)) all_sst[[length(all_sst) + 1]] <- result
    Sys.sleep(0.5)  # Be polite to ERDDAP
  }

  sst_raw <- bind_rows(all_sst)
  cat("  Total SST records:", nrow(sst_raw), "\n")

  # Compute weekly SST climatology per zone (mean across 5 years)
  sst_weekly <- sst_raw %>%
    group_by(StatZone, week) %>%
    summarise(
      sst_mean = mean(sst, na.rm = TRUE),
      sst_sd   = sd(sst, na.rm = TRUE),
      n_obs    = n(),
      .groups = "drop"
    )

  # Monthly climatology
  sst_monthly <- sst_raw %>%
    group_by(StatZone, month) %>%
    summarise(
      sst_mean = mean(sst, na.rm = TRUE),
      sst_sd   = sd(sst, na.rm = TRUE),
      .groups = "drop"
    )

  write.csv(sst_weekly, SST_CACHE, row.names = FALSE)
  write.csv(sst_monthly, file.path(SHAPE_DIR, "zone_sst_monthly.csv"), row.names = FALSE)
  cat("  SST climatology cached to", SST_CACHE, "\n")

} else {
  cat("  Using cached SST climatology:", SST_CACHE, "\n")
  sst_weekly <- read.csv(SST_CACHE)
  sst_monthly_file <- file.path(SHAPE_DIR, "zone_sst_monthly.csv")
  if (file.exists(sst_monthly_file)) {
    sst_monthly <- read.csv(sst_monthly_file)
  }
}

# =============================================================================
# STEP 3: Match SST to tagging peaks
# =============================================================================

# Get SST at peak week for each zone
peak_sst <- zone_summary %>%
  select(tag_zone, peak_week, median_julian, n_tags) %>%
  left_join(
    sst_weekly %>%
      rename(tag_zone = StatZone, peak_sst = sst_mean) %>%
      select(tag_zone, week, peak_sst),
    by = c("tag_zone", "peak_week" = "week")
  )

# Also get SST at the tag week for each individual fish
tag_sst <- mr_data %>%
  filter(!is.na(tag_zone) & !is.na(tag_week)) %>%
  left_join(
    sst_weekly %>%
      rename(tag_zone = StatZone, tag_sst = sst_mean) %>%
      select(tag_zone, week, tag_sst),
    by = c("tag_zone", "tag_week" = "week")
  )

cat("\n  PEAK SST BY ZONE:\n")
cat("  ", paste(rep("-", 60), collapse = ""), "\n")
cat(sprintf("  %4s %6s %7s %8s %8s\n",
            "Zone", "nTags", "PkWeek", "MedJul", "PkSST_C"))
cat("  ", paste(rep("-", 60), collapse = ""), "\n")
for (i in seq_len(nrow(peak_sst))) {
  p <- peak_sst[i, ]
  cat(sprintf("  %4d %6d %7d %8.0f %7.1f\n",
              p$tag_zone, p$n_tags, p$peak_week, p$median_julian,
              ifelse(is.na(p$peak_sst), NA, p$peak_sst)))
}

# =============================================================================
# STEP 4: Correlation analysis — what drives peak timing?
# =============================================================================

cat("\n  CORRELATION ANALYSIS:\n")

valid_peaks <- peak_sst %>% filter(!is.na(peak_sst) & n_tags >= 10)

if (nrow(valid_peaks) >= 5) {
  # A) Cross-zone: are peaks at the same SST (thermal tracking)?
  cor_sst_week <- cor.test(valid_peaks$peak_sst, valid_peaks$peak_week,
                           method = "spearman")
  cor_jul_week <- cor.test(valid_peaks$median_julian, valid_peaks$peak_week,
                           method = "spearman")

  cat("  Cross-zone (n =", nrow(valid_peaks), "zones with >= 10 tags):\n")
  cat("    Peak SST vs Peak Week:  rho =", round(cor_sst_week$estimate, 3),
      " p =", round(cor_sst_week$p.value, 4), "\n")
  cat("    Median Julian vs Peak Week:  rho =", round(cor_jul_week$estimate, 3),
      " p =", round(cor_jul_week$p.value, 4), "\n")

  # B) CV comparison: if fish track isotherms, peak SST should be consistent
  cv_sst <- sd(valid_peaks$peak_sst, na.rm = TRUE) / mean(valid_peaks$peak_sst, na.rm = TRUE)
  cv_jul <- sd(valid_peaks$median_julian, na.rm = TRUE) / mean(valid_peaks$median_julian, na.rm = TRUE)

  cat("\n  Coefficient of Variation (lower = more consistent = likely driver):\n")
  cat("    Peak SST across zones: CV =", round(cv_sst, 3), "\n")
  cat("    Peak Julian day across zones: CV =", round(cv_jul, 3), "\n")
  if (cv_sst < cv_jul) {
    cat("    -> SST is more consistent across zones (supports thermal tracking)\n")
  } else {
    cat("    -> Julian day is more consistent (supports calendar/photoperiod cue)\n")
  }
}

# C) Individual-level SST at tagging
if (nrow(tag_sst) > 0 & sum(!is.na(tag_sst$tag_sst)) > 100) {
  cat("\n  Individual-level SST at tagging:\n")
  cat("    Mean SST:", round(mean(tag_sst$tag_sst, na.rm = TRUE), 1), "C\n")
  cat("    SD:", round(sd(tag_sst$tag_sst, na.rm = TRUE), 1), "C\n")
  cat("    Range:", round(min(tag_sst$tag_sst, na.rm = TRUE), 1), "-",
      round(max(tag_sst$tag_sst, na.rm = TRUE), 1), "C\n")
  q10 <- quantile(tag_sst$tag_sst, 0.10, na.rm = TRUE)
  q90 <- quantile(tag_sst$tag_sst, 0.90, na.rm = TRUE)
  cat("    80% of tagging occurs between",
      round(q10, 1), "-", round(q90, 1), "C\n")
}

# =============================================================================
# STEP 5: SST Visualizations
# =============================================================================

cat("\n  Generating SST plots...\n")
states <- map_data("state")

# --- Plot 7: Peak SST vs Peak Week scatter -----------------------------------
cat("  [7] Peak SST vs peak week scatter...\n")
if (nrow(valid_peaks) >= 3) {
  p7 <- ggplot(valid_peaks, aes(x = peak_sst, y = peak_week)) +
    geom_point(aes(color = median_julian, size = n_tags), alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "gray40", linetype = "dashed") +
    geom_text(aes(label = tag_zone), vjust = -1, size = 3, fontface = "bold") +
    scale_color_viridis_c(option = "plasma", name = "Median\nJulian Day") +
    scale_size_continuous(name = "N Tags", range = c(3, 12)) +
    labs(title = "Peak SST vs Peak Week by Zone",
         subtitle = "Horizontal clustering = thermal tracking; vertical = calendar-driven",
         x = "SST at Peak Week (C)", y = "Peak Week of Year") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(OUTPUT_DIR, "07_peak_sst_vs_week.png"), p7,
         width = 10, height = 7, dpi = 300)
}

# --- Plot 8: SST heatmap by zone and week ------------------------------------
cat("  [8] SST climatology heatmap...\n")
if (nrow(sst_weekly) > 0) {
  p8 <- ggplot(sst_weekly, aes(x = week, y = factor(StatZone), fill = sst_mean)) +
    geom_tile(color = NA) +
    scale_fill_viridis_c(option = "inferno", name = "SST (C)") +
    scale_x_continuous(breaks = seq(1, 52, by = 4),
                       labels = function(x) month.abb[pmin(ceiling(x / 4.33), 12)]) +
    labs(title = "Weekly SST Climatology by NOAA Statistical Zone (2019-2023)",
         subtitle = "Compare with tagging heatmap — do tags follow an isotherm band?",
         x = "Week of Year", y = "Statistical Zone") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), panel.grid = element_blank())
  ggsave(file.path(OUTPUT_DIR, "08_sst_heatmap_zone_week.png"), p8,
         width = 14, height = 8, dpi = 300)
}

# --- Plot 9: SST ridge plot at tagging by zone --------------------------------
cat("  [9] SST at tagging ridge plot...\n")
if (nrow(tag_sst) > 0 & sum(!is.na(tag_sst$tag_sst)) > 100) {
  sst_ridge_data <- tag_sst %>%
    filter(!is.na(tag_sst) & !tag_zone %in% low_n_zones)

  p9 <- ggplot(sst_ridge_data, aes(x = tag_sst, y = factor(tag_zone),
                                    fill = after_stat(x))) +
    geom_density_ridges_gradient(scale = 2.5, rel_min_height = 0.01,
                                  quantile_lines = TRUE, quantiles = 2) +
    scale_fill_viridis_c(option = "inferno", name = "SST (C)") +
    labs(title = "SST at Tagging by NOAA Statistical Zone",
         subtitle = "Overlapping distributions = thermal tracking; separated = calendar-driven",
         x = "Sea Surface Temperature (C)", y = "Statistical Zone") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")
  ggsave(file.path(OUTPUT_DIR, "09_sst_ridgeplot_by_zone.png"), p9,
         width = 10, height = 10, dpi = 300)
}

# --- Plot 10: Peak SST map ---------------------------------------------------
cat("  [10] Peak SST map...\n")
if (nrow(valid_peaks) >= 3) {
  peak_sst_map <- zone_coords %>%
    left_join(valid_peaks, by = c("StatZone" = "tag_zone")) %>%
    filter(!is.na(peak_sst))

  if (nrow(peak_sst_map) > 0) {
    p10 <- ggplot() +
      geom_sf(data = zones, fill = "gray95", color = "gray60", linewidth = 0.3) +
      geom_polygon(data = states, aes(x = long, y = lat, group = group),
                   fill = "gray85", color = "gray40", linewidth = 0.3) +
      geom_point(data = peak_sst_map,
                 aes(x = lon, y = lat, color = peak_sst, size = n_tags),
                 alpha = 0.9) +
      geom_text(data = peak_sst_map,
                aes(x = lon, y = lat, label = StatZone),
                vjust = -1.5, fontface = "bold", size = 3) +
      scale_color_viridis_c(option = "inferno", name = "Peak SST\n(C)") +
      scale_size_continuous(name = "N Tags", range = c(2, 12)) +
      coord_sf(xlim = c(-98, -80), ylim = c(24, 31.5), expand = FALSE) +
      labs(title = "SST at Peak Tagging Week by Statistical Zone",
           subtitle = "Similar colors = fish active at similar temps (thermal tracking)") +
      theme_minimal(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))
    ggsave(file.path(OUTPUT_DIR, "10_peak_sst_map.png"), p10,
           width = 14, height = 7, dpi = 300)
  }
}

# Save summary
write.csv(peak_sst, file.path(OUTPUT_DIR, "peak_sst_by_zone.csv"), row.names = FALSE)

cat("\n05_sst_analysis.R complete.\n\n")
