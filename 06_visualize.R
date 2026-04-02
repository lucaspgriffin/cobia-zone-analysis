# =============================================================================
# 06_visualize.R — Diagnostic and summary visualizations (v2)
# =============================================================================

cat("Generating visualizations...\n")

sf_use_s2(FALSE)
states <- map_data("state")

zones <- st_read(ZONE_CACHE, quiet = TRUE) %>% st_make_valid()
zone_centroids <- st_centroid(zones)
centroid_coords <- st_coordinates(zone_centroids)
zone_label_df <- data.frame(
  StatZone = zones$StatZone,
  lon = centroid_coords[, 1],
  lat = centroid_coords[, 2]
)

# =============================================================================
# 1. Zone map with tag locations (OISST era only)
# =============================================================================
cat("  [1] Zone map with tags...\n")

p1 <- ggplot() +
  geom_sf(data = zones, aes(fill = factor(StatZone)),
          color = "black", linewidth = 0.4, alpha = 0.3) +
  geom_polygon(data = states, aes(x = long, y = lat, group = group),
               fill = "gray85", color = "gray40", linewidth = 0.3) +
  geom_point(data = mr_oisst %>% filter(!is.na(tag_zone)),
             aes(x = tag_lon, y = tag_lat, color = dataset),
             size = 0.6, alpha = 0.4) +
  geom_text(data = zone_label_df, aes(x = lon, y = lat, label = StatZone),
            fontface = "bold", size = 3) +
  scale_fill_viridis_d(option = "turbo", guide = "none") +
  scale_color_manual(values = c("SEFSC_CTC" = "#E66100", "USM_Cooperative" = "#5D3A9B")) +
  coord_sf(xlim = c(-98, -80), ylim = c(24, 31.5), expand = FALSE) +
  labs(title = "Cobia M&R Tag Locations (1982-2025) by NOAA Statistical Zone",
       x = "Longitude", y = "Latitude", color = "Dataset") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR, "01_zone_map.png"), p1, width = 14, height = 7, dpi = 300)

# =============================================================================
# 2. Yearly peak week by zone (the key diagnostic — NOT pooled)
# =============================================================================
cat("  [2] Yearly peak week by zone...\n")

p2 <- ggplot(yearly_peaks, aes(x = tag_year, y = peak_week, color = factor(tag_zone))) +
  geom_point(aes(size = n_tags), alpha = 0.6) +
  geom_line(alpha = 0.3) +
  facet_wrap(~tag_zone, ncol = 4, scales = "free_y",
             labeller = labeller(tag_zone = function(x) paste("Zone", x))) +
  scale_color_viridis_d(option = "turbo", guide = "none") +
  scale_size_continuous(name = "N Tags", range = c(0.5, 3)) +
  labs(title = "Peak Tagging Week by Year and Zone",
       subtitle = "Each point = peak week for that zone-year (not pooled across years)",
       x = "Year", y = "Peak Week") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")
ggsave(file.path(OUTPUT_DIR, "02_yearly_peaks_by_zone.png"), p2,
       width = 16, height = 12, dpi = 300)

# =============================================================================
# 3. SST at peak vs peak week (zone-year level scatter)
# =============================================================================
cat("  [3] SST at peak vs peak week...\n")

p3 <- ggplot(yearly_peaks %>% filter(!is.na(peak_sst)),
             aes(x = peak_sst, y = peak_week)) +
  geom_point(aes(color = factor(tag_zone), size = n_tags), alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 6),
              color = "black", linewidth = 1, se = TRUE) +
  scale_color_viridis_d(option = "turbo", name = "Zone") +
  scale_size_continuous(name = "N Tags", range = c(1, 5)) +
  labs(title = "Peak SST vs Peak Week (per zone-year, not pooled)",
       subtitle = "Tight horizontal band = thermal tracking; scattered = no thermal cue",
       x = "SST at Peak Week (C)", y = "Peak Week of Year") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR, "03_sst_vs_peak_week.png"), p3,
       width = 10, height = 7, dpi = 300)

# =============================================================================
# 4. SST ridgeplot at tagging (individual-level)
# =============================================================================
cat("  [4] SST at tagging ridge plot...\n")

sst_ridge_data <- mr_oisst %>%
  filter(!is.na(sst_mean) & !is.na(tag_zone) & tag_zone %in% valid_zones)

if (nrow(sst_ridge_data) > 100) {
  p4 <- ggplot(sst_ridge_data, aes(x = sst_mean, y = factor(tag_zone),
                                    fill = after_stat(x))) +
    geom_density_ridges_gradient(scale = 2.5, rel_min_height = 0.01,
                                  quantile_lines = TRUE, quantiles = 2) +
    scale_fill_viridis_c(option = "inferno", name = "SST (C)") +
    labs(title = "SST at Tagging by Zone (year-matched, individual level)",
         subtitle = "Overlapping = thermal tracking; separated = calendar-driven",
         x = "Sea Surface Temperature (C)", y = "Statistical Zone") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")
  ggsave(file.path(OUTPUT_DIR, "04_sst_ridgeplot.png"), p4,
         width = 10, height = 10, dpi = 300)
}

# =============================================================================
# 5. Julian day ridgeplot (seasonal timing by zone)
# =============================================================================
cat("  [5] Julian day ridge plot...\n")

p5 <- ggplot(mr_oisst %>% filter(tag_zone %in% valid_zones),
             aes(x = tag_julian, y = factor(tag_zone), fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 2.5, rel_min_height = 0.01,
                                quantile_lines = TRUE, quantiles = 2) +
  scale_fill_viridis_c(option = "plasma", name = "Julian Day") +
  scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                     labels = month.abb) +
  labs(title = "Seasonal Distribution of Cobia Tagging by Zone (1982-2025)",
       x = "Julian Day", y = "Statistical Zone") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")
ggsave(file.path(OUTPUT_DIR, "05_julian_ridgeplot.png"), p5,
       width = 12, height = 10, dpi = 300)

# =============================================================================
# 6. Daylength at tagging ridge plot
# =============================================================================
cat("  [6] Daylength ridge plot...\n")

p6 <- ggplot(mr_oisst %>% filter(tag_zone %in% valid_zones & !is.na(daylength_hrs)),
             aes(x = daylength_hrs, y = factor(tag_zone), fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 2.5, rel_min_height = 0.01,
                                quantile_lines = TRUE, quantiles = 2) +
  scale_fill_viridis_c(option = "mako", name = "Day Length (h)") +
  labs(title = "Daylength at Tagging by Zone",
       subtitle = "Overlapping = photoperiod cue; separated = other driver",
       x = "Day Length (hours)", y = "Statistical Zone") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "none")
ggsave(file.path(OUTPUT_DIR, "06_daylength_ridgeplot.png"), p6,
       width = 10, height = 10, dpi = 300)

# =============================================================================
# 7. Movement matrix
# =============================================================================
cat("  [7] Movement matrix...\n")

recap_zoned <- mr_oisst %>%
  filter(recaptured & !is.na(tag_zone) & !is.na(recap_zone))

if (nrow(recap_zoned) > 0) {
  movement_long <- as.data.frame(table(
    Tag_Zone = recap_zoned$tag_zone,
    Recap_Zone = recap_zoned$recap_zone
  )) %>%
    mutate(Tag_Zone = as.integer(as.character(Tag_Zone)),
           Recap_Zone = as.integer(as.character(Recap_Zone)))

  p7 <- ggplot(movement_long, aes(x = factor(Recap_Zone), y = factor(Tag_Zone), fill = Freq)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(Freq > 0, Freq, "")),
              color = "white", size = 2.5, fontface = "bold") +
    scale_fill_viridis_c(option = "inferno", name = "Count", trans = "sqrt") +
    labs(title = "Tag-Recapture Movement Matrix",
         x = "Recapture Zone", y = "Tag Zone") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), panel.grid = element_blank(),
          aspect.ratio = 1)
  ggsave(file.path(OUTPUT_DIR, "07_movement_matrix.png"), p7,
         width = 10, height = 10, dpi = 300)
}

# =============================================================================
# 8. Bimodality diagnostic — zones flagged by dip test
# =============================================================================
cat("  [8] Bimodality diagnostic...\n")

if (length(bimodal_zones) > 0) {
  bimod_data <- mr_oisst %>%
    filter(tag_zone %in% bimodal_zones)

  p8 <- ggplot(bimod_data, aes(x = tag_julian)) +
    geom_histogram(aes(fill = after_stat(x)), bins = 52, color = "white", linewidth = 0.2) +
    facet_wrap(~tag_zone, ncol = 3, scales = "free_y",
               labeller = labeller(tag_zone = function(x) paste("Zone", x, "(bimodal)"))) +
    scale_fill_viridis_c(option = "plasma") +
    scale_x_continuous(breaks = c(1, 91, 182, 274),
                       labels = c("Jan", "Apr", "Jul", "Oct")) +
    labs(title = "Bimodal Zones: Julian Day Distributions (Hartigan dip test p < 0.05)",
         subtitle = "Spring vs fall peaks — may represent different behavioral modes",
         x = "Julian Day", y = "Count") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "none")
  ggsave(file.path(OUTPUT_DIR, "08_bimodal_zones.png"), p8,
         width = 12, height = 8, dpi = 300)
}

cat("\nAll plots saved to", OUTPUT_DIR, "\n")
cat("06_visualize.R complete.\n\n")
