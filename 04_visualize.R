# =============================================================================
# 04_visualize.R — All visualizations
# Cobia M&R NOAA Zone Analysis Pipeline
# =============================================================================

cat("Generating visualizations...\n")

sf_use_s2(FALSE)

# US states for coastline context
states <- map_data("state")

# Load zone polygons for map plots
zones <- st_read(ZONE_CACHE, quiet = TRUE) %>% st_make_valid()
zone_centroids <- st_centroid(zones)
centroid_coords <- st_coordinates(zone_centroids)
zone_label_df <- data.frame(
  StatZone = zones$StatZone,
  lon = centroid_coords[, 1],
  lat = centroid_coords[, 2]
)

# =============================================================================
# PLOT 1: Zone map with tag locations
# =============================================================================
cat("  [1/5] Zone map with tag locations...\n")

p1 <- ggplot() +
  geom_sf(data = zones, aes(fill = factor(StatZone)),
          color = "black", linewidth = 0.4, alpha = 0.3) +
  geom_polygon(data = states, aes(x = long, y = lat, group = group),
               fill = "gray85", color = "gray40", linewidth = 0.3) +
  geom_point(data = mr_data %>% filter(!is.na(tag_zone)),
             aes(x = tag_lon, y = tag_lat, color = dataset),
             size = 0.8, alpha = 0.5) +
  geom_text(data = zone_label_df, aes(x = lon, y = lat, label = StatZone),
            fontface = "bold", size = 3, color = "black") +
  scale_fill_viridis_d(option = "turbo", guide = "none") +
  scale_color_manual(values = c("SEFSC_CTC" = "#E66100", "USM_Cooperative" = "#5D3A9B"),
                     name = "Dataset") +
  coord_sf(xlim = c(-98, -80), ylim = c(24, 31.5), expand = FALSE) +
  labs(title = "Cobia Mark-and-Recapture Tag Locations by NOAA Statistical Zone",
       subtitle = paste0("N = ", sum(!is.na(mr_data$tag_zone)), " tagged fish in zones; ",
                         sum(is.na(mr_data$tag_zone)), " outside zones"),
       x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUTPUT_DIR, "01_zone_map_tags.png"), p1,
       width = 14, height = 7, dpi = 300)

# =============================================================================
# PLOT 2: Heatmap — week of year x zone (the migration diagnostic)
# =============================================================================
cat("  [2/5] Week x Zone heatmap...\n")

# Filter out low-sample zones
heatmap_data <- weekly_catch %>%
  filter(!tag_zone %in% low_n_zones)

p2 <- ggplot(heatmap_data, aes(x = tag_week, y = factor(tag_zone), fill = n)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis_c(option = "magma", trans = "sqrt",
                       name = "Tag count", na.value = "gray95") +
  scale_x_continuous(breaks = seq(1, 52, by = 4),
                     labels = function(x) {
                       # Convert week number to approximate month
                       month.abb[pmin(ceiling(x / 4.33), 12)]
                     }) +
  labs(title = "Cobia Tagging Activity by Week and NOAA Statistical Zone",
       subtitle = "Diagonal pattern = northward migration; vertical bands = synchronous activity",
       x = "Week of Year (approximate month)", y = "Statistical Zone (S to N/W)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid = element_blank())

ggsave(file.path(OUTPUT_DIR, "02_heatmap_week_zone.png"), p2,
       width = 14, height = 8, dpi = 300)

# =============================================================================
# PLOT 3: Ridge plot — Julian day density per zone
# =============================================================================
cat("  [3/5] Ridge plot (Julian day by zone)...\n")

ridge_data <- zoned_data %>%
  filter(!tag_zone %in% low_n_zones)

p3 <- ggplot(ridge_data, aes(x = tag_julian, y = factor(tag_zone),
                              fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 2.5, rel_min_height = 0.01,
                                quantile_lines = TRUE, quantiles = 2) +
  scale_fill_viridis_c(option = "plasma", name = "Julian Day") +
  scale_x_continuous(breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
                     labels = month.abb) +
  labs(title = "Seasonal Distribution of Cobia Tagging by NOAA Statistical Zone",
       subtitle = "Density ridges with median line; zones ordered S to N/W",
       x = "Julian Day (Month)", y = "Statistical Zone") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "none")

ggsave(file.path(OUTPUT_DIR, "03_ridgeplot_seasonal.png"), p3,
       width = 12, height = 10, dpi = 300)

# =============================================================================
# PLOT 4: Faceted weekly time series per zone
# =============================================================================
cat("  [4/5] Faceted weekly time series...\n")

facet_data <- weekly_catch %>%
  filter(!tag_zone %in% low_n_zones & n > 0)

p4 <- ggplot(facet_data, aes(x = tag_week, y = n)) +
  geom_col(fill = "#2171B5", alpha = 0.7) +
  facet_wrap(~tag_zone, ncol = 4, scales = "free_y",
             labeller = labeller(tag_zone = function(x) paste("Zone", x))) +
  scale_x_continuous(breaks = seq(1, 52, by = 13),
                     labels = c("Jan", "Apr", "Jul", "Oct")) +
  labs(title = "Weekly Tagging Effort by Statistical Zone (Pooled Across Years)",
       x = "Week of Year", y = "Number of Tags") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())

ggsave(file.path(OUTPUT_DIR, "04_faceted_weekly.png"), p4,
       width = 14, height = 12, dpi = 300)

# =============================================================================
# PLOT 5: Movement matrix heatmap (tag zone vs recap zone)
# =============================================================================
cat("  [5/5] Movement matrix heatmap...\n")

if (nrow(movement_long) > 0 && sum(movement_long$Freq) > 0) {
  p5 <- ggplot(movement_long, aes(x = factor(Recap_Zone), y = factor(Tag_Zone),
                                   fill = Freq)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(Freq > 0, Freq, "")),
              color = "white", size = 3, fontface = "bold") +
    scale_fill_viridis_c(option = "inferno", name = "Count",
                         trans = "sqrt", na.value = "gray95") +
    labs(title = "Cobia Tag-Recapture Movement Matrix",
         subtitle = "Diagonal = same-zone recapture; off-diagonal = zone movement",
         x = "Recapture Zone", y = "Tag Zone") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          panel.grid = element_blank(),
          aspect.ratio = 1)

  ggsave(file.path(OUTPUT_DIR, "05_movement_matrix.png"), p5,
         width = 10, height = 10, dpi = 300)
} else {
  cat("    (skipped — no recaptured fish with both zones)\n")
}

# =============================================================================
# PLOT 6 (bonus): Peak timing map — zone centroids colored by peak Julian day
# =============================================================================
cat("  [bonus] Peak timing map...\n")

peak_map_df <- zone_label_df %>%
  left_join(peak_by_zone, by = c("StatZone" = "tag_zone")) %>%
  filter(!is.na(median_julian))

if (nrow(peak_map_df) > 0) {
  p6 <- ggplot() +
    geom_sf(data = zones, fill = "gray95", color = "gray60", linewidth = 0.3) +
    geom_polygon(data = states, aes(x = long, y = lat, group = group),
                 fill = "gray85", color = "gray40", linewidth = 0.3) +
    geom_point(data = peak_map_df,
               aes(x = lon, y = lat, color = median_julian, size = n_tags),
               alpha = 0.9) +
    geom_text(data = peak_map_df, aes(x = lon, y = lat, label = StatZone),
              vjust = -1.5, fontface = "bold", size = 3) +
    scale_color_viridis_c(option = "plasma", name = "Median\nJulian Day",
                          breaks = c(60, 120, 180, 240, 300),
                          labels = c("Mar", "May", "Jul", "Sep", "Nov")) +
    scale_size_continuous(name = "N Tags", range = c(2, 12)) +
    coord_sf(xlim = c(-98, -80), ylim = c(24, 31.5), expand = FALSE) +
    labs(title = "Peak Tagging Season by Statistical Zone",
         subtitle = "Color = median Julian day; size = sample size") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "right")

  ggsave(file.path(OUTPUT_DIR, "06_peak_timing_map.png"), p6,
         width = 14, height = 7, dpi = 300)
}

cat("\nAll plots saved to", OUTPUT_DIR, "\n")
cat("04_visualize.R complete.\n\n")
