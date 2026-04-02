# =============================================================================
# 07_pub_figures.R — Publication-ready figures v2
# Geographic labels on all figures, improved clarity
# Color system: blue = spring, orange/red = fall throughout
# =============================================================================

cat("Generating publication-ready figures (v2)...\n")

sf_use_s2(FALSE)
library(patchwork)

states <- map_data("state")
zones <- st_read(ZONE_CACHE, quiet = TRUE) %>% st_make_valid()
zone_centroids_sf <- st_centroid(zones)
centroid_coords <- st_coordinates(zone_centroids_sf)
zone_label_df <- data.frame(
  StatZone = zones$StatZone, lon = centroid_coords[,1], lat = centroid_coords[,2]
)

# --- Geographic zone labels (used across all figures) ------------------------
geo_labels <- tribble(
  ~tag_zone, ~geo_label,
  1,  "1 - FL Keys",
  2,  "2 - Dry Tortugas",
  3,  "3 - SW Florida",
  4,  "4 - Charlotte Hbr",
  5,  "5 - Tampa Bay",
  6,  "6 - Crystal River",
  7,  "7 - Big Bend FL",
  8,  "8 - FL Panhandle W",
  9,  "9 - FL Panhandle",
  10, "10 - Pensacola/AL",
  11, "11 - MS Sound",
  12, "12 - MS/LA Border",
  13, "13 - SE Louisiana",
  14, "14 - Central LA",
  15, "15 - SW Louisiana",
  16, "16 - W Louisiana",
  17, "17 - LA/TX Border",
  18, "18 - Upper TX",
  19, "19 - Central TX",
  20, "20 - S Texas",
  21, "21 - Lower TX"
)

# Short labels for scatter plot annotations
geo_short <- tribble(
  ~tag_zone, ~short_label,
  1, "Keys", 2, "Tortugas", 3, "SW FL", 4, "Charlotte",
  5, "Tampa", 6, "Crystal R", 7, "Big Bend", 8, "Panhandle W",
  9, "Panhandle", 10, "Pensacola", 11, "MS Sound", 12, "MS/LA",
  13, "SE LA", 14, "Central LA", 15, "SW LA", 16, "W LA",
  17, "LA/TX", 18, "Upper TX", 19, "Central TX", 20, "S TX", 21, "Lower TX"
)

# Merge modality into zones sf
zones_classified <- zones %>%
  left_join(circ_summary %>% select(tag_zone, modality), by = c("StatZone" = "tag_zone"))

zone_n_df <- mr_oisst %>%
  filter(!is.na(tag_zone)) %>%
  count(tag_zone, name = "n_tags") %>%
  left_join(zone_label_df, by = c("tag_zone" = "StatZone"))

# =============================================================================
# FIGURE 1: Study Area — Zone modality map with state labels
# =============================================================================
cat("  Fig 1: Study area map...\n")

fig1 <- ggplot() +
  geom_sf(data = zones_classified,
          aes(fill = modality), color = "gray30", linewidth = 0.4) +
  geom_polygon(data = states, aes(x = long, y = lat, group = group),
               fill = "gray90", color = "gray50", linewidth = 0.3) +
  geom_point(data = zone_n_df, aes(x = lon, y = lat, size = n_tags),
             shape = 21, fill = "white", color = "black", stroke = 0.6) +
  geom_text(data = zone_n_df, aes(x = lon, y = lat, label = tag_zone),
            size = 2.5, fontface = "bold") +
  # State labels
  annotate("text", x = -81.5, y = 27.5, label = "FL", size = 3.5, color = "gray40", fontface = "bold") +
  annotate("text", x = -87.5, y = 31.2, label = "AL", size = 3, color = "gray40", fontface = "bold") +
  annotate("text", x = -89.5, y = 31.2, label = "MS", size = 3, color = "gray40", fontface = "bold") +
  annotate("text", x = -91.5, y = 31.2, label = "LA", size = 3, color = "gray40", fontface = "bold") +
  annotate("text", x = -96.5, y = 29.5, label = "TX", size = 3.5, color = "gray40", fontface = "bold") +
  scale_fill_manual(
    values = c("bimodal" = "#4575B4", "unimodal_spring" = "#91BFDB",
               "unimodal_fall" = "#FC8D59"),
    labels = c("Bimodal (spring + fall)", "Unimodal (spring only)",
               "Unimodal (fall only)"),
    name = NULL, na.value = "gray85"
  ) +
  scale_size_continuous(name = "N tags", range = c(2, 12),
                        breaks = c(100, 500, 2000, 6000)) +
  coord_sf(xlim = c(-98, -80), ylim = c(24, 31.5), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", legend.box = "vertical",
        panel.grid = element_line(color = "gray95"))

ggsave(file.path(OUTPUT_DIR, "Fig1_study_area.png"), fig1,
       width = 170, height = 110, units = "mm", dpi = 300)

# =============================================================================
# FIGURE 2: Dumbbell chart with geographic labels on y-axis
# =============================================================================
cat("  Fig 2: Dumbbell chart...\n")

peak_summary <- seasonal_peaks %>%
  filter(modality == "bimodal") %>%
  group_by(tag_zone, season) %>%
  summarise(median_julian = median(circ_mean_julian, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = season, values_from = median_julian) %>%
  left_join(geo_labels, by = "tag_zone")

fig2 <- ggplot(peak_summary) +
  geom_segment(aes(x = spring, xend = fall,
                   y = reorder(geo_label, spring),
                   yend = reorder(geo_label, spring)),
               color = "gray70", linewidth = 0.8) +
  geom_point(aes(x = spring, y = reorder(geo_label, spring)),
             color = "#2171B5", size = 3.5) +
  geom_point(aes(x = fall, y = reorder(geo_label, spring)),
             color = "#E6550D", size = 3.5) +
  geom_vline(xintercept = 172, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  annotate("text", x = 172, y = 0.6, label = "Solstice", hjust = -0.1,
           size = 2.8, color = "gray50") +
  annotate("text", x = 60, y = 12.5, label = "Spring", color = "#2171B5",
           fontface = "bold", size = 3.5) +
  annotate("text", x = 290, y = 12.5, label = "Fall", color = "#E6550D",
           fontface = "bold", size = 3.5) +
  scale_x_continuous(breaks = c(1, 60, 121, 172, 244, 305),
                     labels = c("Jan", "Mar", "May", "Solstice", "Sep", "Nov"),
                     limits = c(1, 365)) +
  labs(x = "Julian day", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file.path(OUTPUT_DIR, "Fig2_dumbbell.png"), fig2,
       width = 170, height = 110, units = "mm", dpi = 300)

# =============================================================================
# FIGURE 3: Hero dual-cue panel — WIDER, geographic labels, R2 annotations
# =============================================================================
cat("  Fig 3: Hero dual-cue panel...\n")

spring_df <- seasonal_peaks %>% filter(season == "spring" & !is.na(peak_sst))
fall_df <- seasonal_peaks %>% filter(season == "fall" & !is.na(peak_daylength))

# Zone medians with geographic labels
spring_zone_med <- spring_df %>%
  group_by(tag_zone) %>%
  summarise(peak_sst = median(peak_sst, na.rm = TRUE),
            peak_week = median(peak_week), n = n(), .groups = "drop") %>%
  left_join(geo_short, by = "tag_zone")

fall_zone_med <- fall_df %>%
  group_by(tag_zone) %>%
  summarise(peak_daylength = median(peak_daylength, na.rm = TRUE),
            peak_week = median(peak_week), n = n(), .groups = "drop") %>%
  left_join(geo_short, by = "tag_zone")

# Get R2 values from GAMs
sst_r2 <- round(summary(gam(peak_week ~ s(peak_sst, k = 5), data = spring_df, method = "REML"))$r.sq * 100, 1)
dl_r2 <- round(summary(gam(peak_week ~ s(peak_daylength, k = 5), data = fall_df, method = "REML"))$r.sq * 100, 1)

p3a <- ggplot(spring_df, aes(x = peak_sst, y = peak_week)) +
  geom_point(alpha = 0.2, color = "#2171B5", size = 1.5) +
  geom_point(data = spring_zone_med, aes(x = peak_sst, y = peak_week),
             color = "#2171B5", shape = 18, size = 4) +
  geom_text(data = spring_zone_med, aes(x = peak_sst, y = peak_week, label = short_label),
            vjust = -1, size = 2, color = "gray30") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4),
              color = "#2171B5", fill = "#2171B5", alpha = 0.15, linewidth = 1) +
  annotate("text", x = Inf, y = Inf, label = "SPRING", hjust = 1.1, vjust = 1.5,
           fontface = "bold", color = "#2171B5", size = 5) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("R\u00B2 = ", sst_r2, "%"),
           hjust = 1.1, vjust = 3.5, color = "#2171B5", size = 3.5) +
  labs(x = expression("SST at peak ("*degree*"C)"), y = "Peak week of year", tag = "A") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

p3b <- ggplot(fall_df, aes(x = peak_daylength, y = peak_week)) +
  geom_point(alpha = 0.2, color = "#E6550D", size = 1.5) +
  geom_point(data = fall_zone_med, aes(x = peak_daylength, y = peak_week),
             color = "#E6550D", shape = 18, size = 4) +
  geom_text(data = fall_zone_med, aes(x = peak_daylength, y = peak_week, label = short_label),
            vjust = -1, size = 2, color = "gray30") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4),
              color = "#E6550D", fill = "#E6550D", alpha = 0.15, linewidth = 1) +
  annotate("text", x = Inf, y = Inf, label = "FALL", hjust = 1.1, vjust = 1.5,
           fontface = "bold", color = "#E6550D", size = 5) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("R\u00B2 = ", dl_r2, "%"),
           hjust = 1.1, vjust = 3.5, color = "#E6550D", size = 3.5) +
  labs(x = "Daylength at peak (hours)", y = NULL, tag = "B") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

fig3 <- p3a + p3b + plot_layout(axes = "collect_y")

ggsave(file.path(OUTPUT_DIR, "Fig3_dual_cue.png"), fig3,
       width = 200, height = 100, units = "mm", dpi = 300)

# =============================================================================
# FIGURE 4: Peak timing wave maps (spring + fall side by side)
# Does peak timing propagate geographically across zones?
# =============================================================================
cat("  Fig 4: Peak timing wave maps...\n")

tag_centroids <- mr_oisst %>%
  filter(!is.na(tag_zone)) %>%
  group_by(tag_zone) %>%
  summarise(clon = median(tag_lon), clat = median(tag_lat), .groups = "drop") %>%
  left_join(geo_short, by = "tag_zone")

# Zone-level median peak week per season
spring_timing <- seasonal_peaks %>%
  filter(season == "spring") %>%
  group_by(tag_zone) %>%
  summarise(med_peak_week = median(peak_week), n = n(), .groups = "drop") %>%
  left_join(tag_centroids, by = "tag_zone") %>%
  filter(!is.na(clon))

fall_timing <- seasonal_peaks %>%
  filter(season == "fall") %>%
  group_by(tag_zone) %>%
  summarise(med_peak_week = median(peak_week), n = n(), .groups = "drop") %>%
  left_join(tag_centroids, by = "tag_zone") %>%
  filter(!is.na(clon))

# Shared map base
map_base <- function() {
  list(
    geom_sf(data = zones, fill = "gray97", color = "gray80", linewidth = 0.3),
    geom_polygon(data = states, aes(x = long, y = lat, group = group),
                 fill = "gray90", color = "gray50", linewidth = 0.3),
    annotate("text", x = -81.5, y = 27.5, label = "FL", size = 2.5, color = "gray50", fontface = "bold"),
    annotate("text", x = -87.5, y = 31.2, label = "AL", size = 2, color = "gray50", fontface = "bold"),
    annotate("text", x = -89.5, y = 31.2, label = "MS", size = 2, color = "gray50", fontface = "bold"),
    annotate("text", x = -91.5, y = 31.2, label = "LA", size = 2, color = "gray50", fontface = "bold"),
    annotate("text", x = -96.5, y = 29.5, label = "TX", size = 2.5, color = "gray50", fontface = "bold"),
    coord_sf(xlim = c(-98, -80), ylim = c(24, 31.5), expand = FALSE),
    theme_minimal(base_size = 9),
    theme(panel.grid = element_line(color = "gray95"), legend.position = "bottom")
  )
}

p4a <- ggplot() +
  map_base() +
  geom_point(data = spring_timing,
             aes(x = clon, y = clat, color = med_peak_week, size = n), alpha = 0.9) +
  geom_text(data = spring_timing, aes(x = clon, y = clat, label = short_label),
            vjust = -1.3, size = 2, fontface = "bold") +
  scale_color_viridis_c(option = "viridis", name = "Peak\nweek",
                        breaks = c(8, 12, 16, 20, 24),
                        labels = c("Feb", "Mar", "Apr", "May", "Jun")) +
  scale_size_continuous(name = "N years", range = c(2, 8)) +
  labs(x = NULL, y = "Latitude", tag = "A") +
  annotate("text", x = -89, y = 24.5, label = "SPRING", fontface = "bold",
           color = "#2171B5", size = 4)

p4b <- ggplot() +
  map_base() +
  geom_point(data = fall_timing,
             aes(x = clon, y = clat, color = med_peak_week, size = n), alpha = 0.9) +
  geom_text(data = fall_timing, aes(x = clon, y = clat, label = short_label),
            vjust = -1.3, size = 2, fontface = "bold") +
  scale_color_viridis_c(option = "inferno", name = "Peak\nweek",
                        breaks = c(28, 32, 36, 40, 44, 48),
                        labels = c("Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  scale_size_continuous(name = "N years", range = c(2, 8)) +
  labs(x = NULL, y = NULL, tag = "B") +
  annotate("text", x = -89, y = 24.5, label = "FALL", fontface = "bold",
           color = "#E6550D", size = 4)

fig4 <- p4a + p4b
ggsave(file.path(OUTPUT_DIR, "Fig4_timing_wave.png"), fig4,
       width = 340, height = 110, units = "mm", dpi = 300)

# Keep movement arcs as supplementary
cat("  Fig S1: Movement arc map (supplementary)...\n")
recap_zoned <- mr_oisst %>%
  filter(recaptured & !is.na(tag_zone) & !is.na(recap_zone))

if (nrow(recap_zoned) > 0) {
  movement_counts <- recap_zoned %>%
    count(tag_zone, recap_zone, name = "Freq") %>%
    filter(tag_zone != recap_zone, Freq >= 5)

  movement_arcs <- movement_counts %>%
    left_join(tag_centroids %>% select(tag_zone, x = clon, y = clat), by = "tag_zone") %>%
    left_join(tag_centroids %>% select(tag_zone, xend = clon, yend = clat),
              by = c("recap_zone" = "tag_zone")) %>%
    filter(!is.na(x) & !is.na(xend)) %>%
    mutate(
      delta_lon = abs(xend - x), delta_lat = abs(yend - y),
      direction = ifelse(delta_lon > delta_lat, "Along-shelf (E-W)", "Cross-shelf (N-S)")
    )

  fig_s1 <- ggplot() +
    geom_sf(data = zones, fill = "gray97", color = "gray80", linewidth = 0.3) +
    geom_polygon(data = states, aes(x = long, y = lat, group = group),
                 fill = "gray90", color = "gray50", linewidth = 0.3) +
    geom_curve(data = movement_arcs,
               aes(x = x, y = y, xend = xend, yend = yend,
                   linewidth = Freq, color = direction),
               curvature = 0.15, alpha = 0.6,
               arrow = arrow(length = unit(1.5, "mm"), type = "closed")) +
    scale_linewidth_continuous(name = "N recaptures", range = c(0.3, 2.5)) +
    scale_color_manual(values = c("Along-shelf (E-W)" = "#2171B5",
                                   "Cross-shelf (N-S)" = "#E6550D")) +
    geom_point(data = tag_centroids, aes(x = clon, y = clat), size = 1.5, color = "black") +
    geom_text(data = tag_centroids, aes(x = clon, y = clat, label = short_label),
              size = 2, vjust = -1, fontface = "bold") +
    coord_sf(xlim = c(-98, -80), ylim = c(24, 31.5), expand = FALSE) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_line(color = "gray95"),
          legend.position = "bottom", legend.box = "vertical")

  ggsave(file.path(OUTPUT_DIR, "FigS1_movement_arcs.png"), fig_s1,
         width = 170, height = 110, units = "mm", dpi = 300)
}

# =============================================================================
# FIGURE 5: Normalized heatmap with geographic y-axis labels
# =============================================================================
cat("  Fig 5: Normalized heatmap...\n")

# Weekly catch, normalized within each zone (0-1 proportion of zone max)
weekly_catch_all <- mr_oisst %>%
  filter(!is.na(tag_zone) & tag_zone %in% valid_zones) %>%
  count(tag_zone, tag_week, name = "n") %>%
  group_by(tag_zone) %>%
  mutate(n_norm = n / max(n)) %>%
  ungroup() %>%
  left_join(geo_labels, by = "tag_zone")

# Peak annotations
peak_annot <- seasonal_peaks %>%
  group_by(tag_zone, season) %>%
  summarise(med_week = round(median(peak_week)), .groups = "drop") %>%
  left_join(geo_labels, by = "tag_zone")

fig5 <- ggplot(weekly_catch_all, aes(x = tag_week, y = reorder(geo_label, tag_zone),
                                      fill = n_norm)) +
  geom_tile(color = NA) +
  geom_point(data = peak_annot %>% filter(season == "spring"),
             aes(x = med_week, y = reorder(geo_label, tag_zone)),
             inherit.aes = FALSE, shape = 25, size = 2.5,
             fill = "#2171B5", color = "white") +
  geom_point(data = peak_annot %>% filter(season == "fall"),
             aes(x = med_week, y = reorder(geo_label, tag_zone)),
             inherit.aes = FALSE, shape = 25, size = 2.5,
             fill = "#E6550D", color = "white") +
  scale_fill_gradient(low = "white", high = "#2D004B",
                      name = "Relative\nactivity") +
  scale_x_continuous(breaks = c(1, 13, 26, 39, 52),
                     labels = c("Jan", "Apr", "Jul", "Oct", "Jan")) +
  geom_vline(xintercept = 25, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  labs(x = "Week of year", y = NULL) +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 7))

ggsave(file.path(OUTPUT_DIR, "Fig5_heatmap_annotated.png"), fig5,
       width = 170, height = 140, units = "mm", dpi = 300)

cat("  All publication figures saved.\n")
cat("07_pub_figures.R complete.\n\n")
