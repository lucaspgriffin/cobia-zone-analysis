# =============================================================================
# 07_pub_figures.R — Publication-ready figures for the results document
# 5 figures: study area, dumbbell, hero dual-cue, movement arcs, heatmap
# Color system: blue = spring, orange/red = fall throughout
# =============================================================================

cat("Generating publication-ready figures...\n")

sf_use_s2(FALSE)
library(patchwork)

states <- map_data("state")
zones <- st_read(ZONE_CACHE, quiet = TRUE) %>% st_make_valid()
zone_centroids_sf <- st_centroid(zones)
centroid_coords <- st_coordinates(zone_centroids_sf)
zone_label_df <- data.frame(
  StatZone = zones$StatZone, lon = centroid_coords[,1], lat = centroid_coords[,2]
)

# Merge modality into zones sf
zones_classified <- zones %>%
  left_join(circ_summary %>% select(tag_zone, modality), by = c("StatZone" = "tag_zone"))

# Zone sample sizes
zone_n_df <- mr_oisst %>%
  filter(!is.na(tag_zone)) %>%
  count(tag_zone, name = "n_tags") %>%
  left_join(zone_label_df, by = c("tag_zone" = "StatZone"))

# =============================================================================
# FIGURE 1: Study Area — Zone modality map
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
  scale_fill_manual(
    values = c("bimodal" = "#4575B4",
               "unimodal_spring" = "#91BFDB",
               "unimodal_fall" = "#FC8D59"),
    labels = c("Bimodal (spring + fall)",
               "Unimodal (spring only)",
               "Unimodal (fall only)"),
    name = NULL, na.value = "gray85"
  ) +
  scale_size_continuous(name = "N tags", range = c(2, 12),
                        breaks = c(100, 500, 2000, 6000)) +
  coord_sf(xlim = c(-98, -80), ylim = c(24, 31.5), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        panel.grid = element_line(color = "gray95"),
        plot.margin = margin(5, 5, 5, 5))

ggsave(file.path(OUTPUT_DIR, "Fig1_study_area.png"), fig1,
       width = 170, height = 100, units = "mm", dpi = 300)

# =============================================================================
# FIGURE 2: Dumbbell chart — spring vs fall peak timing per zone
# =============================================================================
cat("  Fig 2: Dumbbell chart...\n")

peak_summary <- seasonal_peaks %>%
  filter(modality == "bimodal") %>%
  group_by(tag_zone, season) %>%
  summarise(median_julian = median(circ_mean_julian, na.rm = TRUE),
            median_week = median(peak_week, na.rm = TRUE),
            n_years = n(), .groups = "drop") %>%
  pivot_wider(names_from = season, values_from = c(median_julian, median_week, n_years))

fig2 <- ggplot(peak_summary) +
  geom_segment(aes(x = median_julian_spring, xend = median_julian_fall,
                   y = reorder(factor(tag_zone), median_julian_spring),
                   yend = reorder(factor(tag_zone), median_julian_spring)),
               color = "gray70", linewidth = 0.8) +
  geom_point(aes(x = median_julian_spring,
                 y = reorder(factor(tag_zone), median_julian_spring)),
             color = "#2171B5", size = 3.5) +
  geom_point(aes(x = median_julian_fall,
                 y = reorder(factor(tag_zone), median_julian_spring)),
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
  labs(x = "Julian day", y = "Statistical zone") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file.path(OUTPUT_DIR, "Fig2_dumbbell.png"), fig2,
       width = 170, height = 100, units = "mm", dpi = 300)

# =============================================================================
# FIGURE 3: Hero dual-cue panel (SST spring | Daylength fall)
# =============================================================================
cat("  Fig 3: Hero dual-cue panel...\n")

spring_df <- seasonal_peaks %>% filter(season == "spring" & !is.na(peak_sst))
fall_df <- seasonal_peaks %>% filter(season == "fall" & !is.na(peak_daylength))

# Zone-level medians for emphasis
spring_zone_med <- spring_df %>%
  group_by(tag_zone) %>%
  summarise(peak_sst = median(peak_sst, na.rm = TRUE),
            peak_week = median(peak_week), n = n(), .groups = "drop")
fall_zone_med <- fall_df %>%
  group_by(tag_zone) %>%
  summarise(peak_daylength = median(peak_daylength, na.rm = TRUE),
            peak_week = median(peak_week), n = n(), .groups = "drop")

p3a <- ggplot(spring_df, aes(x = peak_sst, y = peak_week)) +
  geom_point(alpha = 0.25, color = "#2171B5", size = 1.5) +
  geom_point(data = spring_zone_med, aes(size = n),
             color = "#2171B5", shape = 18, size = 3.5) +
  geom_text(data = spring_zone_med, aes(label = tag_zone),
            vjust = -1, size = 2.3, color = "gray30") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4),
              color = "#2171B5", fill = "#2171B5", alpha = 0.15, linewidth = 1) +
  annotate("text", x = Inf, y = Inf, label = "SPRING", hjust = 1.1, vjust = 1.5,
           fontface = "bold", color = "#2171B5", size = 4.5) +
  labs(x = expression("SST at peak ("*degree*"C)"), y = "Peak week of year", tag = "A") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

p3b <- ggplot(fall_df, aes(x = peak_daylength, y = peak_week)) +
  geom_point(alpha = 0.25, color = "#E6550D", size = 1.5) +
  geom_point(data = fall_zone_med, aes(size = n),
             color = "#E6550D", shape = 18, size = 3.5) +
  geom_text(data = fall_zone_med, aes(label = tag_zone),
            vjust = -1, size = 2.3, color = "gray30") +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 4),
              color = "#E6550D", fill = "#E6550D", alpha = 0.15, linewidth = 1) +
  annotate("text", x = Inf, y = Inf, label = "FALL", hjust = 1.1, vjust = 1.5,
           fontface = "bold", color = "#E6550D", size = 4.5) +
  labs(x = "Daylength at peak (hours)", y = NULL, tag = "B") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none", panel.grid.minor = element_blank())

fig3 <- p3a + p3b + plot_layout(axes = "collect_y")

ggsave(file.path(OUTPUT_DIR, "Fig3_dual_cue.png"), fig3,
       width = 170, height = 85, units = "mm", dpi = 300)

# =============================================================================
# FIGURE 4: Movement arc map
# =============================================================================
cat("  Fig 4: Movement arc map...\n")

# Compute median tag locations per zone for arc endpoints
tag_centroids <- mr_oisst %>%
  filter(!is.na(tag_zone)) %>%
  group_by(tag_zone) %>%
  summarise(clon = median(tag_lon), clat = median(tag_lat), .groups = "drop")

recap_zoned <- mr_oisst %>%
  filter(recaptured & !is.na(tag_zone) & !is.na(recap_zone))

if (nrow(recap_zoned) > 0) {
  movement_counts <- recap_zoned %>%
    count(tag_zone, recap_zone, name = "Freq") %>%
    filter(tag_zone != recap_zone, Freq >= 3)  # Only show movements with 3+ fish

  movement_arcs <- movement_counts %>%
    left_join(tag_centroids %>% rename(x = clon, y = clat), by = "tag_zone") %>%
    left_join(tag_centroids %>% rename(xend = clon, yend = clat),
              by = c("recap_zone" = "tag_zone")) %>%
    filter(!is.na(x) & !is.na(xend)) %>%
    mutate(
      delta_lon = abs(xend - x),
      delta_lat = abs(yend - y),
      direction = ifelse(delta_lon > delta_lat, "Along-shelf (E-W)", "Cross-shelf (N-S)")
    )

  fig4 <- ggplot() +
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
                                   "Cross-shelf (N-S)" = "#E6550D"),
                       name = "Movement axis") +
    geom_point(data = tag_centroids, aes(x = clon, y = clat),
               size = 1.5, color = "black") +
    geom_text(data = tag_centroids, aes(x = clon, y = clat, label = tag_zone),
              size = 2, vjust = -1, fontface = "bold") +
    coord_sf(xlim = c(-98, -80), ylim = c(24, 31.5), expand = FALSE) +
    labs(x = "Longitude", y = "Latitude") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_line(color = "gray95"),
          legend.position = "bottom", legend.box = "vertical")

  ggsave(file.path(OUTPUT_DIR, "Fig4_movement_arcs.png"), fig4,
         width = 170, height = 100, units = "mm", dpi = 300)
}

# =============================================================================
# FIGURE 5: Annotated week x zone heatmap
# =============================================================================
cat("  Fig 5: Annotated heatmap...\n")

# Weekly catch pooled across years for the heatmap background
weekly_catch_all <- mr_oisst %>%
  filter(!is.na(tag_zone) & tag_zone %in% valid_zones) %>%
  count(tag_zone, tag_week, name = "n")

# Peak annotations
peak_annot <- seasonal_peaks %>%
  group_by(tag_zone, season) %>%
  summarise(med_week = round(median(peak_week)), .groups = "drop")

fig5 <- ggplot(weekly_catch_all, aes(x = tag_week, y = factor(tag_zone), fill = n)) +
  geom_tile(color = NA) +
  geom_point(data = peak_annot %>% filter(season == "spring"),
             aes(x = med_week, y = factor(tag_zone)),
             inherit.aes = FALSE, shape = 25, size = 2.5,
             fill = "#2171B5", color = "white") +
  geom_point(data = peak_annot %>% filter(season == "fall"),
             aes(x = med_week, y = factor(tag_zone)),
             inherit.aes = FALSE, shape = 25, size = 2.5,
             fill = "#E6550D", color = "white") +
  scale_fill_gradient(low = "white", high = "#2D004B",
                      trans = "sqrt", name = "Tag count") +
  scale_x_continuous(breaks = c(1, 13, 26, 39, 52),
                     labels = c("Jan", "Apr", "Jul", "Oct", "Jan")) +
  geom_vline(xintercept = 25, linetype = "dashed", color = "gray60", linewidth = 0.3) +
  labs(x = "Week of year", y = "Statistical zone") +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank())

ggsave(file.path(OUTPUT_DIR, "Fig5_heatmap_annotated.png"), fig5,
       width = 170, height = 120, units = "mm", dpi = 300)

cat("  All publication figures saved.\n")
cat("07_pub_figures.R complete.\n\n")
