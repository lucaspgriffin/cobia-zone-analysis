# =============================================================================
# 05_gam_analysis.R — Season-specific GAMs: what predicts spring vs fall peaks?
# =============================================================================

cat("Running GAM analysis (spring vs fall)...\n")

# --- Run GAMs for each season separately + combined --------------------------
run_season_gams <- function(data, season_label) {
  d <- data %>%
    filter(!is.na(peak_sst) & !is.na(peak_daylength) & !is.na(circ_mean_julian)) %>%
    mutate(zone_f = factor(tag_zone))

  if (nrow(d) < 15) {
    cat("    ", season_label, ": insufficient data (n=", nrow(d), ")\n")
    return(NULL)
  }

  cat(sprintf("  %s: n=%d zone-year obs, %d zones, %d-%d\n",
              season_label, nrow(d), length(unique(d$tag_zone)),
              min(d$tag_year), max(d$tag_year)))

  # Fit competing models
  n_zones <- length(unique(d$zone_f))

  models <- list()
  models[["SST only"]] <- tryCatch(
    gam(peak_week ~ s(peak_sst, k = min(6, nrow(d) - 1)), data = d, method = "REML"),
    error = function(e) NULL)
  models[["Daylength only"]] <- tryCatch(
    gam(peak_week ~ s(peak_daylength, k = min(6, nrow(d) - 1)), data = d, method = "REML"),
    error = function(e) NULL)
  models[["SST + daylength"]] <- tryCatch(
    gam(peak_week ~ s(peak_sst, k = min(5, nrow(d) - 1)) +
          s(peak_daylength, k = min(5, nrow(d) - 1)), data = d, method = "REML"),
    error = function(e) NULL)
  models[["SST + year"]] <- tryCatch(
    gam(peak_week ~ s(peak_sst, k = min(5, nrow(d) - 1)) +
          s(tag_year, k = min(5, nrow(d) - 1)), data = d, method = "REML"),
    error = function(e) NULL)

  if (n_zones >= 3) {
    models[["Full (SST+DL+yr+zone)"]] <- tryCatch(
      gam(peak_week ~ s(peak_sst, k = min(5, nrow(d) - 1)) +
            s(peak_daylength, k = min(5, nrow(d) - 1)) +
            s(tag_year, k = min(5, nrow(d) - 1)) +
            s(zone_f, bs = "re"), data = d, method = "REML"),
      error = function(e) NULL)
  }

  # Remove NULLs
  models <- models[!sapply(models, is.null)]

  if (length(models) == 0) {
    cat("    All models failed\n")
    return(NULL)
  }

  # Compare
  aic_table <- data.frame(
    Model    = names(models),
    AIC      = sapply(models, AIC),
    R2_adj   = sapply(models, function(m) summary(m)$r.sq),
    Dev_expl = sapply(models, function(m) summary(m)$dev.expl),
    stringsAsFactors = FALSE
  ) %>%
    arrange(AIC) %>%
    mutate(delta_AIC = AIC - min(AIC), rank = row_number())

  cat(sprintf("\n    %-30s %8s %6s %8s\n", "Model", "dAIC", "R2adj", "DevExpl"))
  cat("    ", paste(rep("-", 60), collapse = ""), "\n")
  for (i in seq_len(nrow(aic_table))) {
    a <- aic_table[i, ]
    cat(sprintf("    %-30s %8.1f %5.1f%% %7.1f%%\n",
                a$Model, a$delta_AIC, a$R2_adj * 100, a$Dev_expl * 100))
  }

  list(models = models, aic_table = aic_table, data = d,
       best = models[[aic_table$Model[1]]])
}

# --- Spring GAMs --------------------------------------------------------------
cat("\n  === SPRING PEAKS ===\n")
spring_data <- seasonal_peaks %>% filter(season == "spring")
spring_result <- run_season_gams(spring_data, "SPRING")

# --- Fall GAMs ----------------------------------------------------------------
cat("\n  === FALL PEAKS ===\n")
fall_data <- seasonal_peaks %>% filter(season == "fall")
fall_result <- run_season_gams(fall_data, "FALL")

# --- Combined (all seasons) for comparison ------------------------------------
cat("\n  === COMBINED (all seasons) ===\n")
combined_result <- run_season_gams(seasonal_peaks, "COMBINED")

# --- Cross-season comparison --------------------------------------------------
cat("\n  CROSS-SEASON DRIVER COMPARISON:\n")
for (res_name in c("spring_result", "fall_result")) {
  res <- get(res_name)
  if (is.null(res)) next
  label <- gsub("_result", "", res_name)

  best <- res$aic_table$Model[1]
  r2 <- round(res$aic_table$R2_adj[1] * 100, 1)

  # Get individual model R2s for SST-only and DL-only
  sst_r2 <- if ("SST only" %in% res$aic_table$Model)
    round(res$aic_table$R2_adj[res$aic_table$Model == "SST only"] * 100, 1) else NA
  dl_r2 <- if ("Daylength only" %in% res$aic_table$Model)
    round(res$aic_table$R2_adj[res$aic_table$Model == "Daylength only"] * 100, 1) else NA

  cat(sprintf("  %7s: Best=%s (R2=%.1f%%) | SST-only R2=%s%% | DL-only R2=%s%%\n",
              toupper(label), best, r2,
              ifelse(is.na(sst_r2), "NA", sst_r2),
              ifelse(is.na(dl_r2), "NA", dl_r2)))
}

# --- Generate plots -----------------------------------------------------------
cat("\n  Generating GAM plots...\n")

# Plot: SST vs peak week, colored by season
p_sst <- ggplot(seasonal_peaks %>% filter(!is.na(peak_sst)),
                aes(x = peak_sst, y = peak_week, color = season)) +
  geom_point(aes(size = n_tags), alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("spring" = "#2171B5", "fall" = "#CB181D"),
                     name = "Season") +
  scale_size_continuous(name = "N Tags", range = c(1, 5)) +
  labs(title = "Peak SST vs Peak Week by Season",
       subtitle = "Separate GAM smooths — do spring and fall respond to SST differently?",
       x = "SST at Peak Week (C)", y = "Peak Week of Year") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR, "gam_sst_by_season.png"), p_sst,
       width = 10, height = 7, dpi = 300)

# Plot: Daylength vs peak week, colored by season
p_dl <- ggplot(seasonal_peaks %>% filter(!is.na(peak_daylength)),
               aes(x = peak_daylength, y = peak_week, color = season)) +
  geom_point(aes(size = n_tags), alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("spring" = "#2171B5", "fall" = "#CB181D"),
                     name = "Season") +
  scale_size_continuous(name = "N Tags", range = c(1, 5)) +
  labs(title = "Peak Daylength vs Peak Week by Season",
       x = "Daylength at Peak (hours)", y = "Peak Week of Year") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR, "gam_daylength_by_season.png"), p_dl,
       width = 10, height = 7, dpi = 300)

# Plot: Year trend by season
p_yr <- ggplot(seasonal_peaks, aes(x = tag_year, y = peak_week, color = season)) +
  geom_point(alpha = 0.4, size = 1.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = TRUE, linewidth = 1) +
  scale_color_manual(values = c("spring" = "#2171B5", "fall" = "#CB181D")) +
  labs(title = "Peak Week Over Time by Season (phenological shift?)",
       x = "Year", y = "Peak Week of Year") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(OUTPUT_DIR, "gam_year_trend_by_season.png"), p_yr,
       width = 10, height = 7, dpi = 300)

# Save model comparison tables
all_aic <- bind_rows(
  if (!is.null(spring_result)) spring_result$aic_table %>% mutate(season = "spring"),
  if (!is.null(fall_result)) fall_result$aic_table %>% mutate(season = "fall"),
  if (!is.null(combined_result)) combined_result$aic_table %>% mutate(season = "combined")
)
write.csv(all_aic, file.path(OUTPUT_DIR, "gam_model_comparison.csv"), row.names = FALSE)

cat("\n05_gam_analysis.R complete.\n\n")
