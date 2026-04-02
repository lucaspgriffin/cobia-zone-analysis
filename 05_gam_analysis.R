# =============================================================================
# 05_gam_analysis.R — GAM framework: what predicts peak timing?
# Competing models: SST vs photoperiod vs Julian day (latitude)
# =============================================================================

cat("Running GAM analysis...\n")

# --- Prepare data for GAMs ---------------------------------------------------
# Use yearly peaks as the unit of analysis
gam_data <- yearly_peaks %>%
  filter(!is.na(peak_sst) & !is.na(peak_daylength) & !is.na(circ_mean_julian)) %>%
  mutate(
    zone_f = factor(tag_zone),
    peak_week_scaled = peak_week  # keep original scale for cyclic splines
  )

cat("  GAM dataset:", nrow(gam_data), "zone-year observations\n")
cat("  Zones:", length(unique(gam_data$tag_zone)), "\n")
cat("  Years:", min(gam_data$tag_year), "-", max(gam_data$tag_year), "\n")

# --- Model 1: Peak week ~ SST (thermal tracking) ----------------------------
m_sst <- gam(peak_week ~ s(peak_sst, k = 6),
             data = gam_data, method = "REML")

# --- Model 2: Peak week ~ daylength (photoperiod cue) -----------------------
m_day <- gam(peak_week ~ s(peak_daylength, k = 6),
             data = gam_data, method = "REML")

# --- Model 3: Peak week ~ zone (fixed geographic/calendar effect) ------------
m_zone <- gam(peak_week ~ zone_f,
              data = gam_data, method = "REML")

# --- Model 4: Peak week ~ SST + year trend ----------------------------------
m_sst_yr <- gam(peak_week ~ s(peak_sst, k = 6) + s(tag_year, k = 6),
                data = gam_data, method = "REML")

# --- Model 5: SST + daylength (both environmental drivers) -------------------
m_both <- gam(peak_week ~ s(peak_sst, k = 6) + s(peak_daylength, k = 6),
              data = gam_data, method = "REML")

# --- Model 6: Full model with zone random effect ----------------------------
m_full <- gam(peak_week ~ s(peak_sst, k = 6) + s(peak_daylength, k = 6) +
                s(tag_year, k = 6) + s(zone_f, bs = "re"),
              data = gam_data, method = "REML")

# --- Model comparison (AIC) --------------------------------------------------
models <- list(
  "SST only"         = m_sst,
  "Daylength only"   = m_day,
  "Zone (fixed)"     = m_zone,
  "SST + year"       = m_sst_yr,
  "SST + daylength"  = m_both,
  "Full (SST+DL+yr+zone)" = m_full
)

aic_table <- data.frame(
  Model    = names(models),
  AIC      = sapply(models, AIC),
  BIC      = sapply(models, BIC),
  R2_adj   = sapply(models, function(m) summary(m)$r.sq),
  Dev_expl = sapply(models, function(m) summary(m)$dev.expl)
) %>%
  arrange(AIC) %>%
  mutate(
    delta_AIC = AIC - min(AIC),
    rank = row_number()
  )

cat("\n  MODEL COMPARISON:\n")
cat("  ", paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("  %4s %-28s %8s %8s %6s %8s\n",
            "Rank", "Model", "AIC", "dAIC", "R2adj", "DevExpl"))
cat("  ", paste(rep("-", 80), collapse = ""), "\n")
for (i in seq_len(nrow(aic_table))) {
  a <- aic_table[i, ]
  cat(sprintf("  %4d %-28s %8.1f %8.1f %5.1f%% %7.1f%%\n",
              a$rank, a$Model, a$AIC, a$delta_AIC,
              a$R2_adj * 100, a$Dev_expl * 100))
}

# --- Best model summary ------------------------------------------------------
best_name <- aic_table$Model[1]
best_model <- models[[best_name]]
cat("\n  Best model:", best_name, "\n")
cat("  Summary:\n")
print(summary(best_model))

# --- SST effect visualization ------------------------------------------------
cat("\n  Generating GAM plots...\n")

# Plot: SST partial effect from full model
p_sst_effect <- ggplot(gam_data, aes(x = peak_sst, y = peak_week)) +
  geom_point(aes(color = factor(tag_zone), size = n_tags), alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 6),
              color = "black", linewidth = 1.2, se = TRUE) +
  scale_color_viridis_d(option = "turbo", name = "Zone") +
  scale_size_continuous(name = "N Tags", range = c(1, 5)) +
  labs(title = "Peak Week vs SST at Peak (zone-year level)",
       subtitle = paste0("GAM SST-only: R2=", round(summary(m_sst)$r.sq, 3),
                         ", Dev.expl=", round(summary(m_sst)$dev.expl * 100, 1), "%"),
       x = "SST at Peak Week (C)", y = "Peak Week of Year") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUTPUT_DIR, "gam_sst_vs_peak_week.png"), p_sst_effect,
       width = 10, height = 7, dpi = 300)

# Plot: Daylength partial effect
p_day_effect <- ggplot(gam_data, aes(x = peak_daylength, y = peak_week)) +
  geom_point(aes(color = factor(tag_zone), size = n_tags), alpha = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 6),
              color = "black", linewidth = 1.2, se = TRUE) +
  scale_color_viridis_d(option = "turbo", name = "Zone") +
  scale_size_continuous(name = "N Tags", range = c(1, 5)) +
  labs(title = "Peak Week vs Daylength at Peak (zone-year level)",
       subtitle = paste0("GAM Daylength-only: R2=", round(summary(m_day)$r.sq, 3),
                         ", Dev.expl=", round(summary(m_day)$dev.expl * 100, 1), "%"),
       x = "Daylength at Peak (hours)", y = "Peak Week of Year") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUTPUT_DIR, "gam_daylength_vs_peak_week.png"), p_day_effect,
       width = 10, height = 7, dpi = 300)

# Plot: Year trend (has peak timing shifted?)
if (summary(m_sst_yr)$s.table["s(tag_year)", "p-value"] < 0.1) {
  cat("  Year trend is significant — peak timing may have shifted over time\n")
}

p_year <- ggplot(gam_data, aes(x = tag_year, y = peak_week)) +
  geom_point(aes(color = factor(tag_zone)), alpha = 0.4, size = 1.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 6),
              color = "black", linewidth = 1.2, se = TRUE) +
  scale_color_viridis_d(option = "turbo", name = "Zone") +
  labs(title = "Peak Week Over Time (phenological shift?)",
       subtitle = "GAM smooth ± SE across all zone-year observations",
       x = "Year", y = "Peak Week of Year") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(OUTPUT_DIR, "gam_year_trend.png"), p_year,
       width = 10, height = 7, dpi = 300)

# Save model comparison
write.csv(aic_table, file.path(OUTPUT_DIR, "gam_model_comparison.csv"), row.names = FALSE)

cat("\n05_gam_analysis.R complete.\n\n")
