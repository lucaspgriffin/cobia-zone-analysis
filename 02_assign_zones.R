# =============================================================================
# 02_assign_zones.R — Assign NOAA statistical zones via spatial join
# Cobia M&R NOAA Zone Analysis Pipeline
# =============================================================================

cat("Assigning NOAA statistical zones...\n")

sf_use_s2(FALSE)

# --- Load zone polygons ------------------------------------------------------
zones <- st_read(ZONE_CACHE, quiet = TRUE) %>%
  st_make_valid() %>%
  select(StatZone, Region, geometry)

cat("  Loaded", nrow(zones), "zone polygons\n")

# --- Create sf points for tag locations --------------------------------------
tag_pts <- mr_data %>%
  filter(!is.na(tag_lat) & !is.na(tag_lon)) %>%
  st_as_sf(coords = c("tag_lon", "tag_lat"), crs = 4326, remove = FALSE)

# Spatial join: tag locations -> zones
tag_joined <- st_join(tag_pts, zones, join = st_intersects, left = TRUE)

# Extract zone assignments back to the main data
# Use distinct() to avoid many-to-many joins from duplicate animal_ids
tag_zone_lookup <- tag_joined %>%
  st_drop_geometry() %>%
  select(animal_id, tag_zone = StatZone) %>%
  distinct(animal_id, .keep_all = TRUE)

mr_data <- mr_data %>%
  left_join(tag_zone_lookup, by = "animal_id")

# --- Assign zones to recapture locations too ---------------------------------
recap_subset <- mr_data %>%
  filter(recaptured & !is.na(recap_lat) & !is.na(recap_lon))

if (nrow(recap_subset) > 0) {
  recap_pts <- recap_subset %>%
    st_as_sf(coords = c("recap_lon", "recap_lat"), crs = 4326, remove = FALSE)

  recap_joined <- st_join(recap_pts, zones, join = st_intersects, left = TRUE)

  recap_zone_lookup <- recap_joined %>%
    st_drop_geometry() %>%
    select(animal_id, recap_zone = StatZone) %>%
    distinct(animal_id, .keep_all = TRUE)

  mr_data <- mr_data %>%
    left_join(recap_zone_lookup, by = "animal_id")
} else {
  mr_data$recap_zone <- NA_integer_
}

# --- Summary -----------------------------------------------------------------
n_total <- nrow(mr_data)
n_assigned <- sum(!is.na(mr_data$tag_zone))
n_na <- sum(is.na(mr_data$tag_zone))
pct_na <- round(100 * n_na / n_total, 1)

cat("\n  ZONE ASSIGNMENT SUMMARY\n")
cat("  -----------------------\n")
cat("  Total records:", n_total, "\n")
cat("  Assigned to zone:", n_assigned, "(", round(100 * n_assigned / n_total, 1), "%)\n")
cat("  No zone (NA):", n_na, "(", pct_na, "%)\n")

# Tags per zone
zone_counts <- mr_data %>%
  filter(!is.na(tag_zone)) %>%
  count(tag_zone) %>%
  arrange(tag_zone)
cat("\n  Tags per zone:\n")
for (i in seq_len(nrow(zone_counts))) {
  cat("    Zone", sprintf("%2d", zone_counts$tag_zone[i]), ":",
      sprintf("%5d", zone_counts$n[i]), "\n")
}

# --- USM zone concordance validation -----------------------------------------
usm_subset <- mr_data %>%
  filter(dataset == "USM_Cooperative" & !is.na(usm_original_zone) & !is.na(tag_zone))

if (nrow(usm_subset) > 0) {
  cat("\n  USM ZONE CONCORDANCE (usm_original_zone vs NOAA tag_zone):\n")
  concordance <- table(USM_Zone = usm_subset$usm_original_zone,
                       NOAA_Zone = usm_subset$tag_zone)
  print(concordance)
}

# Recapture zone summary
n_recap <- sum(mr_data$recaptured)
n_recap_zoned <- sum(!is.na(mr_data$recap_zone))
cat("\n  Recaptures:", n_recap, "| With zone:", n_recap_zoned, "\n")

# Make tag_zone an ordered factor for plotting (S to N, then W)
mr_data$tag_zone_f <- factor(mr_data$tag_zone, levels = 1:21)

# Update the OISST-era subset with zone assignments
mr_oisst <- mr_data %>% filter(!qa_pre_oisst)

cat("\n02_assign_zones.R complete.\n\n")
