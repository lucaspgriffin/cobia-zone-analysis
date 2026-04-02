# =============================================================================
# 03_sst_extract.R — Extract year-matched SST via NOAA OISST v2.1 ERDDAP
#
# Strategy: For each zone, find a valid OISST grid cell near the median tag
# location (some cells are NaN/land). Download monthly Gulf-wide composites
# for the full time period, then match to individual tags by zone+year+month.
#
# Using monthly resolution (not weekly) for reliability — SST doesn't change
# drastically week-to-week and this avoids 900+ individual ERDDAP queries.
# =============================================================================

cat("Extracting year-matched SST...\n")

sf_use_s2(FALSE)

SST_CACHE <- file.path(SHAPE_DIR, "zone_monthly_sst.csv")

# --- OISST grid helpers -------------------------------------------------------
snap_to_oisst <- function(x) {
  floor(x) + (round((x - floor(x)) * 4) / 4) + 0.125
}

# --- Compute median tag location per zone ------------------------------------
zone_tag_locs <- mr_oisst %>%
  filter(!is.na(tag_zone)) %>%
  group_by(tag_zone) %>%
  summarise(med_lat = median(tag_lat), med_lon = median(tag_lon),
            n_tags = n(), .groups = "drop")

if (!file.exists(SST_CACHE)) {
  cat("  Downloading monthly SST from NOAA OISST v2.1 ERDDAP...\n")

  erddap_base <- "https://coastwatch.pfeg.noaa.gov/erddap/griddap/ncdcOisst21Agg.csv"
  options(timeout = 120)

  # For each zone, find a valid ocean grid cell by testing nearby points
  # Then query the full time series at that point (1 year at a time, monthly stride)
  all_sst <- list()

  for (i in seq_len(nrow(zone_tag_locs))) {
    z <- zone_tag_locs[i, ]
    cat(sprintf("    Zone %2d (%.1f, %.1f)...", z$tag_zone, z$med_lat, z$med_lon))

    # Try grid cells: median, then nudge south/offshore in 0.25-deg steps
    candidates <- expand.grid(
      lat_offset = c(0, -0.25, -0.5, 0.25, -0.75),
      lon_offset = c(0, -0.25, 0.25, -0.5)
    )

    found <- FALSE
    for (ci in seq_len(nrow(candidates))) {
      test_lat <- snap_to_oisst(z$med_lat + candidates$lat_offset[ci])
      test_lon360 <- snap_to_oisst(z$med_lon + 360 + candidates$lon_offset[ci])

      # Quick test: get 1 day
      test_url <- paste0(erddap_base,
                         "?sst[(2020-06-15T12:00:00Z):1:(2020-06-15T12:00:00Z)]",
                         "[(0.0):1:(0.0)]",
                         "[(", test_lat, "):1:(", test_lat, ")]",
                         "[(", test_lon360, "):1:(", test_lon360, ")]")
      test_ok <- tryCatch({
        tmp <- tempfile(fileext = ".csv")
        download.file(test_url, tmp, quiet = TRUE, method = "libcurl")
        d <- read.csv(tmp, skip = 1, stringsAsFactors = FALSE)
        unlink(tmp)
        !is.na(as.numeric(d[1, ncol(d)]))
      }, error = function(e) FALSE)

      if (test_ok) {
        # Found a valid cell — now download monthly data year by year
        # Use stride of 30 to get ~monthly samples (reduces data volume)
        zone_data <- list()
        for (yr in 1982:2024) {
          url <- paste0(erddap_base,
                        "?sst[(", yr, "-01-01T12:00:00Z):30:(",
                        yr, "-12-31T12:00:00Z)]",
                        "[(0.0):1:(0.0)]",
                        "[(", test_lat, "):1:(", test_lat, ")]",
                        "[(", test_lon360, "):1:(", test_lon360, ")]")
          result <- tryCatch({
            tmp <- tempfile(fileext = ".csv")
            download.file(url, tmp, quiet = TRUE, method = "libcurl")
            d <- read.csv(tmp, skip = 1, stringsAsFactors = FALSE)
            names(d) <- c("time", "zlev", "latitude", "longitude", "sst")
            d$sst <- as.numeric(d$sst)
            d$date <- as.Date(substr(d$time, 1, 10))
            d$tag_zone <- z$tag_zone
            d$month <- month(d$date)
            d$year <- year(d$date)
            d <- d %>% select(tag_zone, date, year, month, sst) %>% filter(!is.na(sst))
            unlink(tmp)
            d
          }, error = function(e) NULL)
          if (!is.null(result) && nrow(result) > 0) {
            zone_data[[length(zone_data) + 1]] <- result
          }
          Sys.sleep(0.15)
        }

        zone_combined <- bind_rows(zone_data)
        if (nrow(zone_combined) > 0) {
          all_sst[[length(all_sst) + 1]] <- zone_combined
          cat(sprintf(" OK (%d obs, grid=%.3f,%.3f)\n",
                      nrow(zone_combined), test_lat, test_lon360))
          found <- TRUE
        }
        break
      }
      Sys.sleep(0.2)
    }

    if (!found) cat(" FAILED (no valid ocean cell found)\n")
  }

  sst_monthly_raw <- bind_rows(all_sst)

  # Compute monthly SST per zone-year
  sst_monthly <- sst_monthly_raw %>%
    group_by(tag_zone, year, month) %>%
    summarise(sst_mean = mean(sst, na.rm = TRUE), .groups = "drop")

  write.csv(sst_monthly, SST_CACHE, row.names = FALSE)
  cat("  Cached", nrow(sst_monthly), "monthly SST records\n")
  cat("  Zones with SST:", length(unique(sst_monthly$tag_zone)), "\n")

} else {
  cat("  Using cached SST:", SST_CACHE, "\n")
  sst_monthly <- read.csv(SST_CACHE)
}

# --- Join SST to each tag record by zone + year + month -----------------------
mr_oisst <- mr_oisst %>%
  left_join(sst_monthly, by = c("tag_zone", "tag_year" = "year", "tag_month" = "month"))

n_with_sst <- sum(!is.na(mr_oisst$sst_mean))
n_zoned <- sum(!is.na(mr_oisst$tag_zone))
cat("  Tags with year-matched SST:", n_with_sst, "of", n_zoned, "zoned tags\n")

if (n_with_sst > 0) {
  cat("  Mean SST at tagging:", round(mean(mr_oisst$sst_mean, na.rm = TRUE), 1), "C\n")
  q10 <- quantile(mr_oisst$sst_mean, 0.10, na.rm = TRUE)
  q90 <- quantile(mr_oisst$sst_mean, 0.90, na.rm = TRUE)
  cat("  80% of tagging:", round(q10, 1), "-", round(q90, 1), "C\n")
}

cat("03_sst_extract.R complete.\n\n")
