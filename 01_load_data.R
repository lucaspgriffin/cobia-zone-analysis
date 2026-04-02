# =============================================================================
# 01_load_data.R — Load, standardize, QA flag both M&R databases
# Restrict to OISST era (1982+) for year-matched SST analysis
# =============================================================================

cat("Loading mark-and-recapture data...\n")

# --- SEFSC CTC ---------------------------------------------------------------
sefsc_file <- file.path(DATA_DIR,
                        "250722 SEFSC CTC Cobia Data File NCEAS GEI.xlsx")
sefsc_data <- tryCatch({
  raw <- read_excel(sefsc_file, sheet = "COBIA")
  raw %>%
    filter(!is.na(LATITUDE_1) & !is.na(LONGITUDE_1)) %>%
    mutate(
      animal_id    = as.character(ANIMAL_ID),
      dataset      = "SEFSC_CTC",
      tag_lat      = as.numeric(LATITUDE_1),
      tag_lon      = as.numeric(LONGITUDE_1),
      tag_date     = as.Date(TAG_DATE_1),
      recap_lat    = as.numeric(LATITUDE_2),
      recap_lon    = as.numeric(LONGITUDE_2),
      recap_date   = as.Date(TAG_DATE_2),
      days_at_large = as.numeric(DAYS_AT_LARGE),
      length_cm    = ifelse(as.numeric(LENGTH_1) > 300,
                            as.numeric(LENGTH_1) / 10,
                            as.numeric(LENGTH_1)),
      recaptured   = !is.na(LATITUDE_2) & !is.na(LONGITUDE_2),
      usm_original_zone = NA_integer_
    ) %>%
    select(animal_id, dataset, tag_lat, tag_lon, tag_date,
           recap_lat, recap_lon, recap_date, days_at_large,
           length_cm, recaptured, usm_original_zone)
}, error = function(e) { warning("SEFSC load failed: ", e$message); data.frame() })

cat("  SEFSC CTC:", nrow(sefsc_data), "records\n")

# --- USM Cooperative ----------------------------------------------------------
usm_file <- file.path(DATA_DIR,
                      "USM Cooperative Sport Fish Tag and Release Cobia Database for GEI.xlsx")
usm_data <- tryCatch({
  raw <- read_excel(usm_file, sheet = "Cobia")
  raw %>%
    filter(!is.na(`Tag Lat`) & !is.na(`Tag Long`)) %>%
    mutate(
      animal_id    = paste0("USM-", as.character(TagNo)),
      dataset      = "USM_Cooperative",
      tag_lat      = as.numeric(`Tag Lat`),
      tag_lon      = as.numeric(`Tag Long`),
      tag_date     = as.Date(TagDate),
      recap_lat    = as.numeric(`Rec Lat`),
      recap_lon    = as.numeric(`Rec Long`),
      recap_date   = as.Date(RecDate),
      days_at_large = as.numeric(DaysOut),
      length_cm    = as.numeric(TagTL),
      recaptured   = !is.na(`Rec Lat`) & !is.na(`Rec Long`),
      usm_original_zone = as.integer(TagZN)
    ) %>%
    select(animal_id, dataset, tag_lat, tag_lon, tag_date,
           recap_lat, recap_lon, recap_date, days_at_large,
           length_cm, recaptured, usm_original_zone)
}, error = function(e) { warning("USM load failed: ", e$message); data.frame() })

cat("  USM Cooperative:", nrow(usm_data), "records\n")

# --- Combine, filter, add temporal + photoperiod + QA columns ----------------
mr_data <- bind_rows(sefsc_data, usm_data) %>%
  filter(
    tag_lat >= 24 & tag_lat <= 31,
    tag_lon >= -98 & tag_lon <= -80
  ) %>%
  mutate(
    tag_year   = year(tag_date),
    tag_month  = month(tag_date),
    tag_week   = isoweek(tag_date),
    tag_julian = yday(tag_date),
    recap_julian = ifelse(recaptured, yday(recap_date), NA_integer_),
    recap_week   = ifelse(recaptured, isoweek(recap_date), NA_integer_),
    recap_year   = ifelse(recaptured, year(recap_date), NA_integer_),
    distance_km  = ifelse(recaptured,
                          calculate_distance_km(tag_lat, tag_lon, recap_lat, recap_lon),
                          NA_real_),
    # Photoperiod at tag location and date
    daylength_hrs = calc_daylength(tag_lat, tag_julian),
    # QA flags
    qa_pre_oisst       = tag_year < 1982,
    qa_length_suspect  = (!is.na(length_cm) & (length_cm < 10 | length_cm > 200)),
    qa_coord_imprecise = (round(tag_lat, 1) == tag_lat & round(tag_lon, 1) == tag_lon),
    # Snap to OISST 0.25-deg grid for SST extraction
    sst_grid_lat = round(tag_lat * 4) / 4,
    sst_grid_lon = round(tag_lon * 4) / 4
  )

# Restrict primary analysis to OISST era
mr_oisst <- mr_data %>% filter(!qa_pre_oisst)

n_total <- nrow(mr_data)
n_oisst <- nrow(mr_oisst)
cat("  Combined (Gulf-filtered):", n_total, "\n")
cat("  OISST era (1982+):", n_oisst, "(primary analysis set)\n")
cat("  Pre-1982 (excluded from SST):", sum(mr_data$qa_pre_oisst), "\n")
cat("  Recaptured:", sum(mr_oisst$recaptured), "\n")
cat("  Date range:", as.character(min(mr_oisst$tag_date, na.rm = TRUE)),
    "to", as.character(max(mr_oisst$tag_date, na.rm = TRUE)), "\n")

cat("01_load_data.R complete.\n\n")
