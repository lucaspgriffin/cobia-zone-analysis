# =============================================================================
# 00_setup.R — Packages, paths, utilities, zone data download
# Cobia M&R NOAA Zone Analysis Pipeline v2
# =============================================================================

# --- Packages ----------------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(ggridges)
library(scales)
library(viridis)
library(sf)
library(circular)
library(mgcv)

# --- Paths -------------------------------------------------------------------
DATA_DIR <- file.path(dirname(getwd()), "florida-fish-tagging",
                      "data", "GEI", "Cobia Mark and Recapture")
if (!dir.exists(DATA_DIR)) {
  DATA_DIR <- Sys.getenv("COBIA_DATA_DIR", unset = DATA_DIR)
  if (!dir.exists(DATA_DIR)) stop("Cannot find data directory: ", DATA_DIR,
    "\n  Set COBIA_DATA_DIR env var or place data alongside this repo.")
}
SHAPE_DIR  <- "data"
OUTPUT_DIR <- "output"

dir.create(SHAPE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Haversine distance ------------------------------------------------------
calculate_distance_km <- function(lat1, lon1, lat2, lon2) {
  lat1_rad <- lat1 * pi / 180; lon1_rad <- lon1 * pi / 180
  lat2_rad <- lat2 * pi / 180; lon2_rad <- lon2 * pi / 180
  dlat <- lat2_rad - lat1_rad; dlon <- lon2_rad - lon1_rad
  a <- sin(dlat/2)^2 + cos(lat1_rad) * cos(lat2_rad) * sin(dlon/2)^2
  c <- 2 * atan2(sqrt(a), sqrt(1-a))
  return(6371 * c)
}

# --- Photoperiod calculation -------------------------------------------------
# Compute day length (hours) from latitude and Julian day
# Uses CBM model (Forsythe et al. 1995)
calc_daylength <- function(lat, julian_day) {
  P <- asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (julian_day - 186)))))
  D <- 24 - (24 / pi) * acos(
    (sin(0.8333 * pi / 180) + sin(lat * pi / 180) * sin(P)) /
    (cos(lat * pi / 180) * cos(P))
  )
  return(D)
}

# --- Download & cache NOAA statistical zone polygons -------------------------
ZONE_CACHE <- file.path(SHAPE_DIR, "noaa_stat_zones.geojson")
FEATURESERVER_BASE <- "https://services1.arcgis.com/qr14biwnHA6Vis6l/arcgis/rest/services/Southeast_US_Statistical_Zones/FeatureServer"

if (!file.exists(ZONE_CACHE)) {
  cat("Downloading NOAA statistical zone polygons...\n")
  zone_layers <- list()
  for (layer_id in 0:3) {
    url <- paste0(FEATURESERVER_BASE, "/", layer_id,
                  "/query?where=1%3D1&outFields=*&f=geojson&outSR=4326")
    tmp_file <- file.path(SHAPE_DIR, paste0("zones_layer_", layer_id, ".geojson"))
    download.file(url, tmp_file, quiet = TRUE)
    zone_layers[[layer_id + 1]] <- st_read(tmp_file, quiet = TRUE)
  }
  all_zones <- do.call(rbind, zone_layers)
  sf_use_s2(FALSE)
  all_zones <- st_make_valid(all_zones)
  st_write(all_zones, ZONE_CACHE, delete_dsn = TRUE, quiet = TRUE)
  cat("  Cached", nrow(all_zones), "zones\n")
} else {
  cat("Using cached zone polygons\n")
}

cat("00_setup.R complete.\n\n")
