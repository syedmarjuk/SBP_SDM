# =============================================================================
# 04b_clip_test_to_study_area.R
# -----------------------------------------------------------------------------
# Purpose : Clip test occurrence points to the Tamil Nadu + 10km coastal
#           buffer study area shapefile.
#
# Input   :
#   data/processed/test_thinned_with_year.csv   (n=872, from Script 01b)
#   FINAL_STUDY_AREAWGS84/FINAL_STUDY_AREAWGS84/
#     study_area_tn_coast_10km_wgs84.shp
#
# Output  :
#   data/processed/test_thinned_TN.csv   (points inside study area only)
#
# Note    : This filtered file is used in Script 07 evaluation only.
#           test_thinned_with_year.csv is NOT overwritten.
# =============================================================================


# ---------------------------------------------------------------------------
# 0. PACKAGES
# ---------------------------------------------------------------------------

required_pkgs <- c("sf", "dplyr", "here")
new_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)

library(sf)
library(dplyr)
library(here)


# ---------------------------------------------------------------------------
# 1. PATHS
# ---------------------------------------------------------------------------

shapefile_path <- here(
  "FINAL_STUDY_AREAWGS84", "FINAL_STUDY_AREAWGS84",
  "study_area_tn_coast_10km_wgs84.shp"
)

test_input  <- here("data", "processed", "test_thinned_with_year.csv")
test_output <- here("data", "processed", "test_thinned_TN.csv")


# ---------------------------------------------------------------------------
# 2. INPUT VALIDATION
# ---------------------------------------------------------------------------

if (!file.exists(shapefile_path)) stop("Shapefile not found:\n", shapefile_path)
if (!file.exists(test_input))     stop("Missing: ", test_input,
                                       "\nRun 01b first.")


# ---------------------------------------------------------------------------
# 3. LOAD SHAPEFILE
# ---------------------------------------------------------------------------

message("Loading study area shapefile...")
study_area <- sf::st_read(shapefile_path, quiet = TRUE)

if (sf::st_crs(study_area)$epsg != 4326) {
  message("  Reprojecting shapefile to EPSG:4326...")
  study_area <- sf::st_transform(study_area, 4326)
}

message("  CRS: EPSG:4326")
message("  Features: ", nrow(study_area))


# ---------------------------------------------------------------------------
# 4. LOAD TEST POINTS
# ---------------------------------------------------------------------------

message("Loading test occurrences...")
test_occ <- read.csv(test_input, stringsAsFactors = FALSE)
message("  Points loaded: ", nrow(test_occ))

if (!all(c("decimalLongitude", "decimalLatitude") %in% names(test_occ))) {
  stop("CSV must have decimalLongitude and decimalLatitude columns.")
}


# ---------------------------------------------------------------------------
# 5. SPATIAL FILTER
# ---------------------------------------------------------------------------

message("Clipping to study area...")

test_sf <- sf::st_as_sf(
  test_occ,
  coords = c("decimalLongitude", "decimalLatitude"),
  crs    = 4326,
  remove = FALSE
)

inside     <- sf::st_intersects(test_sf, study_area, sparse = FALSE)
inside_any <- apply(inside, 1, any)

test_clipped <- test_occ[inside_any, ]

message("  Before : ", nrow(test_occ))
message("  After  : ", nrow(test_clipped))
message("  Removed: ", nrow(test_occ) - nrow(test_clipped))


# ---------------------------------------------------------------------------
# 6. VALIDATION
# ---------------------------------------------------------------------------

if (nrow(test_clipped) == 0) {
  stop("No points remain. Check CRS alignment between shapefile and CSV.")
}

if (nrow(test_clipped) == nrow(test_occ)) {
  warning("All points passed — none were removed. Verify the shapefile boundary.")
}

message("  Year distribution:")
print(table(test_clipped$year))


# ---------------------------------------------------------------------------
# 7. WRITE OUTPUT
# ---------------------------------------------------------------------------

write.csv(test_clipped, test_output, row.names = FALSE)
message("\nSaved: data/processed/test_thinned_TN.csv  [n = ",
        nrow(test_clipped), "]")
message("Script 04b complete. Next: 04_VIF_predictor_selection.R")