# =============================================================================
# 03c_process_test_background.R
# -----------------------------------------------------------------------------
# Purpose : Read the massive 6GB+ GBIF Aves download, extract 100,000 unique
#           spatial locations, and generate a Kernel Density (KDE) bias
#           surface for the 2011-2021 testing period.
# =============================================================================

# ---------------------------------------------------------------------------
# 0. PACKAGES
# ---------------------------------------------------------------------------
required_pkgs <- c("data.table", "dplyr", "terra", "MASS", "here")
new_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)

library(data.table)
library(dplyr)
library(terra)
library(MASS)
library(here)

set.seed(42)

# ---------------------------------------------------------------------------
# 1. CONFIGURATION
# ---------------------------------------------------------------------------
# Path to the large GBIF Aves background CSV (project-relative via here())
gbif_massive_csv <- here("GBIF_Background_Aves_India_2011_2021",
                         "GBIF_Background_Aves_India_2011_2021.csv")

# Reference raster (to make sure the bias surface matches your bioclim extent/resolution)
ref_raster_file <- here("data", "processed", "bioclim", "bioclim_test_2011_2021.tif")

out_dir <- here("data", "processed", "background")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 2. FAST READ AND THIN (THE 6GB BEAST)
# ---------------------------------------------------------------------------
message("Reading the 6GB GBIF file... (this will take 10-30 seconds)")

# fread is insanely fast. We only select the lat/lon columns to save RAM.
aves_massive <- fread(
  file = gbif_massive_csv, 
  select = c("decimalLongitude", "decimalLatitude"),
  na.strings = c("", "NA")
)

message("  Initial rows loaded: ", format(nrow(aves_massive), big.mark = ","))

# Filter NAs and get distinct spatial points
aves_unique <- aves_massive %>%
  filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
  distinct(decimalLongitude, decimalLatitude)

message("  Unique spatial locations: ", format(nrow(aves_unique), big.mark = ","))

# Sample 100,000 points
if(nrow(aves_unique) > 100000) {
  aves_sampled <- aves_unique %>% sample_n(100000)
} else {
  aves_sampled <- aves_unique
  message("  Note: Less than 100k unique points found, using all of them.")
}

# Save the thinned CSV
sampled_csv_path <- file.path(out_dir, "gbif_100k_test_2011_2021.csv")
write.csv(aves_sampled, sampled_csv_path, row.names = FALSE)
message("  Saved 100k sample to: ", sampled_csv_path)

# Free up RAM
rm(aves_massive, aves_unique)
gc()

# ---------------------------------------------------------------------------
# 3. GENERATE BIAS SURFACE (KERNEL DENSITY)
# ---------------------------------------------------------------------------
message("\nGenerating Bias Surface Raster...")

ref_raster <- terra::rast(ref_raster_file)[[1]]
ext_ref    <- terra::ext(ref_raster)

# Run 2D Kernel Density Estimation
# Using n=200 for a smooth, high-resolution density calculation
kde <- MASS::kde2d(
  x = aves_sampled$decimalLongitude,
  y = aves_sampled$decimalLatitude,
  n = 200, 
  lims = c(ext_ref[1], ext_ref[2], ext_ref[3], ext_ref[4])
)

# Convert KDE to a terra raster
bias_rast <- terra::rast(list(x = kde$x, y = kde$y, z = kde$z))
terra::crs(bias_rast) <- terra::crs(ref_raster)

# Resample to match your reference Bioclim grid perfectly
bias_aligned <- terra::resample(bias_rast, ref_raster, method = "bilinear")

# Normalize (Divide by Max to keep values relative to highest effort)
bias_norm <- bias_aligned / terra::global(bias_aligned, "max", na.rm=TRUE)[1,1]
names(bias_norm) <- "bias_weight"

# Save the final bias surface
bias_tif_path <- file.path(out_dir, "bias_surface_test_2011_2021.tif")
terra::writeRaster(bias_norm, bias_tif_path, overwrite = TRUE, datatype = "FLT4S")

message("  Saved normalized bias surface to: ", bias_tif_path)

message("\n", strrep("=", 60))
message("SUCCESS. You are ready to update Script 07 with the new bias surface.")
message(strrep("=", 60))