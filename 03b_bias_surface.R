# =============================================================================
# 03b_bias_surface.R
# -----------------------------------------------------------------------------
# Purpose : Create a kernel density bias surface from Aves observer locations
#           in the training AOI (TN+KA+KL). Used in Script 05 to weight
#           background point sampling proportional to observer effort.
#
# RATIONALE:
#   GBIF occurrence records are spatially biased toward areas with more
#   observer effort (bird sanctuaries, accessible wetlands, urban areas).
#   If background points are sampled uniformly at random, the model learns
#   "species prefers observer-accessible areas" rather than actual habitat.
#
#   Solution: sample background points proportional to observer effort,
#   so background represents "where observers looked" rather than
#   "everywhere equally". This directly improves CBI because background
#   and presence share the same observer effort distribution.
#
# Method:
#   1. Load bg_optionC.csv — 3,639 unique Aves observer locations in AOI
#   2. Create kernel density raster using terra::rasterize + focal smoothing
#   3. Normalise to 0-1 probability surface
#   4. Save as bias_surface.tif for use in Script 05
#
# Input   :
#   background/bg_optionC.csv
#   data/processed/bioclim/bioclim_train_2000_2010.tif  (for grid template)
#
# Output  :
#   data/processed/background/bias_surface.tif
#
# Author  : [Your name]
# Date    : 2026-03-28
# =============================================================================


# ---------------------------------------------------------------------------
# 0. PACKAGES
# ---------------------------------------------------------------------------

required_pkgs <- c("terra", "here")
new_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[,"Package"]]
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)

library(terra)
library(here)

options(scipen = 999)


# ---------------------------------------------------------------------------
# 1. PATHS
# ---------------------------------------------------------------------------

aves_bg_file   <- here("background", "bg_optionC.csv")
bio_train_file <- here("data", "processed", "bioclim",
                       "bioclim_train_2000_2010.tif")
out_dir        <- here("data", "processed", "background")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ---------------------------------------------------------------------------
# 2. LOAD INPUTS
# ---------------------------------------------------------------------------

message("Loading Aves observer locations...")
aves_pts <- read.csv(aves_bg_file, stringsAsFactors = FALSE)
message("  Observer locations: ", nrow(aves_pts))
message("  Lon range: ", round(min(aves_pts$decimalLongitude), 3),
        " to ", round(max(aves_pts$decimalLongitude), 3))
message("  Lat range: ", round(min(aves_pts$decimalLatitude),  3),
        " to ", round(max(aves_pts$decimalLatitude),  3))

message("Loading BIO template raster...")
bio_template <- terra::rast(bio_train_file)[[1]]
message("  Resolution: ", terra::res(bio_template)[1], " deg")
message("  Extent    : ", as.character(terra::ext(bio_template)))


# ---------------------------------------------------------------------------
# 3. RASTERIZE OBSERVER LOCATIONS
#    Count number of observer locations per 1km pixel.
#    This gives raw observation density before smoothing.
# ---------------------------------------------------------------------------

message("\nRasterizing observer locations...")

aves_v <- terra::vect(aves_pts,
                      geom = c("decimalLongitude", "decimalLatitude"),
                      crs  = "EPSG:4326")

# Count points per cell
density_raw <- terra::rasterize(aves_v, bio_template,
                                 fun = "count", background = 0)
names(density_raw) <- "observer_count"

n_cells_with_obs <- sum(terra::values(density_raw) > 0, na.rm = TRUE)
message("  Cells with at least one observation: ", n_cells_with_obs,
        " / ", terra::ncell(density_raw))


# ---------------------------------------------------------------------------
# 4. KERNEL DENSITY SMOOTHING
#    Apply Gaussian focal smoothing to spread observer effort to
#    neighbouring pixels. This prevents sampling being restricted only
#    to exact observer location pixels and creates a smooth probability
#    surface representing observer accessibility.
#
#    Kernel size: 5x5 window (~5km at 1km resolution)
#    This is a conservative smoothing — large enough to fill gaps between
#    nearby observation locations but not so large that it washes out
#    the spatial structure of observer effort.
# ---------------------------------------------------------------------------

message("\nApplying kernel density smoothing (5x5 Gaussian kernel)...")

# Gaussian weights for 5x5 window
gauss_weights <- matrix(
  c(1, 2, 3, 2, 1,
    2, 4, 6, 4, 2,
    3, 6, 9, 6, 3,
    2, 4, 6, 4, 2,
    1, 2, 3, 2, 1),
  nrow = 5, ncol = 5
)

bias_smooth <- terra::focal(
  density_raw,
  w      = gauss_weights,
  fun    = "sum",
  na.rm  = TRUE,
  pad    = TRUE,
  padValue = 0
)

message("  Smoothing complete.")


# ---------------------------------------------------------------------------
# 5. NORMALISE TO 0-1 PROBABILITY SURFACE
#    Add small offset (1) before normalising so zero-observation cells
#    still have a small non-zero probability of being sampled.
#    This prevents complete exclusion of unobserved areas while still
#    strongly upweighting high observer effort areas.
# ---------------------------------------------------------------------------

message("\nNormalising to probability surface...")

# Add offset so zero-obs cells aren't completely excluded
# Divide by max only (not min-max) to preserve the offset —
# min-max normalisation would map the offset back to zero, defeating its purpose
bias_offset <- bias_smooth + 1

bias_max  <- terra::global(bias_offset, "max", na.rm = TRUE)[[1]]
bias_norm <- bias_offset / bias_max
names(bias_norm) <- "bias_probability"

message("  Bias surface range: ",
        round(terra::global(bias_norm, "min", na.rm=TRUE)[[1]], 4),
        " to ",
        round(terra::global(bias_norm, "max", na.rm=TRUE)[[1]], 4))


# ---------------------------------------------------------------------------
# 6. MASK TO LAND SURFACE
#    Use GLC 2005 as land mask — same as used for background sampling
#    in Script 05. Ocean/water pixels should have zero probability.
# ---------------------------------------------------------------------------

message("\nMasking to land surface...")

glc_2005 <- terra::rast(
  here("GEE_Exports", "SDM_GLC_Yearly", "GLC_2005.tif")
)

# Align GLC to BIO grid
glc_aligned <- terra::resample(glc_2005, bio_template, method = "near")
land_mask   <- terra::ifel(glc_aligned > 0, 1, NA)

bias_masked <- terra::mask(bias_norm, land_mask)
names(bias_masked) <- "bias_probability"

message("  Non-NA cells in bias surface: ",
        sum(!is.na(terra::values(bias_masked))))


# ---------------------------------------------------------------------------
# 7. SAVE
# ---------------------------------------------------------------------------

out_path <- file.path(out_dir, "bias_surface.tif")
terra::writeRaster(bias_masked, out_path,
                   overwrite = TRUE, datatype = "FLT4S")

message("\nSaved: data/processed/background/bias_surface.tif")

# Quick summary
message("\nBias surface summary:")
message("  High probability areas = dense observer coverage")
message("  Low probability areas  = sparse observer coverage")
message("  All land pixels included (offset prevents zero probability)")
message("  Background sampling in Script 05 will be weighted by this surface")

message("\nScript 04d complete. Next: 05_ENMeval_tuning.R")
