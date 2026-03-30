# =============================================================================
# 02_bioclim_computation.R
# -----------------------------------------------------------------------------
# Purpose : Compute period-mean BIO1-19 bioclimatic variables from monthly
#           CHELSA data for training (2000-2010) and testing (2011-2021).
#
# METHOD — Period climatology (monthly means → single biovars call):
#   For each of the 12 calendar months, compute the pixel-wise mean across
#   all years in the period (e.g. mean of all January values 2000-2010).
#   Pass those 12 period-mean monthly rasters into dismo::biovars() ONCE
#   per period. This produces a true period climatology.
#
# DESIGN DECISION — why period means for climate, annual for non-climatic:
#   BIO variables are climatological summaries by definition. BIO1 = annual
#   mean temperature, BIO12 = annual precipitation — these are meaningless
#   for a single year and are appropriately averaged across the period.
#   This matches how CHELSA itself computes bioclim products (30-year means)
#   and is what published dynamic SDM papers do.
#   The dynamic innovation is in the NON-CLIMATIC predictors (GLC, GWL,
#   JRC, DistCoastalWet) which ARE matched to each year individually.
#   Defense: "Bioclimatic variables represent long-term climatic conditions
#   and are summarised as period means, while land cover, wetland, and
#   surface water predictors were temporally matched to each year to capture
#   landscape dynamics."
#
# WHY dismo::biovars() NOT terra::biovars():
#   terra does NOT export biovars() — it does not exist in terra.
#   dismo::biovars() is the correct function.
#   SpatRasters are converted to raster::Brick before calling dismo,
#   then converted back to SpatRaster for writing.
#
# CHELSA V2.1 scaling:
#   tasmin, tasmax : K×10 → divide by 10, subtract 273.15 → °C
#   pr             : mm×100 → divide by 100 → mm/month
#   Pre-run diagnostic detects whether terra auto-applied scaling
#   (double-scaling trap) and sets SCALING_NEEDED flag accordingly.
#
# Inputs  :
#   Chelsa/CHELSA_monthly/
#     TRAINING_2000_2010/{tasmin,tasmax,pr}/  (132 files each)
#     TESTING_2011_2021/{tasmin,tasmax,pr}/   (132 files each)
#   GEE_Exports/SDM_GLC_Yearly/GLC_2000.tif
#     (land mask — REQUIRED, hard stop if missing)
#
# Outputs :
#   data/processed/bioclim/
#     bioclim_train_2000_2010.tif   — 19-band period mean, training
#     bioclim_test_2011_2021.tif    — 19-band period mean, testing
#
# Author  : [Your name]
# Date    : 2026-03-27
# =============================================================================


# ---------------------------------------------------------------------------
# 0. PACKAGES
# ---------------------------------------------------------------------------

required_pkgs <- c("terra", "dismo", "raster", "here")
new_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[,"Package"]]
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)

# Load raster BEFORE terra — terra must load last so its functions
# (mask, resample, project, crs, values, writeRaster etc.) shadow raster's.
library(raster)
library(terra)
library(dismo)
library(here)

# Terra temp file management — prevents C: drive filling silently
terraOptions(memfrac = 0.6)
terra::tmpFiles(current = FALSE, orphan = TRUE, old = TRUE, remove = TRUE)
message("terra temp files cleared.")

# Raster temp file management — dismo::biovars() uses raster engine internally
# and writes to a SEPARATE temp folder that terra::tmpFiles() cannot see.
rasterOptions(tmpdir = tempdir())
raster::removeTmpFiles(h = 0)
message("raster temp files cleared.")


# ---------------------------------------------------------------------------
# 1. PATHS AND CONSTANTS
# ---------------------------------------------------------------------------

chelsa_root <- here("Chelsa", "CHELSA_monthly")
gee_glc_dir <- here("GEE_Exports", "SDM_GLC_Yearly")

out_root <- here("data", "processed", "bioclim")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

MONTHS    <- sprintf("%02d", 1:12)
BIO_NAMES <- paste0("BIO", 1:19)

periods <- list(
  list(
    label         = "Training 2000-2010",
    folder        = "TRAINING_2000_2010",
    years         = 2000:2010,
    expected_yrs  = 11L,
    out_file      = file.path(out_root, "bioclim_train_2000_2010.tif")
  ),
  list(
    label         = "Testing 2011-2021",
    folder        = "TESTING_2011_2021",
    years         = 2011:2021,
    expected_yrs  = 11L,
    out_file      = file.path(out_root, "bioclim_test_2011_2021.tif")
  )
)


# ---------------------------------------------------------------------------
# 2. PRE-RUN SCALING DIAGNOSTICS — temperature and precipitation separately
#    CHELSA V2.1 files are consistently encoded, but terra may auto-apply
#    scale/offset metadata differently per variable. Diagnosing separately
#    prevents silently wrong precipitation if terra handles metadata
#    differently for temp vs precip files.
# ---------------------------------------------------------------------------

message("\n", strrep("=", 60))
message("PRE-RUN SCALING DIAGNOSTICS")
message(strrep("=", 60))

# --- Temperature diagnostic (tasmax, January — any month works; January chosen
#     as it is consistently available and values are well within expected range) ---
temp_diag_file <- file.path(chelsa_root, "TRAINING_2000_2010", "tasmax",
                             "CHELSA_tasmax_01_2000_AOI.tif")
if (!file.exists(temp_diag_file)) {
  stop("Temperature diagnostic file not found:\n", temp_diag_file,
       "\nCheck CHELSA downloads completed.")
}

t_vals <- as.vector(values(rast(temp_diag_file), na.rm = TRUE))
t_min  <- round(min(t_vals), 1)
t_max  <- round(max(t_vals), 1)
message("tasmax Jan 2000 raw values: ", t_min, " to ", t_max)

# Three-way detection — terra may apply scale factor only (÷10) without
# the offset (-273.15), yielding pure Kelvin (~270-320) rather than K×10
# (~2700-3200) or °C (~15-50). All three cases must be handled.
if (t_min > 2000 & t_max < 4000) {
  TEMP_SCALE_TYPE <- "K10"
  message("TEMP: K×10 detected — (r/10) - 273.15 will be applied.")
} else if (t_min > 250 & t_max < 350) {
  TEMP_SCALE_TYPE <- "K"
  message("TEMP: Pure Kelvin detected — r - 273.15 will be applied.")
} else if (t_min > -60 & t_max < 70) {
  TEMP_SCALE_TYPE <- "C"
  message("TEMP: Already in °C — no scaling needed.")
} else {
  stop("TEMP: Unexpected values ", t_min, " to ", t_max,
       "\nExpected K×10 (~2700-3200), pure Kelvin (~250-350), or °C (~5-50).",
       "\nCheck CHELSA files for corruption.")
}
rm(t_vals); gc()

# --- Precipitation diagnostic (pr, July — peak monsoon, highest values) ---
prec_diag_file <- file.path(chelsa_root, "TRAINING_2000_2010", "pr",
                             "CHELSA_pr_07_2000_AOI.tif")
if (!file.exists(prec_diag_file)) {
  stop("Precipitation diagnostic file not found:\n", prec_diag_file,
       "\nCheck CHELSA downloads completed.")
}

p_vals <- as.vector(values(rast(prec_diag_file), na.rm = TRUE))
p_vals <- p_vals[p_vals >= 0]   # precip cannot be negative

if (length(p_vals) == 0) {
  stop("PREC: No valid (non-negative, non-NA) precipitation values found.\n",
       "File may be entirely corrupt or masked: ", prec_diag_file)
}

p_max  <- round(max(p_vals), 1)
message("pr July 2000 raw max: ", p_max)

# Thresholds account for Western Ghats extreme rainfall:
# Raw mm×100: July peak monsoon ~20,000-50,000 (200-500mm × 100)
# Already mm: July peak monsoon W.Ghats can legitimately reach 1500-2000mm
# Gap between 2000 and 5000 is ambiguous — warn but assume mm (safer default)
if (p_max > 5000) {
  PREC_SCALING_NEEDED <- TRUE
  message("PREC: Raw mm×100 detected — r/100 will be applied.")
} else if (p_max <= 2000) {
  PREC_SCALING_NEEDED <- FALSE
  message("PREC: Already in mm — no scaling needed.")
} else if (p_max > 2000 & p_max <= 5000) {
  PREC_SCALING_NEEDED <- FALSE
  message("PREC: WARNING — ambiguous max value ", p_max, " mm.")
  message("      Assuming already in mm (plausible for W.Ghats peak monsoon).")
  message("      If BIO12 sanity check fails, this assumption may be wrong.")
} else {
  stop("PREC: Negative max or zero values detected: ", p_max,
       "\nCheck CHELSA pr files for corruption.")
}
rm(p_vals); gc()


# ---------------------------------------------------------------------------
# 3. SCALING FUNCTIONS — use independent flags for temp and precip
# ---------------------------------------------------------------------------

scale_temp <- function(r) {
  if      (TEMP_SCALE_TYPE == "K10") (r / 10) - 273.15
  else if (TEMP_SCALE_TYPE == "K")   r - 273.15
  else                               r
}
scale_prec <- function(r) if (PREC_SCALING_NEEDED) r / 100 else r


# ---------------------------------------------------------------------------
# 4. LAND MASK — hard stop if missing, resampled ONCE
# ---------------------------------------------------------------------------

message("\n", strrep("=", 60))
message("LOADING LAND MASK")
message(strrep("=", 60))

land_mask_file <- file.path(gee_glc_dir, "GLC_2000.tif")

if (!file.exists(land_mask_file)) {
  stop("LAND MASK NOT FOUND:\n", land_mask_file,
       "\nHard stop — ocean pixels cannot be included for a coastal species.",
       "\nCheck GEE_Exports path and ensure GLC_2000.tif is present.")
}

chelsa_template <- rast(file.path(chelsa_root, "TRAINING_2000_2010",
                                   "tasmax", "CHELSA_tasmax_01_2000_AOI.tif"))
raw_mask <- rast(land_mask_file) > 0

message("GLC mask CRS  : ", crs(raw_mask,        describe = TRUE)$name)
message("CHELSA CRS    : ", crs(chelsa_template,  describe = TRUE)$name)

if (!terra::same.crs(raw_mask, chelsa_template)) {
  message("CRS mismatch — reprojecting with project() before resample.")
  raw_mask <- project(raw_mask, chelsa_template, method = "near")
} else {
  message("CRS match confirmed.")
}

mask_rs <- resample(raw_mask, chelsa_template, method = "near")
message("Land mask resampled to CHELSA grid (reused for all periods).")

rm(raw_mask, chelsa_template)
gc()


# ---------------------------------------------------------------------------
# 5. CORE FUNCTION: compute period climatology BIO1-19
#
#    Step 1 — For each month (1-12), load that month across all years in
#             the period and compute pixel-wise mean → 12 monthly means.
#    Step 2 — Feed 12 monthly means into dismo::biovars() once.
#
#    Hard stop if any expected file is missing.
#    Hard stop if year count doesn't match expected (partial period).
# ---------------------------------------------------------------------------

compute_period_bioclim <- function(period) {

  message("\n", strrep("=", 60))
  message("Period: ", period$label)
  message("Years : ", min(period$years), "-", max(period$years),
          " (", length(period$years), " years)")
  message(strrep("=", 60))

  # HARD STOP: verify expected year count
  if (length(period$years) != period$expected_yrs) {
    stop("Year count mismatch for ", period$label,
         ": expected ", period$expected_yrs,
         ", got ", length(period$years))
  }

  # HARD STOP: verify all files exist before starting
  message("Checking files...")
  for (var in c("tasmin", "tasmax", "pr")) {
    for (yr in period$years) {
      for (mo in MONTHS) {
        f <- file.path(chelsa_root, period$folder, var,
                       paste0("CHELSA_", var, "_", mo, "_", yr, "_AOI.tif"))
        if (!file.exists(f)) {
          stop("MISSING: ", f,
               "\nRe-run CHELSA_download_", var, ".R")
        }
      }
    }
  }
  message("All ", length(period$years) * 12 * 3,
          " input files confirmed present.")

  # Step 1: monthly period means — one variable at a time to manage RAM
  get_monthly_means <- function(var, scale_fn) {
    monthly_means <- vector("list", 12)
    for (mo_idx in seq_along(MONTHS)) {
      mo    <- MONTHS[mo_idx]
      files <- file.path(chelsa_root, period$folder, var,
                         paste0("CHELSA_", var, "_", mo, "_",
                                period$years, "_AOI.tif"))
      stk   <- scale_fn(rast(files))
      stk   <- mask(stk, mask_rs, maskvalues = 0)
      monthly_means[[mo_idx]] <- app(stk, mean, na.rm = TRUE)
      rm(stk); gc()
      message("  ", var, " month ", mo, " (", mo_idx, "/12) done")
    }
    rast(monthly_means)
  }

  message("\nComputing monthly means for tasmin...")
  tasmin_means <- get_monthly_means("tasmin", scale_temp)
  message("Computing monthly means for tasmax...")
  tasmax_means <- get_monthly_means("tasmax", scale_temp)
  message("Computing monthly means for pr...")
  pr_means     <- get_monthly_means("pr",     scale_prec)

  # Step 2: compute BIO1-19 from period monthly means
  message("\nRunning dismo::biovars()...")

  tasmin_b <- raster::brick(tasmin_means)
  tasmax_b <- raster::brick(tasmax_means)
  pr_b     <- raster::brick(pr_means)
  rm(tasmin_means, tasmax_means, pr_means); gc()

  bio_brick <- dismo::biovars(prec = pr_b, tmin = tasmin_b, tmax = tasmax_b)
  rm(tasmin_b, tasmax_b, pr_b); gc()

  bio <- rast(bio_brick)
  names(bio) <- BIO_NAMES
  rm(bio_brick); gc()

  return(bio)
}


# ---------------------------------------------------------------------------
# 6. RUN BOTH PERIODS
# ---------------------------------------------------------------------------

for (period in periods) {

  if (file.exists(period$out_file)) {
    # Validate existing file: 19 bands AND ecological ranges are plausible
    # A corrupt file may have 19 bands but wrong values — range check catches it
    check <- tryCatch({
      r_check    <- rast(period$out_file)
      n_bands    <- nlyr(r_check)
      bio1_range <- minmax(r_check[["BIO1"]])
      bio12_range<- minmax(r_check[["BIO12"]])
      rm(r_check)
      n_bands == 19 &&
        bio1_range[1]  > 5   && bio1_range[2]  < 40   &&   # mean temp °C (>5 accounts for W.Ghats high elevation)
        bio12_range[1] >= 0  && bio12_range[2] < 8000       # annual precip mm (<8000 accounts for W.Ghats high rainfall)
    }, error = function(e) FALSE)

    if (check) {
      message("\nSkipping (exists, 19 bands + ranges confirmed): ",
              basename(period$out_file))
      next
    } else {
      message("\nExisting file invalid or out of range — recomputing: ",
              basename(period$out_file))
    }
  }

  bio <- compute_period_bioclim(period)

  writeRaster(bio, period$out_file, overwrite = TRUE,
              datatype = "FLT4S",           # explicit float32 — correct for BIO vars
              gdal = c("COMPRESS=DEFLATE"))
  message("Saved: ", basename(period$out_file))

  rm(bio); gc()

  # Clear both temp file engines
  terra::tmpFiles(current = FALSE, orphan = TRUE, old = TRUE, remove = TRUE)
  raster::removeTmpFiles(h = 0)
}


# ---------------------------------------------------------------------------
# 7. SANITY CHECK
#    BIO1  (mean annual temp)  South India: expect 5-40 °C
#                              (W.Ghats high elevation ~12-14°C, coast ~28-30°C)
#    BIO12 (annual precip)     South India: expect 0-8000 mm
#                              (W.Ghats windward slopes can exceed 5000mm)
# ---------------------------------------------------------------------------

message("\n", strrep("=", 60))
message("SANITY CHECK")
message(strrep("=", 60))

for (period in periods) {
  if (!file.exists(period$out_file)) {
    message(basename(period$out_file), " — NOT FOUND")
    next
  }
  bio         <- rast(period$out_file)
  bio1_range  <- round(minmax(bio[["BIO1"]]),  2)
  bio12_range <- round(minmax(bio[["BIO12"]]), 1)
  bio1_flag   <- ifelse(bio1_range[1]  > 5  & bio1_range[2]  < 40,
                        " OK", " !! CHECK SCALING")
  bio12_flag  <- ifelse(bio12_range[1] >= 0  & bio12_range[2] < 8000,
                        " OK", " !! CHECK SCALING")

  message("\n", period$label)
  message("  Bands : ", nlyr(bio),
          if (nlyr(bio) == 19) " (19 OK)" else " !! EXPECTED 19")
  message("  BIO1  (mean temp °C)    : ",
          bio1_range[1],  " to ", bio1_range[2],  bio1_flag)
  message("  BIO12 (annual precip mm): ",
          bio12_range[1], " to ", bio12_range[2], bio12_flag)
  rm(bio)
}

message("\n", strrep("=", 60))
message("COMPLETE.")
message("  ", file.path(out_root, "bioclim_train_2000_2010.tif"))
message("  ", file.path(out_root, "bioclim_test_2011_2021.tif"))
message("Next: 04_annual_extraction_VIF.R")
message(strrep("=", 60))
