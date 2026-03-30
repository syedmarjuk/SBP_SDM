# =============================================================================
# 04_ENMeval_tuning.R
# -----------------------------------------------------------------------------
# Purpose : Tune MaxEnt regularisation multiplier (RM) and feature classes (FC)
#           using ENMeval with the maxnet algorithm.
#
# SNAPSHOT DESIGN:
#   All training occurrence and background points use 2005 as the
#   representative mid-period snapshot for non-climatic predictors.
#   BIO variables use the period-mean stack (2000-2010).
#   DistCoast is static throughout.
#   Snapshot year: 2005 (mid-point of 2000-2010 training period)
#
# BIAS-CORRECTED BACKGROUND:
#   Background points are sampled proportional to observer effort using
#   a kernel density bias surface derived from 3,639 unique Aves observer
#   locations in the training AOI (from GBIF Aves India 2000-2010).
#   This corrects for spatial sampling bias in GBIF occurrence data,
#   directly improving CBI by ensuring background and presence points
#   share the same observer effort distribution.
#
# Inputs  :
#   data/processed/vif/vif_selected_vars.rds
#   data/processed/bioclim/bioclim_train_2000_2010.tif
#   data/processed/train_thinned_with_year.csv
#   data/processed/background/bias_surface.tif  (from Script 04d)
#   GEE_Exports/ (2005 snapshot rasters + DistCoast static)
#
# Outputs :
#   data/v3/enmeval/
#     occ_swd.csv, bg_swd.csv, enmeval_object.rds
#     enmeval_results.csv, enmeval_results_partitions.csv
#     best_tuning_settings.csv, tuning_or10p_cbi.png
#     tuning_aucdiff.png, enmeval_report.txt
# =============================================================================


# ---------------------------------------------------------------------------
# 0. PACKAGES
# ---------------------------------------------------------------------------

required_pkgs <- c("terra", "here", "ENMeval", "maxnet", "dplyr", "ggplot2")
new_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[,"Package"]]
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)

if (!"ecospat" %in% installed.packages()[,"Package"]) {
  message("Installing ecospat from R-Universe...")
  install.packages("ecospat",
                   repos = c("https://ropensci.r-universe.dev",
                             "https://cloud.r-project.org"),
                   dependencies = TRUE)
}

library(terra)
library(here)
library(ENMeval)
library(maxnet)
library(ecospat)
library(dplyr)
library(ggplot2)

options(scipen = 999)
set.seed(42)


# ---------------------------------------------------------------------------
# 1. CONFIGURATION
# ---------------------------------------------------------------------------

gee_root <- here("GEE_Exports")

bioclim_train_file <- here("data", "processed", "bioclim",
                           "bioclim_train_2000_2010.tif")
vif_rds_file       <- here("data", "processed", "vif", "vif_selected_vars.rds")
train_file         <- here("data", "processed", "train_thinned_with_year.csv")
bias_surface_file  <- here("data", "processed", "background", "bias_surface.tif")

out_dir <- here("data", "v3", "enmeval")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

SNAPSHOT_YEAR <- 2005
N_BG          <- 10000L

TUNE_FC <- c("LQ", "LQH", "H")
TUNE_RM <- c(0.5, 1, 1.5, 2, 3, 4)


# ---------------------------------------------------------------------------
# 2. INPUT VALIDATION
# ---------------------------------------------------------------------------

for (f in c(bioclim_train_file, vif_rds_file, train_file, bias_surface_file)) {
  if (!file.exists(f)) stop("Missing: ", f)
}

snapshot_files <- list(
  GLC            = file.path(gee_root, "SDM_GLC_Yearly",
                             paste0("GLC_", SNAPSHOT_YEAR, ".tif")),
  GWL            = file.path(gee_root, "SDM_GWL_Yearly",
                             paste0("GWL_", SNAPSHOT_YEAR, ".tif")),
  JRC            = file.path(gee_root, "SDM_JRC_Yearly",
                             paste0("JRC_", SNAPSHOT_YEAR, ".tif")),
  DistCoastalWet = file.path(gee_root, "SDM_DistCoastalWet_Yearly",
                             paste0("DistCoastalWet_", SNAPSHOT_YEAR, ".tif")),
  DistCoast      = file.path(gee_root, "SDM_Predictors_South_India",
                             "00_DistCoast_static.tif")
)
for (nm in names(snapshot_files)) {
  if (!file.exists(snapshot_files[[nm]])) {
    stop("Missing snapshot raster (", nm, "): ", snapshot_files[[nm]])
  }
}

message("All input files confirmed.")
message("Snapshot year     : ", SNAPSHOT_YEAR)
message("Background method : bias-corrected (kernel density from Aves observers)")


# ---------------------------------------------------------------------------
# 3. LOAD INPUTS
# ---------------------------------------------------------------------------

message("Loading inputs...")

vars     <- readRDS(vif_rds_file)
bio_vars <- vars$bio_selected
cat_vars <- vars$categorical
all_vars <- vars$all
message("  Predictors: ", paste(all_vars, collapse=", "))

bio_stack <- terra::rast(bioclim_train_file)
if (terra::nlyr(bio_stack) != 19) {
  stop("bioclim_train_2000_2010.tif must have exactly 19 bands, found: ",
       terra::nlyr(bio_stack))
}
names(bio_stack) <- paste0("BIO", 1:19)
template <- bio_stack[[1]]

# Verify all expected bio_vars exist in the named stack
missing_bio <- setdiff(bio_vars, names(bio_stack))
if (length(missing_bio) > 0) {
  stop("BIO variables not found in stack: ", paste(missing_bio, collapse=", "))
}

train_occ <- read.csv(train_file, stringsAsFactors = FALSE)
message("  Training occurrences: ", nrow(train_occ))

bias_r <- terra::rast(bias_surface_file)
message("  Bias surface loaded.")


# ---------------------------------------------------------------------------
# 4. LOAD AND ALIGN SNAPSHOT RASTERS
# ---------------------------------------------------------------------------

message("\nLoading and aligning ", SNAPSHOT_YEAR, " snapshot rasters...")

snap_jrc <- terra::resample(
  terra::rast(file.path(gee_root, "SDM_JRC_Yearly",
                        paste0("JRC_", SNAPSHOT_YEAR, ".tif"))),
  template, method = "bilinear")

snap_dcw <- terra::resample(
  terra::rast(file.path(gee_root, "SDM_DistCoastalWet_Yearly",
                        paste0("DistCoastalWet_", SNAPSHOT_YEAR, ".tif"))),
  template, method = "bilinear")

snap_glc <- terra::resample(
  terra::rast(file.path(gee_root, "SDM_GLC_Yearly",
                        paste0("GLC_", SNAPSHOT_YEAR, ".tif"))),
  template, method = "near")

snap_gwl <- terra::resample(
  terra::rast(file.path(gee_root, "SDM_GWL_Yearly",
                        paste0("GWL_", SNAPSHOT_YEAR, ".tif"))),
  template, method = "near")

snap_dc <- terra::resample(
  terra::rast(file.path(gee_root, "SDM_Predictors_South_India",
                        "00_DistCoast_static.tif")),
  template, method = "bilinear")

names(snap_jrc) <- "JRC"
names(snap_dcw) <- "DistCoastalWet"
names(snap_glc) <- "GLC"
names(snap_gwl) <- "GWL"
names(snap_dc)  <- "DistCoast"

message("  All snapshot rasters aligned.")


# ---------------------------------------------------------------------------
# 5. EXTRACT PREDICTORS AT OCCURRENCE POINTS
# ---------------------------------------------------------------------------

message("\n--- OCCURRENCE EXTRACTION ---")

pts_all <- terra::vect(train_occ,
                       geom = c("decimalLongitude", "decimalLatitude"),
                       crs  = "EPSG:4326")

bio_vals  <- terra::extract(bio_stack[[bio_vars]], pts_all, ID = FALSE)
dc_vals   <- terra::extract(snap_dc,  pts_all, ID = FALSE); names(dc_vals)  <- "DistCoast"
jrc_vals  <- terra::extract(snap_jrc, pts_all, ID = FALSE); names(jrc_vals) <- "JRC"
dcw_vals  <- terra::extract(snap_dcw, pts_all, ID = FALSE); names(dcw_vals) <- "DistCoastalWet"
glc_vals  <- terra::extract(snap_glc, pts_all, ID = FALSE); names(glc_vals) <- "GLC"
gwl_vals  <- terra::extract(snap_gwl, pts_all, ID = FALSE); names(gwl_vals) <- "GWL"

glc_vals$GLC <- as.integer(glc_vals$GLC)
gwl_vals$GWL <- as.integer(gwl_vals$GWL)

occ_swd <- cbind(
  train_occ[, c("decimalLongitude", "decimalLatitude")],
  bio_vals, dc_vals, jrc_vals, dcw_vals, glc_vals, gwl_vals
)
names(occ_swd)[1:2] <- c("x", "y")
occ_swd <- occ_swd[complete.cases(occ_swd[, all_vars]), ]
message("  Occurrence SWD complete rows: ", nrow(occ_swd))

write.csv(occ_swd, file.path(out_dir, "occ_swd.csv"), row.names = FALSE)


# ---------------------------------------------------------------------------
# 6. SAMPLE BIAS-CORRECTED BACKGROUND POINTS
#
#    Sample N_BG background points from the land surface using the bias
#    surface as sampling probability weights. Pixels with higher observer
#    density get proportionally more background points sampled from them.
#
#    This ensures background represents "where observers looked" rather
#    than "everywhere equally", directly correcting for sampling bias.
# ---------------------------------------------------------------------------

message("\n--- BIAS-CORRECTED BACKGROUND SAMPLING ---")
message("  Sampling ", N_BG, " background points weighted by observer effort...")

# bias_surface.tif is already aligned to BIO grid and masked to land
# in Script 04d — load and use directly, no re-processing needed
bias_aligned <- terra::resample(bias_r, template, method = "bilinear")
land_mask    <- terra::ifel(snap_glc > 0, 1, NA)
bias_ready   <- terra::mask(bias_aligned, land_mask)

# Sample proportional to bias surface values
bg_pts_v <- terra::spatSample(
  bias_ready,
  size      = N_BG * 3,
  method    = "weights",
  as.points = TRUE,
  values    = FALSE,
  na.rm     = TRUE,
  warn      = FALSE
)

message("  Candidate background points sampled: ", nrow(bg_pts_v))

# Extract predictors directly from SpatVector — no dataframe round-trip needed
bg_bio  <- terra::extract(bio_stack[[bio_vars]], bg_pts_v, ID = FALSE)
bg_dc   <- terra::extract(snap_dc,  bg_pts_v, ID = FALSE); names(bg_dc)  <- "DistCoast"
bg_jrc  <- terra::extract(snap_jrc, bg_pts_v, ID = FALSE); names(bg_jrc) <- "JRC"
bg_dcw  <- terra::extract(snap_dcw, bg_pts_v, ID = FALSE); names(bg_dcw) <- "DistCoastalWet"
bg_glc  <- terra::extract(snap_glc, bg_pts_v, ID = FALSE); names(bg_glc) <- "GLC"
bg_gwl  <- terra::extract(snap_gwl, bg_pts_v, ID = FALSE); names(bg_gwl) <- "GWL"

bg_glc$GLC <- as.integer(bg_glc$GLC)
bg_gwl$GWL <- as.integer(bg_gwl$GWL)

# Get coordinates from SpatVector directly
bg_coords        <- as.data.frame(terra::crds(bg_pts_v))
names(bg_coords) <- c("x", "y")

bg_swd <- cbind(bg_coords, bg_bio, bg_dc, bg_jrc, bg_dcw, bg_glc, bg_gwl)

bg_swd <- bg_swd[complete.cases(bg_swd[, all_vars]), ]
message("  Background complete rows: ", nrow(bg_swd))

if (nrow(bg_swd) < N_BG) {
  warning("Only ", nrow(bg_swd), " complete background rows — using all.")
} else {
  set.seed(42)
  bg_swd <- bg_swd[sample(seq_len(nrow(bg_swd)), N_BG), ]
}

write.csv(bg_swd, file.path(out_dir, "bg_swd.csv"), row.names = FALSE)
message("  Background SWD rows used: ", nrow(bg_swd))


# ---------------------------------------------------------------------------
# 7. RUN ENMevaluate
# ---------------------------------------------------------------------------

message("\n", strrep("=", 60))
message("RUNNING ENMevaluate (maxnet)")
message("  Snapshot year : ", SNAPSHOT_YEAR)
message("  Background    : bias-corrected (observer effort weighted)")
message("  FC            : ", paste(TUNE_FC, collapse=", "))
message("  RM            : ", paste(TUNE_RM, collapse=", "))
message("  Models to fit : ", length(TUNE_FC) * length(TUNE_RM))
message(strrep("=", 60))

e.mx <- ENMeval::ENMevaluate(
  occs               = occ_swd[, c("x", "y", all_vars)],
  bg                 = bg_swd[,  c("x", "y", all_vars)],
  algorithm          = "maxnet",
  tune.args          = list(fc = TUNE_FC, rm = TUNE_RM),
  partitions         = "block",
  partition.settings = list(orientation = "lat_lon"),
  categoricals       = cat_vars,
  raster.preds       = FALSE,
  taxon.name         = "Pelecanus_philippensis"
)

saveRDS(e.mx, file.path(out_dir, "enmeval_object.rds"))
message("Saved: enmeval_object.rds")

res      <- ENMeval::eval.results(e.mx)
res_part <- ENMeval::eval.results.partitions(e.mx)
write.csv(res,      file.path(out_dir, "enmeval_results.csv"),            row.names = FALSE)
write.csv(res_part, file.path(out_dir, "enmeval_results_partitions.csv"), row.names = FALSE)


# ---------------------------------------------------------------------------
# 8. SELECT BEST MODEL
# ---------------------------------------------------------------------------

message("\nSelecting best model...")

cbi_available  <- "cbi.val.avg"  %in% names(res) && any(!is.na(res$cbi.val.avg))
or10p_available <- "or.10p.avg" %in% names(res) && any(!is.na(res$or.10p.avg))

message("  CBI available   : ", cbi_available)
message("  or.10p available: ", or10p_available)

res_clean <- if (or10p_available) res[!is.na(res$or.10p.avg), ] else
                                  res[!is.na(res$auc.val.avg), ]

if (nrow(res_clean) == 0) stop("No models with valid metrics.")

if ("delta.AICc" %in% names(res_clean) && any(!is.na(res_clean$delta.AICc))) {
  candidate <- res_clean[!is.na(res_clean$delta.AICc) &
                           res_clean$delta.AICc <= 2, ]
  if (nrow(candidate) > 0) res_clean <- candidate
}

cbi_sort   <- if (cbi_available)   -res_clean$cbi.val.avg  else rep(0, nrow(res_clean))
or10p_sort <- if (or10p_available)  res_clean$or.10p.avg   else rep(0, nrow(res_clean))

best_row <- res_clean[order(
  or10p_sort, cbi_sort,
  ifelse(is.na(res_clean$auc.diff.avg), Inf, res_clean$auc.diff.avg),
  ifelse(is.na(res_clean$ncoef),        Inf, res_clean$ncoef),
  -ifelse(is.na(res_clean$rm),         -Inf, res_clean$rm)
)[1], , drop = FALSE]

write.csv(best_row, file.path(out_dir, "best_tuning_settings.csv"),
          row.names = FALSE)

message("  Best FC        : ", best_row$fc)
message("  Best RM        : ", best_row$rm)
message("  or.10p.avg     : ", round(best_row$or.10p.avg,  4))
message("  cbi.val.avg    : ", round(best_row$cbi.val.avg, 4))
message("  auc.val.avg    : ", round(best_row$auc.val.avg, 4))
message("  auc.diff.avg   : ", round(best_row$auc.diff.avg,4))
message("  ncoef          : ", best_row$ncoef)


# ---------------------------------------------------------------------------
# 9. PLOTS
# ---------------------------------------------------------------------------

p1 <- ENMeval::evalplot.stats(
  e=e.mx, stats=c("or.10p","cbi.val"),
  x.var="rm", color.var="fc", error.bars=FALSE)
ggplot2::ggsave(file.path(out_dir, "tuning_or10p_cbi.png"),
                p1, width=10, height=6, dpi=200)

p2 <- ENMeval::evalplot.stats(
  e=e.mx, stats=c("auc.diff","auc.val"),
  x.var="rm", color.var="fc", error.bars=FALSE)
ggplot2::ggsave(file.path(out_dir, "tuning_aucdiff.png"),
                p2, width=10, height=6, dpi=200)

message("Plots saved.")


# ---------------------------------------------------------------------------
# 10. REPORT
# ---------------------------------------------------------------------------

report <- c(
  "ENMeval TUNING REPORT",
  paste0("Date             : ", Sys.time()),
  strrep("=", 60),
  "",
  paste0("SNAPSHOT YEAR    : ", SNAPSHOT_YEAR),
  "BACKGROUND       : Bias-corrected (kernel density from 3,639 Aves",
  "                   observer locations, GBIF India 2000-2010)",
  paste0("Predictors       : ", paste(all_vars, collapse=", ")),
  paste0("Categorical      : ", paste(cat_vars, collapse=", ")),
  paste0("Occurrences      : ", nrow(occ_swd)),
  paste0("Background       : ", nrow(bg_swd)),
  paste0("FC tested        : ", paste(TUNE_FC, collapse=", ")),
  paste0("RM tested        : ", paste(TUNE_RM, collapse=", ")),
  paste0("Models fit       : ", length(TUNE_FC) * length(TUNE_RM)),
  "",
  "BEST MODEL",
  paste0("  FC             : ", best_row$fc),
  paste0("  RM             : ", best_row$rm),
  paste0("  or.10p.avg     : ", round(best_row$or.10p.avg,  4)),
  paste0("  cbi.val.avg    : ", round(best_row$cbi.val.avg, 4)),
  paste0("  auc.val.avg    : ", round(best_row$auc.val.avg, 4)),
  paste0("  auc.diff.avg   : ", round(best_row$auc.diff.avg,4)),
  paste0("  ncoef          : ", best_row$ncoef),
  "",
  paste0("Next: 06_maxent_final_run.R"),
  paste0("  Use FC=", best_row$fc, "  RM=", best_row$rm, " with maxent.jar")
)

writeLines(report, file.path(out_dir, "enmeval_report.txt"))
message("Saved: enmeval_report.txt")

message("\n", strrep("=", 60))
message("ENMeval TUNING COMPLETE")
message("Best: FC=", best_row$fc, "  RM=", best_row$rm)
message("Next: 06_maxent_final_run.R")
message(strrep("=", 60))
