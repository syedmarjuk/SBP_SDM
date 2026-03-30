# =============================================================================
# 03_VIF_predictor_selection.R
# -----------------------------------------------------------------------------
# Purpose : Remove collinear BIO predictors from BIO1-19 using Pearson
#           correlation and Variance Inflation Factor (VIF) analysis.
#
# SCOPE (updated):
#   VIF filtering is applied to CLIMATIC PREDICTORS ONLY (BIO1-19).
#   Non-climatic predictors (GLC, GWL, JRC, DistCoast, DistCoastalWet)
#   are NOT filtered here — all are carried forward to MaxEnt as-is.
#   GLC and GWL are categorical; the others are continuous but ecologically
#   distinct enough from BIO variables that multicollinearity with them
#   is not a concern for this study.
#
# METHOD:
#   1. Load train_thinned.csv (237 training points — no year needed)
#   2. Extract BIO1-19 values from the period-mean training stack
#   3. Pearson correlation pass: drop one from any pair with |r| > 0.7
#      (the variable with higher VIF in that pair is dropped)
#   4. VIF pass: iteratively drop variable with highest VIF until all < 10
#   5. Save surviving BIO names + non-climatic predictor list to .rds
#
# Inputs  :
#   data/processed/train_thinned.csv              (n=237, from Script 01)
#   data/processed/bioclim/bioclim_train_2000_2010.tif  (19 bands)
#
# Outputs :
#   data/processed/vif/
#     vif_selected_vars.rds        — named list (see structure below)
#     vif_correlation_matrix.csv   — full Pearson matrix (for reporting)
#     vif_results.csv              — final VIF scores (for reporting)
#     correlation_matrix.png       — heatmap visualisation
#     vif_report.txt               — human-readable summary
#
# Output RDS structure:
#   vars <- readRDS("vif_selected_vars.rds")
#   vars$bio_selected    — character vector: surviving BIO variable names
#   vars$nonclimatic     — character vector: all non-climatic predictor names
#   vars$categorical     — character vector: c("GLC", "GWL")
#   vars$all             — character vector: full predictor list for MaxEnt
#
# Author  : [Your name]
# Date    : 2026-03-28
# =============================================================================


# ---------------------------------------------------------------------------
# 0. PACKAGES
# ---------------------------------------------------------------------------

required_pkgs <- c("terra", "here", "usdm", "corrplot", "dplyr")
new_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)

library(terra)
library(here)
library(usdm)
library(corrplot)
library(dplyr)

options(scipen = 999)


# ---------------------------------------------------------------------------
# 1. PATHS
# ---------------------------------------------------------------------------

bioclim_train_file <- here("data", "processed", "bioclim",
                           "bioclim_train_2000_2010.tif")
train_file         <- here("data", "processed", "train_thinned.csv")

out_dir <- here("data", "processed", "vif")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ---------------------------------------------------------------------------
# 2. INPUT VALIDATION
# ---------------------------------------------------------------------------

for (f in c(bioclim_train_file, train_file)) {
  if (!file.exists(f)) stop("MISSING INPUT FILE: ", f)
}
message("Input files confirmed.")


# ---------------------------------------------------------------------------
# 3. LOAD OCCURRENCE POINTS AND BIO STACK
# ---------------------------------------------------------------------------

message("Loading training occurrences...")
train_occ <- read.csv(train_file)
message("  Points loaded: ", nrow(train_occ))

if (!all(c("decimalLongitude", "decimalLatitude") %in% names(train_occ))) {
  stop("train_thinned.csv must have columns: decimalLongitude, decimalLatitude")
}

message("Loading BIO1-19 stack...")
bio_stack <- rast(bioclim_train_file)
message("  Bands: ", nlyr(bio_stack), "  (expected 19)")

if (nlyr(bio_stack) != 19) {
  stop("Expected 19 bands in bioclim_train_2000_2010.tif, found: ", nlyr(bio_stack))
}

# Ensure band names are BIO1...BIO19
names(bio_stack) <- paste0("BIO", 1:19)


# ---------------------------------------------------------------------------
# 4. EXTRACT BIO1-19 AT ALL TRAINING POINTS
# ---------------------------------------------------------------------------

message("\nExtracting BIO1-19 at ", nrow(train_occ), " training points...")

pts <- vect(train_occ,
            geom = c("decimalLongitude", "decimalLatitude"),
            crs  = "EPSG:4326")

bio_vals <- terra::extract(bio_stack, pts, ID = FALSE)

n_complete <- sum(complete.cases(bio_vals))
message("  Complete rows: ", n_complete, " / ", nrow(bio_vals))

if (n_complete < nrow(bio_vals) * 0.9) {
  stop("More than 10% of training points returned NA BIO values. ",
       "Check that the BIO stack covers the full AOI.")
}

# Drop any NA rows before VIF
bio_vals <- bio_vals[complete.cases(bio_vals), , drop = FALSE]
message("  Using ", nrow(bio_vals), " points for VIF analysis.")


# ---------------------------------------------------------------------------
# 5. PEARSON CORRELATION PASS  (|r| > 0.7)
#    usdm::vifcor() computes pairwise Pearson correlations and removes
#    one variable from each highly correlated pair (the one with higher VIF).
# ---------------------------------------------------------------------------

message("\n", strrep("=", 55))
message("STEP 1: Pearson correlation filter (|r| > 0.7)")
message(strrep("=", 55))

# Save full correlation matrix for reporting
cor_matrix <- cor(bio_vals, use = "complete.obs", method = "pearson")

write.csv(
  round(cor_matrix, 3),
  file.path(out_dir, "vif_correlation_matrix.csv")
)

# Correlation heatmap
png(file.path(out_dir, "correlation_matrix.png"),
    width = 1400, height = 1400, res = 130)
corrplot(
  cor_matrix,
  method      = "color",
  type        = "upper",
  tl.cex      = 0.8,
  tl.col      = "black",
  addCoef.col = "black",
  number.cex  = 0.5,
  title       = "Pearson Correlation — BIO1-19 (training period mean)",
  mar         = c(0, 0, 2, 0)
)
dev.off()
message("  Correlation matrix saved: vif_correlation_matrix.csv + correlation_matrix.png")

# Run vifcor
cor_result     <- usdm::vifcor(bio_vals, th = 0.7)
vars_after_cor <- cor_result@results$Variables

message("  BIO variables after correlation filter: ", length(vars_after_cor))
message("  Remaining: ", paste(vars_after_cor, collapse = ", "))
message("  Dropped  : ",
        paste(setdiff(paste0("BIO", 1:19), vars_after_cor), collapse = ", "))


# ---------------------------------------------------------------------------
# 6. VIF PASS  (VIF > 10)
#    Iteratively drops the variable with the highest VIF until all < 10.
# ---------------------------------------------------------------------------

message("\n", strrep("=", 55))
message("STEP 2: VIF filter (threshold = 10)")
message(strrep("=", 55))

df_after_cor <- bio_vals[, vars_after_cor, drop = FALSE]
vif_result   <- usdm::vifstep(df_after_cor, th = 10)
bio_selected <- vif_result@results$Variables

message("  BIO variables after VIF filter: ", length(bio_selected))
message("  Remaining: ", paste(bio_selected, collapse = ", "))
message("  Dropped  : ",
        paste(setdiff(vars_after_cor, bio_selected), collapse = ", "))

write.csv(
  vif_result@results,
  file.path(out_dir, "vif_results.csv"),
  row.names = FALSE
)
message("  VIF scores saved: vif_results.csv")


# ---------------------------------------------------------------------------
# 7. DEFINE FINAL PREDICTOR SET
#    Non-climatic predictors are all kept — not subject to VIF here.
#    GLC and GWL are flagged as categorical for MaxEnt.
# ---------------------------------------------------------------------------

nonclimatic_all  <- c("JRC", "DistCoast", "DistCoastalWet", "GLC", "GWL")
nonclimatic_cont <- c("JRC", "DistCoast", "DistCoastalWet")
categorical      <- c("GLC", "GWL")

final_predictors <- list(
  bio_selected  = bio_selected,                                 # VIF-filtered BIOs
  nonclimatic   = nonclimatic_all,                              # all non-climatic kept
  categorical   = categorical,                                  # GLC + GWL = factors
  all           = c(bio_selected, nonclimatic_cont, categorical) # full MaxEnt list
)

message("\n", strrep("=", 55))
message("FINAL PREDICTOR SET FOR MAXENT")
message(strrep("=", 55))
message("  BIO (VIF-selected, ", length(bio_selected), "): ",
        paste(bio_selected, collapse = ", "))
message("  Non-climatic continuous (all kept): ",
        paste(nonclimatic_cont, collapse = ", "))
message("  Categorical (always kept): GLC, GWL")
message("  TOTAL: ", length(final_predictors$all), " predictors")


# ---------------------------------------------------------------------------
# 8. SAVE OUTPUTS
# ---------------------------------------------------------------------------

saveRDS(final_predictors, file.path(out_dir, "vif_selected_vars.rds"))
message("\nSaved: vif_selected_vars.rds")

report <- c(
  "VIF PREDICTOR SELECTION REPORT",
  paste0("Date    : ", Sys.time()),
  paste0("Script  : 04_VIF_predictor_selection.R"),
  strrep("=", 55),
  "",
  "SCOPE: VIF applied to BIO1-19 only.",
  "  Non-climatic predictors (JRC, DistCoast, DistCoastalWet, GLC, GWL)",
  "  are NOT filtered — all carried forward to MaxEnt.",
  "",
  paste0("Training points used : ", nrow(bio_vals)),
  paste0("Starting BIO vars    : 19 (BIO1 to BIO19)"),
  paste0("Pearson threshold    : |r| > 0.7"),
  paste0("VIF threshold        : VIF > 10"),
  "",
  paste0("After correlation filter (", length(vars_after_cor), " BIO vars):"),
  paste0("  Kept   : ", paste(vars_after_cor, collapse = ", ")),
  paste0("  Dropped: ", paste(setdiff(paste0("BIO", 1:19), vars_after_cor), collapse = ", ")),
  "",
  paste0("After VIF filter (", length(bio_selected), " BIO vars):"),
  paste0("  Kept   : ", paste(bio_selected, collapse = ", ")),
  paste0("  Dropped: ", paste(setdiff(vars_after_cor, bio_selected), collapse = ", ")),
  "",
  "FINAL PREDICTOR SET:",
  paste0("  BIO (selected)  : ", paste(bio_selected, collapse = ", ")),
  paste0("  Non-climatic    : JRC, DistCoast, DistCoastalWet (continuous)"),
  paste0("  Categorical     : GLC, GWL (treat as factors in MaxEnt)"),
  paste0("  TOTAL           : ", length(final_predictors$all), " predictors"),
  "",
  "RDS structure: readRDS('vif_selected_vars.rds')",
  "  $bio_selected  — VIF-surviving BIO names",
  "  $nonclimatic   — all non-climatic names (incl. GLC, GWL)",
  "  $categorical   — c('GLC', 'GWL')",
  "  $all           — full predictor list for MaxEnt"
)

writeLines(report, file.path(out_dir, "vif_report.txt"))
message("Saved: vif_report.txt")

message("\n", strrep("=", 55))
message("Script 04 complete.")
message("Next: 05_ENMeval_tuning.R")
message(strrep("=", 55))
