# =============================================================================
# 05_maxent_final_run.R (VERSION 8 — VARIABLE PURGE)
# -----------------------------------------------------------------------------
# PURPOSE:
#   Re-train MaxEnt after surgically removing four "traitor" predictors
#   identified from the v5 Jackknife analysis.
#
# VARIABLES PURGED (and why):
#   BIO15 — Precip Seasonality  : 0.00% permutation importance, delta=0.0000.
#            Removing it changes training gain by literally zero. Fully
#            absorbed by BIO12. Dead weight.
#   BIO14 — Precip Driest Month : 0.00% permutation importance, delta=0.0009.
#            Near-identical to BIO15 diagnosis. Redundant with BIO12.
#   BIO5  — Max Temp Warmest Mo.: 0.12% contribution, delta=0.0014. Already
#            captured by BIO3 (Isothermality). Creates temperature collinearity.
#   DistCoast — Dist. to coast  : 0.12% contribution, delta=0.0013. Signal
#            fully absorbed by DistCoastalWet (13.35% contribution). Redundant.
#
# V8 VARIABLE SET (8 predictors, down from 12):
#   BIO3, BIO4, BIO12, BIO18, JRC, DistCoastalWet, GLC, GWL
#
# NON-NEGOTIABLE SETTINGS (project guardrails):
#   Feature Classes : LQ  (Linear + Quadratic)
#   Regularization  : RM = 3.0
#   Training data   : data/v3/enmeval/  (locked baseline)
#   Jackknife       : TRUE  (always — needed for future purge rounds)
#
# Output: data/v8/maxent/
# =============================================================================

rm(list = ls()); gc()
library(here)
library(dplyr)
options(scipen = 999)

# ---------------------------------------------------------------------------
# 1. PATHS
# ---------------------------------------------------------------------------

maxent_jar   <- here("maxent", "maxent", "maxent.jar")

# LOCKED TRAINING PATHS
occ_swd_file <- here("data", "v3", "enmeval", "occ_swd.csv")
bg_swd_file  <- here("data", "v3", "enmeval", "bg_swd.csv")
vif_rds_file <- here("data", "processed", "vif", "vif_selected_vars.rds")

out_dir <- here("data", "v8", "maxent")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 2. SETTINGS
# ---------------------------------------------------------------------------

SPECIES_NAME <- "Pelecanus_philippensis"
best_fc      <- "LQ"
best_rm      <- 3.0

# V8 PURGE LIST — variables confirmed dead weight by Jackknife
PURGE_VARS <- c("BIO14", "BIO15", "BIO5", "DistCoast")

# ---------------------------------------------------------------------------
# 3. LOAD & FILTER VARIABLES
# ---------------------------------------------------------------------------

vars_full <- readRDS(vif_rds_file)

# Remove purged variables from all relevant lists
all_vars <- vars_full$all[!vars_full$all %in% PURGE_VARS]
cat_vars <- vars_full$categorical[!vars_full$categorical %in% PURGE_VARS]
bio_vars <- vars_full$bio_selected[!vars_full$bio_selected %in% PURGE_VARS]

message("Original predictors  : ", paste(vars_full$all, collapse = ", "))
message("Purged               : ", paste(PURGE_VARS, collapse = ", "))
message("V8 predictor set     : ", paste(all_vars, collapse = ", "))
message("V8 categorical vars  : ", paste(cat_vars, collapse = ", "))

# ---------------------------------------------------------------------------
# 4. LOAD TRAINING DATA
# ---------------------------------------------------------------------------

occ_swd <- read.csv(occ_swd_file, stringsAsFactors = FALSE)
bg_swd  <- read.csv(bg_swd_file,  stringsAsFactors = FALSE)

message("\nTraining occurrences : ", nrow(occ_swd))
message("Training background  : ", nrow(bg_swd))

# Standardise coordinate column names
names(occ_swd)[names(occ_swd) %in% c("x", "longitude")] <- "longitude"
names(occ_swd)[names(occ_swd) %in% c("y", "latitude")]  <- "latitude"
names(bg_swd)[names(bg_swd)   %in% c("x", "longitude")] <- "longitude"
names(bg_swd)[names(bg_swd)   %in% c("y", "latitude")]  <- "latitude"

# ---------------------------------------------------------------------------
# 5. BUILD SWD TABLES (V8 VARIABLE SET ONLY)
# ---------------------------------------------------------------------------

occ_jar <- data.frame(
  species   = SPECIES_NAME,
  longitude = occ_swd$longitude,
  latitude  = occ_swd$latitude,
  occ_swd[, all_vars, drop = FALSE]
)

bg_jar <- data.frame(
  species   = "background",
  longitude = bg_swd$longitude,
  latitude  = bg_swd$latitude,
  bg_swd[, all_vars, drop = FALSE]
)

for (cv in cat_vars) {
  if (cv %in% names(occ_jar)) occ_jar[[cv]] <- as.integer(occ_jar[[cv]])
  if (cv %in% names(bg_jar))  bg_jar[[cv]]  <- as.integer(bg_jar[[cv]])
}

occ_jar_path <- file.path(out_dir, "maxent_swd_occ.csv")
bg_jar_path  <- file.path(out_dir, "maxent_swd_bg.csv")
write.csv(occ_jar, occ_jar_path, row.names = FALSE)
write.csv(bg_jar,  bg_jar_path,  row.names = FALSE)
message("SWD CSVs written to  : ", out_dir)

# ---------------------------------------------------------------------------
# 6. MAXENT JAR ARGUMENTS
# ---------------------------------------------------------------------------

fc_flags  <- c("linear=true", "quadratic=true",
                "hinge=false", "product=false", "threshold=false")
cat_flags <- paste0("togglelayertype=", cat_vars)

out_dir_java <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)

args <- c(
  "-Xmx2048m", "-jar", normalizePath(maxent_jar, winslash = "/"),
  paste0("samplesfile=",         normalizePath(occ_jar_path, winslash = "/")),
  paste0("environmentallayers=", normalizePath(bg_jar_path,  winslash = "/")),
  paste0("outputdirectory=",     out_dir_java),
  paste0("betamultiplier=",      best_rm),
  fc_flags, cat_flags,
  "replicates=1",
  "randomtestpoints=0",
  "writebackgroundpredictions=true",
  "writeclampgrid=false",
  "writemess=false",
  "plots=true",
  "appendtoresultsfile=false",
  "askoverwrite=false",
  "autorun=true",
  "jackknife=true"
)

# ---------------------------------------------------------------------------
# 7. RUN
# ---------------------------------------------------------------------------

message("\nLaunching MaxEnt JAR for v8 (variable purge)...")
message("  FC       : ", best_fc)
message("  RM       : ", best_rm)
message("  Vars (", length(all_vars), ")  : ", paste(all_vars, collapse = ", "))
message("  Out      : ", out_dir)

system2(
  "java",
  args   = args,
  wait   = TRUE,
  stdout = file.path(out_dir, "maxent_stdout.txt"),
  stderr = file.path(out_dir, "maxent_stderr.txt")
)

message("\nSCRIPT 06 COMPLETE FOR v8")
message("  Next: run v8_07_evaluation.R")
