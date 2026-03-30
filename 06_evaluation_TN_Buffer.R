# =============================================================================
# 06_evaluation_TN_Buffer.R (VERSION 8c - 50KM ECOLOGICAL BUFFER)
# -----------------------------------------------------------------------------
# PURPOSE:
#   Evaluate the v8 model (8 variables, RM 3.0) using an Ecological Buffer.
#   The background is sampled from TN + a 50km outward buffer, while
#   maintaining the 5km exclusion zone around all pelican sightings.
#
# OUTPUTS (all saved to data/v8c/evaluation/):
#   evaluation_metrics.csv
#   evaluation_report.txt
#   roc_curve.png
#   jackknife_test_auc.csv
#   jackknife_test_auc.png
# =============================================================================

rm(list = ls()); gc()
library(terra); library(here); library(sf); library(dplyr)
library(ggplot2); library(maxnet); library(ecospat)

options(scipen = 999)
set.seed(42)

# ---------------------------------------------------------------------------
# 1. CONFIGURATION
# ---------------------------------------------------------------------------
gee_root       <- here("GEE_Exports")
shapefile_path <- here("FINAL_STUDY_AREAWGS84", "FINAL_STUDY_AREAWGS84",
                       "study_area_tn_coast_10km_wgs84.shp")
PURGE_VARS     <- c("BIO14", "BIO15", "BIO5", "DistCoast")

test_occ_file     <- here("data", "processed", "test_thinned_TN.csv")
bioclim_test_file <- here("data", "processed", "bioclim", "bioclim_test_2011_2021.tif")
vif_rds_file      <- here("data", "processed", "vif", "vif_selected_vars.rds")
bias_test_file    <- here("data", "processed", "background", "bias_surface_test_2011_2021.tif")
maxent_v8_dir     <- here("data", "v8", "maxent")
occ_swd_file      <- here("data", "v3", "enmeval", "occ_swd.csv")
bg_swd_train_file <- here("data", "v3", "enmeval", "bg_swd.csv")

out_dir <- here("data", "v8c", "evaluation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 2. LOAD TN MASK & CREATE 50KM ECOLOGICAL BUFFER
# ---------------------------------------------------------------------------
message("Creating 50km Ecological Buffer...")
study_area <- sf::st_read(shapefile_path, quiet = TRUE)
if (sf::st_crs(study_area)$epsg != 4326)
  study_area <- sf::st_transform(study_area, 4326)
study_area_v  <- terra::vect(study_area)
eco_buffer_50km <- terra::buffer(study_area_v, width = 50000)

# ---------------------------------------------------------------------------
# 3. FIT V8 MODEL (IN-MEMORY)
# ---------------------------------------------------------------------------
message("Fitting v8 RM 3.0 model from SWD CSVs...")
occ_v8 <- read.csv(file.path(maxent_v8_dir, "maxent_swd_occ.csv"), stringsAsFactors = FALSE)
bg_v8  <- read.csv(file.path(maxent_v8_dir, "maxent_swd_bg.csv"),  stringsAsFactors = FALSE)

vars_full <- readRDS(vif_rds_file)
all_vars  <- vars_full$all[!vars_full$all %in% PURGE_VARS]
cat_vars  <- vars_full$categorical[!vars_full$categorical %in% PURGE_VARS]
bio_vars  <- vars_full$bio_selected[!vars_full$bio_selected %in% PURGE_VARS]

# Establish authoritative factor levels from the v8 training data
cat_levels <- list()
for (cv in cat_vars) {
  lvls             <- sort(unique(c(as.integer(occ_v8[[cv]]), as.integer(bg_v8[[cv]]))))
  cat_levels[[cv]] <- lvls
  occ_v8[[cv]]    <- factor(occ_v8[[cv]], levels = lvls)
  bg_v8[[cv]]     <- factor(bg_v8[[cv]],  levels = lvls)
}

p_vec  <- c(rep(1, nrow(occ_v8)), rep(0, nrow(bg_v8)))
all_df <- rbind(occ_v8[, all_vars, drop = FALSE], bg_v8[, all_vars, drop = FALSE])

best_model <- maxnet::maxnet(
  p       = p_vec,
  data    = all_df,
  f       = maxnet::maxnet.formula(p_vec, all_df, classes = "lq"),
  regmult = 3.0
)

# ---------------------------------------------------------------------------
# 4. LOAD RASTERS & EXTRACT TEST DATA (TN ONLY)
# ---------------------------------------------------------------------------
message("Extracting TN-only test data...")
bio_test      <- terra::rast(bioclim_test_file)
names(bio_test) <- paste0("BIO", 1:19)
bias_test_r   <- terra::rast(bias_test_file)
test_template <- bio_test[[1]]

snap_jrc <- terra::resample(terra::rast(file.path(gee_root, "SDM_JRC_Yearly",           "JRC_2016.tif")),            test_template, method = "bilinear")
snap_dcw <- terra::resample(terra::rast(file.path(gee_root, "SDM_DistCoastalWet_Yearly", "DistCoastalWet_2016.tif")), test_template, method = "bilinear")
snap_glc <- terra::resample(terra::rast(file.path(gee_root, "SDM_GLC_Yearly",           "GLC_2016.tif")),            test_template, method = "near")
snap_gwl <- terra::resample(terra::rast(file.path(gee_root, "SDM_GWL_Yearly",           "GWL_2016.tif")),            test_template, method = "near")

test_occ    <- read.csv(test_occ_file, stringsAsFactors = FALSE)
pts_strip_v <- terra::vect(test_occ, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")

test_env <- data.frame(terra::extract(bio_test[[bio_vars]], pts_strip_v, ID = FALSE))
test_env$DistCoastalWet <- terra::extract(snap_dcw, pts_strip_v, ID = FALSE)[, 1]
test_env$JRC            <- terra::extract(snap_jrc, pts_strip_v, ID = FALSE)[, 1]
test_env$GLC            <- as.integer(terra::extract(snap_glc, pts_strip_v, ID = FALSE)[, 1])
test_env$GWL            <- as.integer(terra::extract(snap_gwl, pts_strip_v, ID = FALSE)[, 1])
test_swd <- test_env[complete.cases(test_env[, all_vars]), ]

# ---------------------------------------------------------------------------
# 5. ECOLOGICAL BUFFER + EXCLUSION ZONE BACKGROUND SAMPLING
# ---------------------------------------------------------------------------
message("\n--- ECOLOGICAL BUFFER PROTOCOL ---")
occ_swd_train <- read.csv(occ_swd_file,      stringsAsFactors = FALSE)
bg_swd_train  <- read.csv(bg_swd_train_file, stringsAsFactors = FALSE)

all_lon <- c(occ_swd_train$longitude, test_occ$decimalLongitude)
all_lat <- c(occ_swd_train$latitude,  test_occ$decimalLatitude)
pts_all_v      <- terra::vect(data.frame(lon = all_lon, lat = all_lat),
                               geom = c("lon", "lat"), crs = "EPSG:4326")
occ_buffer_5km <- terra::buffer(pts_all_v, width = 5000)

# Mask bias to 50km ecological zone, then punch out 5km exclusion holes
bias_eco_masked  <- terra::mask(bias_test_r, eco_buffer_50km)
bias_final_safe  <- terra::mask(bias_eco_masked, occ_buffer_5km, inverse = TRUE)

bg_pts   <- terra::spatSample(bias_final_safe, size = 10000,
                               method = "weights", xy = TRUE, na.rm = TRUE)
pts_bg_v <- terra::vect(bg_pts, geom = c("x", "y"), crs = "EPSG:4326")

bg_env <- data.frame(terra::extract(bio_test[[bio_vars]], pts_bg_v, ID = FALSE))
bg_env$DistCoastalWet <- terra::extract(snap_dcw, pts_bg_v, ID = FALSE)[, 1]
bg_env$JRC            <- terra::extract(snap_jrc, pts_bg_v, ID = FALSE)[, 1]
bg_env$GLC            <- as.integer(terra::extract(snap_glc, pts_bg_v, ID = FALSE)[, 1])
bg_env$GWL            <- as.integer(terra::extract(snap_gwl, pts_bg_v, ID = FALSE)[, 1])
bg_eval_swd <- bg_env[complete.cases(bg_env[, all_vars]), ]

# ---------------------------------------------------------------------------
# 6. ALIGN FACTOR LEVELS & PREDICT
# ---------------------------------------------------------------------------
for (cv in cat_vars) {
  lvls               <- cat_levels[[cv]]
  occ_swd_train[[cv]] <- factor(occ_swd_train[[cv]], levels = lvls)
  bg_eval_swd[[cv]]   <- factor(bg_eval_swd[[cv]],   levels = lvls)
  test_swd[[cv]]      <- factor(test_swd[[cv]],       levels = lvls)
}

test_pred  <- predict(best_model, newdata = test_swd[, all_vars],      type = "cloglog", clamp = TRUE)
bg_pred    <- predict(best_model, newdata = bg_eval_swd[, all_vars],   type = "cloglog", clamp = TRUE)
train_pred <- predict(best_model, newdata = occ_swd_train[, all_vars], type = "cloglog", clamp = TRUE)

# ---------------------------------------------------------------------------
# 7. COMPUTE METRICS
# ---------------------------------------------------------------------------
message("\nCalculating metrics...")

auc_val    <- mean(sapply(test_pred, function(p) mean(p > bg_pred, na.rm = TRUE)))

thresholds  <- seq(0, 1, by = 0.005)
tss_vals    <- sapply(thresholds, function(t) {
  sens <- mean(test_pred >= t, na.rm = TRUE)
  spec <- mean(bg_pred   <  t, na.rm = TRUE)
  sens + spec - 1
})
best_thresh <- thresholds[which.max(tss_vals)]
best_tss    <- max(tss_vals, na.rm = TRUE)

boyce_result <- ecospat::ecospat.boyce(
  fit = c(test_pred, bg_pred), obs = test_pred,
  nclass = 0, PEplot = FALSE)
cbi_val <- boyce_result$cor

thresh_10p <- quantile(train_pred, 0.10, na.rm = TRUE)
or_10p     <- mean(test_pred < thresh_10p, na.rm = TRUE)
thresh_mtp <- min(train_pred, na.rm = TRUE)
or_mtp     <- mean(test_pred < thresh_mtp, na.rm = TRUE)

message("\n=============================================================")
message("  v8c TN (50km ECOLOGICAL BUFFER) — EVALUATION RESULTS")
message("=============================================================")
message("  AUC    : ", round(auc_val,  4))
message("  TSS    : ", round(best_tss, 4))
message("  CBI    : ", round(cbi_val,  4))
message("  OR10p  : ", round(or_10p,   4))
message("  OR_MTP : ", round(or_mtp,   4))
message("=============================================================")

# ---------------------------------------------------------------------------
# 8. SAVE METRICS CSV
# ---------------------------------------------------------------------------
metrics_df <- data.frame(
  metric = c("AUC", "TSS", "CBI", "OR10p", "OR_MTP",
             "Threshold_TSS", "Threshold_10p", "Threshold_MTP",
             "n_test", "n_train", "n_background"),
  value  = c(round(auc_val, 4),    round(best_tss, 4),  round(cbi_val, 4),
             round(or_10p, 4),     round(or_mtp, 4),
             round(best_thresh, 4), round(thresh_10p, 4), round(thresh_mtp, 4),
             nrow(test_swd), nrow(occ_swd_train), nrow(bg_eval_swd))
)
write.csv(metrics_df, file.path(out_dir, "evaluation_metrics.csv"), row.names = FALSE)
message("Saved: evaluation_metrics.csv")

# ---------------------------------------------------------------------------
# 9. ROC CURVE
# ---------------------------------------------------------------------------
message("Plotting ROC curve...")

roc_df <- do.call(rbind, lapply(thresholds, function(t) {
  data.frame(FPR = mean(bg_pred   >= t, na.rm = TRUE),
             TPR = mean(test_pred >= t, na.rm = TRUE))
}))

p_roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  geom_line(colour = "#378ADD", linewidth = 1) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60") +
  annotate("text", x = 0.62, y = 0.18,
           label = paste0("AUC = ",  round(auc_val,  3),
                          "\nCBI = ", round(cbi_val,  3),
                          "\nTSS = ", round(best_tss, 3)),
           size = 4, colour = "#0C447C", hjust = 0, lineheight = 1.5) +
  labs(title    = "ROC Curve — v8c (50km Ecological Buffer, TN Only)",
       subtitle = "Independent test data 2011-2021 | Tamil Nadu AOI",
       x = "False Positive Rate", y = "True Positive Rate") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggplot2::ggsave(file.path(out_dir, "roc_curve.png"),
                p_roc, width = 7, height = 6, dpi = 200, bg = "white")
message("Saved: roc_curve.png")

# ---------------------------------------------------------------------------
# 10. TRUE JACKKNIFE OF TEST AUC
#
#    HOW THE BUG WAS FIXED:
#    The old version used a single geom_col() with aes(fill = Type) and
#    position = "identity". This caused ggplot to draw bars layered in
#    factor-level order, with each bar's LEFT EDGE at 0 but drawn on top
#    of each other in an unpredictable way — the blue "With only" bars
#    appeared to start at ~0.35 because the teal "Without" bar was drawn
#    first and the blue bar was partially hidden behind it.
#
#    FIX: Split into three separate data frames and draw THREE independent
#    geom_col() layers. Each layer covers ALL rows but only for its own
#    type. All bars genuinely start from x = 0. Z-order is explicit:
#      1. Teal  (Without)  — drawn first, goes to back
#      2. Blue  (Only)     — drawn second, overlaps teal
#      3. Red   (Full)     — drawn last, only visible on the " " row
#    This exactly replicates the MaxEnt Java jackknife visual.
# ---------------------------------------------------------------------------
message("\nRunning True Jackknife of Test AUC (this takes a while)...")

calc_auc <- function(pred_p, pred_bg) {
  mean(sapply(pred_p, function(p) mean(p > pred_bg, na.rm = TRUE)))
}

# Training matrix (factor levels already set above)
train_env <- all_df   # combined occ + bg used to fit best_model

fit_and_eval <- function(vars_subset, fc = "lq", rm = 3.0) {
  d_train <- train_env[, vars_subset, drop = FALSE]
  d_test  <- test_swd[,  vars_subset, drop = FALSE]
  d_bg    <- bg_eval_swd[, vars_subset, drop = FALSE]

  f <- tryCatch(maxnet::maxnet.formula(p_vec, d_train, classes = fc),
                error = function(e) NULL)
  if (is.null(f)) return(0.5)

  fit <- tryCatch(maxnet::maxnet(p = p_vec, data = d_train, f = f, regmult = rm),
                  error = function(e) NULL)
  if (is.null(fit)) return(0.5)

  p_t <- tryCatch(predict(fit, newdata = d_test, type = "cloglog", clamp = TRUE),
                  error = function(e) NULL)
  p_b <- tryCatch(predict(fit, newdata = d_bg,   type = "cloglog", clamp = TRUE),
                  error = function(e) NULL)
  if (is.null(p_t) || is.null(p_b)) return(0.5)

  calc_auc(p_t, p_b)
}

auc_full  <- auc_val
jack_rows <- list()

for (v in all_vars) {
  message("  Jackknifing: ", v)
  auc_only <- fit_and_eval(v)
  auc_wo   <- fit_and_eval(all_vars[all_vars != v])

  jack_rows[[length(jack_rows) + 1]] <-
    data.frame(Variable = v, Type = "With only variable", AUC = auc_only,
               stringsAsFactors = FALSE)
  jack_rows[[length(jack_rows) + 1]] <-
    data.frame(Variable = v, Type = "Without variable",   AUC = auc_wo,
               stringsAsFactors = FALSE)
}

results <- do.call(rbind, jack_rows)
results  <- rbind(results,
  data.frame(Variable = " ", Type = "With all variables",
             AUC = auc_full, stringsAsFactors = FALSE)
)

write.csv(results, file.path(out_dir, "jackknife_test_auc.csv"), row.names = FALSE)
message("Saved: jackknife_test_auc.csv")

# --- PLOT ---
var_order <- c(" ", rev(sort(all_vars)))

df_wo   <- results[results$Type == "Without variable",   ]
df_only <- results[results$Type == "With only variable", ]
df_full <- results[results$Type == "With all variables", ]

df_wo$Variable   <- factor(df_wo$Variable,   levels = var_order)
df_only$Variable <- factor(df_only$Variable, levels = var_order)
df_full$Variable <- factor(df_full$Variable, levels = var_order)

p_jack <- ggplot() +
  # Layer 1 — TEAL "Without variable" (drawn to back)
  geom_col(data = df_wo,
           aes(x = AUC, y = Variable),
           fill = "#00A0A0", width = 0.6) +
  # Layer 2 — BLUE "With only variable" (overlaps teal, both start from 0)
  geom_col(data = df_only,
           aes(x = AUC, y = Variable),
           fill = "#0000FF", width = 0.6) +
  # Layer 3 — RED "With all variables" (only the " " row)
  geom_col(data = df_full,
           aes(x = AUC, y = Variable),
           fill = "#FF0000", width = 0.6) +
  coord_cartesian(xlim = c(0.3, 1.0)) +
  labs(
    title = "Jackknife of Test AUC for Pelecanus_philippensis (v8c)",
    x     = "Test AUC (2011-2021)",
    y     = "Environmental Variable"
  ) +
  annotate("rect",
           xmin = c(0.88, 0.88, 0.88),
           xmax = c(0.91, 0.91, 0.91),
           ymin = c(1.6, 2.1, 2.6),
           ymax = c(1.9, 2.4, 2.9),
           fill = c("#FF0000", "#00A0A0", "#0000FF")) +
  annotate("text",
           x     = c(0.92, 0.92, 0.92),
           y     = c(1.75, 2.25, 2.75),
           label = c("With all variables", "Without variable", "With only variable"),
           hjust = 0, size = 3.0) +
  theme_bw(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", hjust = 0.5),
    panel.grid.major.y = element_line(colour = "grey80"),
    panel.grid.minor   = element_blank(),
    legend.position    = "none"
  )

ggplot2::ggsave(file.path(out_dir, "jackknife_test_auc.png"),
                p_jack, width = 9, height = 5, dpi = 200, bg = "white")
message("Saved: jackknife_test_auc.png")

# ---------------------------------------------------------------------------
# 11. EVALUATION REPORT
# ---------------------------------------------------------------------------
report <- c(
  "MODEL EVALUATION REPORT (v8c — 50km Ecological Buffer, TN Only)",
  paste0("Date        : ", Sys.time()),
  paste0("Script      : v8c_07_evaluation_TN_Buffer.R"),
  strrep("=", 55),
  "",
  paste0("PREDICTORS (", length(all_vars), "): ", paste(all_vars, collapse = ", ")),
  "FC = LQ   |   RM = 3.0",
  "",
  "BACKGROUND STRATEGY:",
  "  50km ecological buffer around TN study area",
  "  5km exclusion buffer around all known occurrences",
  "  10,000 pts sampled from masked KDE bias surface (2011-2021)",
  "",
  "SAMPLE SIZES:",
  paste0("  Test occurrences : ", nrow(test_swd)),
  paste0("  Train occurrences: ", nrow(occ_swd_train)),
  paste0("  Background       : ", nrow(bg_eval_swd)),
  "",
  "METRICS:",
  paste0("  AUC    : ", round(auc_val,   4)),
  paste0("  TSS    : ", round(best_tss,  4)),
  paste0("  CBI    : ", round(cbi_val,   4)),
  paste0("  OR10p  : ", round(or_10p,    4)),
  paste0("  OR_MTP : ", round(or_mtp,    4)),
  "",
  "OUTPUTS:",
  "  evaluation_metrics.csv",
  "  evaluation_report.txt",
  "  roc_curve.png",
  "  jackknife_test_auc.csv",
  "  jackknife_test_auc.png"
)
writeLines(report, file.path(out_dir, "evaluation_report.txt"))
message("Saved: evaluation_report.txt")
message("\n  v8c COMPLETE — all outputs in data/v8c/evaluation/")
