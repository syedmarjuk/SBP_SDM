# =============================================================================
# 01_occurrence_processing.R
# -----------------------------------------------------------------------------
# Purpose : Clean, split, and spatially thin GBIF occurrence records of
#           Pelecanus philippensis for MaxEnt species distribution modelling.
#
# Inputs  : GBIF_Occurrence_2000-2021/gbif_2000_2021.csv (tab-separated GBIF download)
# Outputs : data/processed/train_thinned.csv   (MaxEnt-ready, 1 km thin)
#           data/processed/test_thinned.csv    (MaxEnt-ready, 5 km thin)
#           outputs/figures/occurrence_map.png
#
# Author  : [Your name]
# Date    : 2026-03-27
# R version tested : 4.3.x
# =============================================================================


# ---------------------------------------------------------------------------
# 0. PACKAGES
#    Install once, load every run.
#    spThin  : Aiello-Lammens et al. (2015) — the canonical SDM thinning method
#    here    : robust, project-relative file paths
#    readr/dplyr/stringr : lightweight data wrangling (avoids loading full tidyverse)
#    sf      : spatial operations (CRS handling)
# ---------------------------------------------------------------------------

required_pkgs <- c("here", "readr", "dplyr", "stringr", "spThin", "sf",
                   "ggplot2", "rnaturalearth", "rnaturalearthdata", "remotes")

new_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)

# rnaturalearthhires is NOT on CRAN — required by ne_states() for state boundaries
if (!"rnaturalearthhires" %in% installed.packages()[, "Package"]) {
  message("Installing rnaturalearthhires from GitHub...")
  remotes::install_github("ropensci/rnaturalearthhires")
}

library(here)
library(readr)
library(dplyr)
library(stringr)
library(spThin)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)


# ---------------------------------------------------------------------------
# 1. PROJECT PATHS  (using here:: — works regardless of working directory)
#    Expected folder structure inside your RStudio project:
#      project/
#        data/raw/          <- place your GBIF CSV here
#        data/processed/    <- thinned outputs written here
#        outputs/figures/   <- map written here
#        R/                 <- R scripts live here
# ---------------------------------------------------------------------------

dir.create(here("data", "processed"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("outputs", "figures"),showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# 2. STUDY CONSTANTS
#    Pinned from the exact min/max of gbif_2000_2021.csv — do not change.
# ---------------------------------------------------------------------------

AOI <- list(
  south = 7.510477,
  north = 13.964805,
  west  = 76.00262,
  east  = 80.999954
)

SPECIES_NAME  <- "Pelecanus philippensis"
TRAIN_YEARS   <- 2000:2010
TEST_YEARS    <- 2011:2021
THIN_DIST_KM_TRAIN <- 1   # 1 km  — fine grain, maximises training records
THIN_DIST_KM_TEST  <- 5   # 5 km  — coarser, ensures spatial independence for validation
THIN_REPS          <- 20  # 20 reps is sufficient; 100 is overkill for most datasets


# ---------------------------------------------------------------------------
# 3. LOAD RAW GBIF DATA
#    GBIF downloads are tab-separated (sep = "\t"), not comma-separated.
#    cols_only() reads just the 5 needed columns — saves RAM on wide files.
# ---------------------------------------------------------------------------

message("Loading raw GBIF data...")

raw <- read_delim(
  here("GBIF_Occurrence_2000-2021", "gbif_2000_2021.csv"),
  delim     = "\t",
  col_types = cols_only(                    # read ONLY needed columns
    species            = col_character(),
    decimalLongitude   = col_double(),
    decimalLatitude    = col_double(),
    year               = col_integer(),
    occurrenceStatus   = col_character()
  ),
  quote     = ""                             # GBIF rows sometimes contain quotes
)

message(sprintf("  Raw records loaded: %d", nrow(raw)))


# ---------------------------------------------------------------------------
# 4. CLEAN
#    a) Drop records missing coordinates or year
#    b) Clip to AOI bounding box (safety net — should all be inside already)
#    c) Keep only PRESENT occurrences
#    (Column types already set at read time via cols_only)
# ---------------------------------------------------------------------------

message("Cleaning data...")

cleaned <- raw |>
  filter(
    !is.na(decimalLongitude),
    !is.na(decimalLatitude),
    !is.na(year),
    str_to_upper(occurrenceStatus) == "PRESENT",
    decimalLatitude  >= AOI$south,
    decimalLatitude  <= AOI$north,
    decimalLongitude >= AOI$west,
    decimalLongitude <= AOI$east
  ) |>
  # Standardise species name (remove authorship if present)
  mutate(species = SPECIES_NAME) |>
  select(species, decimalLongitude, decimalLatitude, year)

message(sprintf("  Records after cleaning: %d", nrow(cleaned)))

# Save row counts before freeing objects — needed for summary report later
n_raw     <- nrow(raw)
n_cleaned <- nrow(cleaned)

# Free raw data from memory — no longer needed
rm(raw); gc(verbose = FALSE)


# ---------------------------------------------------------------------------
# 5. TEMPORAL SPLIT
#    Training  : 2000-2010 (used to build the MaxEnt model)
#    Testing   : 2011-2021 (independent temporal block for validation)
#    This is a TEMPORAL split, not random — critical for avoiding
#    temporal autocorrelation between train and test sets.
# ---------------------------------------------------------------------------

message("Splitting into training and testing sets...")

train_raw <- cleaned |> filter(year %in% TRAIN_YEARS)
test_raw  <- cleaned |> filter(year %in% TEST_YEARS)

# Free cleaned from memory — split into train/test now
rm(cleaned); gc(verbose = FALSE)

message(sprintf("  Training records (pre-thin): %d", nrow(train_raw)))
message(sprintf("  Testing  records (pre-thin): %d", nrow(test_raw)))

# Remove exact spatial duplicates BEFORE thinning — identical coordinates
# will be removed by spThin anyway, but they bloat the pairwise distance
# matrix and slow down every rep for no benefit.
train_raw <- train_raw |> distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)
test_raw  <- test_raw  |> distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)

message(sprintf("  Training unique locations:   %d", nrow(train_raw)))
message(sprintf("  Testing  unique locations:   %d", nrow(test_raw)))

# Guard against empty sets — fail early with a clear message
if (nrow(train_raw) == 0) stop("No training records remain after cleaning. Check filters/years.")
if (nrow(test_raw)  == 0) stop("No testing records remain after cleaning. Check filters/years.")

# Convert to base data.frame — spThin was written for base R and can
# choke on tibbles when matching column names internally.
train_raw <- as.data.frame(train_raw)
test_raw  <- as.data.frame(test_raw)


# ---------------------------------------------------------------------------
# 6. SPATIAL THINNING  —  spThin::thin()
#
#    WHY spThin?
#    The thin() function runs `reps` random orderings of the dataset and
#    greedily removes points that violate the minimum distance threshold.
#    Across all reps it selects the run that RETAINS THE MOST POINTS.
#    This is superior to simple grid-cell rarefaction because:
#      i)  it is not sensitive to arbitrary grid alignment
#      ii) it is reproducible across platforms
#      iii) it is the method cited in the SDM literature (Aiello-Lammens 2015)
#
#    DISTANCES:
#      Train 1 km  — preserve fine-scale habitat variation in training signal
#      Test  5 km  — ensure test points are spatially independent from each
#                    other within the test set (note: this does NOT enforce
#                    a minimum distance between train and test points)
# ---------------------------------------------------------------------------

message("Spatially thinning training data (1 km)...")

set.seed(42)  # reproducibility

thin_train <- thin(
  loc.data                  = train_raw,
  lat.col                   = "decimalLatitude",
  long.col                  = "decimalLongitude",
  spec.col                  = "species",
  thin.par                  = THIN_DIST_KM_TRAIN,
  reps                      = THIN_REPS,
  locs.thinned.list.return  = TRUE,
  write.files               = FALSE,
  verbose                   = FALSE
)

# spThin returns a list of length `reps`; pick the rep with most points
best_train_idx  <- which.max(sapply(thin_train, nrow))
train_thinned   <- thin_train[[best_train_idx]]

message(sprintf("  Training records after 1 km thinning: %d", nrow(train_thinned)))

# Free thinning list from memory
rm(thin_train); gc(verbose = FALSE)


message("Spatially thinning testing data (5 km)...")

thin_test <- thin(
  loc.data                  = test_raw,
  lat.col                   = "decimalLatitude",
  long.col                  = "decimalLongitude",
  spec.col                  = "species",
  thin.par                  = THIN_DIST_KM_TEST,
  reps                      = THIN_REPS,
  locs.thinned.list.return  = TRUE,
  write.files               = FALSE,
  verbose                   = FALSE
)

best_test_idx <- which.max(sapply(thin_test, nrow))
test_thinned  <- thin_test[[best_test_idx]]

message(sprintf("  Testing  records after 5 km thinning: %d", nrow(test_thinned)))

# Free thinning list from memory
rm(thin_test); gc(verbose = FALSE)


# ---------------------------------------------------------------------------
# 7. FORMAT FOR MAXENT
#    MaxEnt requires: species | longitude | latitude  (in that column order)
#    Column names must be exactly: species, decimalLongitude, decimalLatitude
# ---------------------------------------------------------------------------

format_for_maxent <- function(df, species_name) {
  # spThin typically returns columns named "Longitude" and "Latitude",
  # but this is not formally guaranteed. Detect them robustly.
  cols <- names(df)
  lon_col <- grep("^lon", cols, ignore.case = TRUE, value = TRUE)[1]
  lat_col <- grep("^lat", cols, ignore.case = TRUE, value = TRUE)[1]

  if (is.na(lon_col) || is.na(lat_col)) {
    stop("Cannot find longitude/latitude columns in spThin output. Found: ",
         paste(cols, collapse = ", "))
  }

  data.frame(
    species          = species_name,
    decimalLongitude = df[[lon_col]],
    decimalLatitude  = df[[lat_col]]
  )
}

train_out <- format_for_maxent(train_thinned, SPECIES_NAME)
test_out  <- format_for_maxent(test_thinned,  SPECIES_NAME)


# ---------------------------------------------------------------------------
# 8. WRITE OUTPUTS
# ---------------------------------------------------------------------------

message("Writing output CSVs...")

write_csv(train_out, here("data", "processed", "train_thinned.csv"))
write_csv(test_out,  here("data", "processed", "test_thinned.csv"))

message(sprintf("  Saved: data/processed/train_thinned.csv  [n = %d]", nrow(train_out)))
message(sprintf("  Saved: data/processed/test_thinned.csv   [n = %d]", nrow(test_out)))


# ---------------------------------------------------------------------------
# 9. SUMMARY REPORT
# ---------------------------------------------------------------------------

cat("\n", strrep("=", 55), "\n")
cat("  OCCURRENCE PROCESSING SUMMARY\n")
cat(strrep("=", 55), "\n")
cat(sprintf("  Species          : %s\n", SPECIES_NAME))
cat(sprintf("  AOI              : %.4f–%.4f°N, %.4f–%.4f°E\n",
            AOI$south, AOI$north, AOI$west, AOI$east))
cat(strrep("-", 55), "\n")
cat(sprintf("  %-30s %6d\n", "Raw records:",           n_raw))
cat(sprintf("  %-30s %6d\n", "After cleaning:",        n_cleaned))
cat(strrep("-", 55), "\n")
cat(sprintf("  %-30s %6d\n", "Train raw (2000-2010):", nrow(train_raw)))
cat(sprintf("  %-30s %6d\n",
            sprintf("Train thinned (%d km):", THIN_DIST_KM_TRAIN), nrow(train_out)))
cat(strrep("-", 55), "\n")
cat(sprintf("  %-30s %6d\n", "Test raw  (2011-2021):", nrow(test_raw)))
cat(sprintf("  %-30s %6d\n",
            sprintf("Test thinned  (%d km):", THIN_DIST_KM_TEST),  nrow(test_out)))
cat(strrep("=", 55), "\n\n")


# ---------------------------------------------------------------------------
# 10. DIAGNOSTIC MAP
#     Visualise both thinned datasets over South India.
#     Train = dark blue filled circles
#     Test  = red filled circles
# ---------------------------------------------------------------------------

message("Generating occurrence map...")

world <- ne_countries(scale = "small", returnclass = "sf")
india <- ne_states(country = "india", returnclass = "sf")

aoi_box <- st_as_sfc(
  st_bbox(c(xmin = AOI$west, xmax = AOI$east,
            ymin = AOI$south, ymax = AOI$north),
          crs = 4326)
)

p <- ggplot() +
  geom_sf(data = world, fill = "#e8e0d4", colour = "#b0a898", linewidth = 0.3) +
  geom_sf(data = india, fill = "#dbd3c4", colour = "#a09080", linewidth = 0.25) +
  geom_sf(data = aoi_box, fill = NA, colour = "#D85A30",
          linewidth = 0.8, linetype = "dashed") +
  geom_point(data = test_out,
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#C0392B", size = 1.2, alpha = 0.55, shape = 16) +
  geom_point(data = train_out,
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "#1A3A5C", size = 2, alpha = 0.85, shape = 16) +
  coord_sf(xlim = c(AOI$west - 1, AOI$east + 1),
           ylim = c(AOI$south - 0.5, AOI$north + 0.5),
           expand = FALSE) +
  labs(
    title    = expression(italic("Pelecanus philippensis") ~ "— Thinned Occurrence Records"),
    subtitle = sprintf("Train 2000–2010 (1 km thin, n = %d)  |  Test 2011–2021 (5 km thin, n = %d)",
                       nrow(train_out), nrow(test_out)),
    x = "Longitude (°E)", y = "Latitude (°N)",
    caption = "Source: GBIF (2000–2021). Dashed rectangle = study AOI."
  ) +
  annotate("point", x = -Inf, y = Inf, size = 0) +   # force axis
  annotate("text", x = AOI$west + 0.1, y = AOI$north - 0.15,
           label = "AOI", colour = "#D85A30", size = 3, hjust = 0) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12),
    plot.subtitle    = element_text(size = 9, colour = "grey40"),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
    legend.position  = "none"
  ) +
  # Manual legend via guide_legend workaround
  geom_point(aes(x = AOI$west + 0.2, y = AOI$south + 0.5),
             colour = "#1A3A5C", size = 2.5, shape = 16) +
  geom_point(aes(x = AOI$west + 0.2, y = AOI$south + 0.2),
             colour = "#C0392B", size = 2.5, shape = 16) +
  annotate("text", x = AOI$west + 0.45, y = AOI$south + 0.5,
           label = "Train (2000-2010)", size = 2.8, hjust = 0) +
  annotate("text", x = AOI$west + 0.45, y = AOI$south + 0.2,
           label = "Test (2011-2021)", size = 2.8, hjust = 0)

ggsave(
  here("outputs", "figures", "occurrence_map.png"),
  plot   = p,
  width  = 8,
  height = 9,
  dpi    = 200,
  bg     = "white"
)

message("  Saved: outputs/figures/occurrence_map.png")
message("\nDone. All outputs written successfully.")
