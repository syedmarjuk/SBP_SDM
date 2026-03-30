# =============================================================================
# 01b_train_test_with_year.R
# -----------------------------------------------------------------------------
# Purpose : Attach year column to BOTH thinned occurrence CSVs by joining
#           back to the raw GBIF data — WITHOUT re-thinning.
#
# Strategy: The thinned CSVs (train_thinned.csv, test_thinned.csv) are fixed
#           outputs from Script 01. This script does NOT touch spThin.
#           It simply looks up each thinned point's coordinates in the raw
#           data and retrieves the year, preserving the exact thinned sets.
#
# Inputs  :
#   GBIF_Occurrence_2000-2021/gbif_2000_2021.csv   (raw, tab-separated)
#   data/processed/train_thinned.csv               (n=237, no year)
#   data/processed/test_thinned.csv                (n=872, no year)
#
# Outputs :
#   data/processed/train_thinned_with_year.csv     (n=237, WITH year)
#   data/processed/test_thinned_with_year.csv      (n=872, WITH year)
#
# Note on multi-year coordinates:
#   If the same coordinate was observed in multiple years in the raw data,
#   this script assigns the EARLIEST year. This is deterministic and
#   consistent for both datasets.
#
# Author  : [Your name]
# Date    : 2026-03-28
# =============================================================================


# ---------------------------------------------------------------------------
# 0. PACKAGES
# ---------------------------------------------------------------------------

required_pkgs <- c("dplyr", "here", "readr", "stringr")
new_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)

library(dplyr)
library(here)
library(readr)
library(stringr)


# ---------------------------------------------------------------------------
# 1. CONSTANTS
# ---------------------------------------------------------------------------

SPECIES_NAME <- "Pelecanus philippensis"

AOI <- list(
  south = 7.510477,
  north = 13.964805,
  west  = 76.00262,
  east  = 80.999954
)

TRAIN_YEARS <- 2000:2010
TEST_YEARS  <- 2011:2021


# ---------------------------------------------------------------------------
# 2. LOAD RAW GBIF DATA
# ---------------------------------------------------------------------------

message("Loading raw GBIF data...")

raw <- readr::read_delim(
  here::here("GBIF_Occurrence_2000-2021", "gbif_2000_2021.csv"),
  delim     = "\t",
  col_types = readr::cols_only(
    species            = readr::col_character(),
    decimalLongitude   = readr::col_double(),
    decimalLatitude    = readr::col_double(),
    year               = readr::col_integer(),
    occurrenceStatus   = readr::col_character()
  ),
  quote = ""
)

message("  Raw records: ", nrow(raw))


# ---------------------------------------------------------------------------
# 3. BUILD YEAR-LOOKUP POOLS
#    One pool per split, filtered to the same AOI + PRESENT criteria used
#    in Script 01. Using arrange(year) + distinct ensures that when a
#    coordinate appears in multiple years, the EARLIEST year is kept.
# ---------------------------------------------------------------------------

message("Building year-lookup pools...")

cleaned_all <- raw |>
  dplyr::filter(
    !is.na(decimalLongitude),
    !is.na(decimalLatitude),
    !is.na(year),
    stringr::str_to_upper(occurrenceStatus) == "PRESENT",
    decimalLatitude  >= AOI$south,
    decimalLatitude  <= AOI$north,
    decimalLongitude >= AOI$west,
    decimalLongitude <= AOI$east
  ) |>
  dplyr::mutate(species = SPECIES_NAME) |>
  dplyr::select(species, decimalLongitude, decimalLatitude, year)

# Train pool: 2000-2010, one row per unique coordinate (earliest year kept)
train_pool <- cleaned_all |>
  dplyr::filter(year %in% TRAIN_YEARS) |>
  dplyr::arrange(year) |>
  dplyr::distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)

# Test pool: 2011-2021, one row per unique coordinate (earliest year kept)
test_pool <- cleaned_all |>
  dplyr::filter(year %in% TEST_YEARS) |>
  dplyr::arrange(year) |>
  dplyr::distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)

message("  Train pool unique coords: ", nrow(train_pool))
message("  Test  pool unique coords: ", nrow(test_pool))

rm(raw, cleaned_all); gc(verbose = FALSE)


# ---------------------------------------------------------------------------
# 4. HELPER — join year to a thinned CSV and validate
# ---------------------------------------------------------------------------

attach_year <- function(thinned_file, pool_df, split_name, expected_years) {
  
  thinned <- utils::read.csv(thinned_file, stringsAsFactors = FALSE)
  
  message("\n--- ", split_name, " ---")
  message("  Thinned points loaded: ", nrow(thinned))
  
  result <- thinned |>
    dplyr::left_join(
      pool_df,
      by = c("species", "decimalLongitude", "decimalLatitude"),
      relationship = "many-to-one"
    )
  
  # Validation 1: row count unchanged
  if (nrow(result) != nrow(thinned)) {
    stop(split_name, ": row count changed after join (",
         nrow(thinned), " -> ", nrow(result), "). Check for duplicates in pool.")
  }
  
  # Validation 2: no NA years
  na_count <- sum(is.na(result$year))
  if (na_count > 0) {
    stop(split_name, ": ", na_count, " thinned point(s) did not match any raw record.\n",
         "This likely means coordinate precision changed between thinning and raw data.\n",
         "Inspect: result[is.na(result$year), ]")
  }
  
  # Validation 3: no out-of-range years
  bad_years <- result$year[!result$year %in% expected_years]
  if (length(bad_years) > 0) {
    stop(split_name, ": unexpected years found: ",
         paste(sort(unique(bad_years)), collapse = ", "))
  }
  
  # Validation 4: no duplicate coordinates
  dup_check <- anyDuplicated(result[c("decimalLongitude", "decimalLatitude")])
  if (dup_check > 0) {
    stop(split_name, ": duplicate coordinates detected after year join.")
  }
  
  message("  Year range: ", min(result$year), " to ", max(result$year))
  message("  Year distribution:")
  print(table(result$year))
  
  result
}


# ---------------------------------------------------------------------------
# 5. ATTACH YEARS
# ---------------------------------------------------------------------------

train_with_year <- attach_year(
  thinned_file   = here::here("data", "processed", "train_thinned.csv"),
  pool_df        = train_pool,
  split_name     = "TRAINING (2000-2010)",
  expected_years = TRAIN_YEARS
)

test_with_year <- attach_year(
  thinned_file   = here::here("data", "processed", "test_thinned.csv"),
  pool_df        = test_pool,
  split_name     = "TESTING (2011-2021)",
  expected_years = TEST_YEARS
)


# ---------------------------------------------------------------------------
# 6. WRITE OUTPUTS
# ---------------------------------------------------------------------------

message("\nWriting outputs...")

utils::write.csv(
  train_with_year,
  here::here("data", "processed", "train_thinned_with_year.csv"),
  row.names = FALSE
)

utils::write.csv(
  test_with_year,
  here::here("data", "processed", "test_thinned_with_year.csv"),
  row.names = FALSE
)

message("  Saved: data/processed/train_thinned_with_year.csv  [n = ", nrow(train_with_year), "]")
message("  Saved: data/processed/test_thinned_with_year.csv   [n = ", nrow(test_with_year), "]")
message("\n01b complete. Both year files ready.")