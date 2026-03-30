# =============================================================================
# Exact Point-by-Point Validation: Thinned vs Raw GBIF
# =============================================================================

required_pkgs <- c("dplyr", "here")
new_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(new_pkgs)) install.packages(new_pkgs, dependencies = TRUE)

if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(dplyr)
library(here)

# 1. DEFINE FILE PATHS (project-relative via here())
raw_file     <- here("GBIF_Occurrence_2000-2021", "gbif_2000_2021.csv")
test_file    <- here("data", "processed", "test_thinned_with_year.csv")
test_tn_file <- here("data", "processed", "test_thinned_with_year_TN.csv")
train_file   <- here("data", "processed", "train_thinned_with_year.csv")

# 2. LOAD THE DATA
message("Loading datasets for validation...")
raw_data <- read.delim(raw_file, quote = "") 
test     <- read.csv(test_file)
test_tn  <- read.csv(test_tn_file)
train    <- read.csv(train_file)

# 3. STANDARDIZE PRECISION (Crucial for exact matching in R)
# We round to 5 decimal places (~1 meter precision) to avoid floating point errors
raw_data <- raw_data %>%
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
  mutate(
    decimalLongitude = round(decimalLongitude, 5),
    decimalLatitude  = round(decimalLatitude, 5)
  )

# 4. VALIDATION FUNCTION
verify_dataset <- function(thinned_df, raw_df, dataset_name) {
  message("\n", strrep("-", 50))
  message("CHECKING: ", dataset_name)
  
  # Standardize thinned data precision
  thinned_df <- thinned_df %>%
    mutate(
      decimalLongitude = round(decimalLongitude, 5),
      decimalLatitude  = round(decimalLatitude, 5)
    )
  
  n_total <- nrow(thinned_df)
  
  # A. Check for Duplicates
  n_unique <- thinned_df %>% select(decimalLongitude, decimalLatitude) %>% distinct() %>% nrow()
  if (n_total == n_unique) {
    message("  [PASS] No duplicate coordinates detected.")
  } else {
    message("  [FAIL] Found ", n_total - n_unique, " duplicate coordinates!")
  }
  
  # B. Check Exact Match in Raw Data
  # semi_join keeps only rows in thinned_df that have a perfect match in raw_df
  matches <- semi_join(thinned_df, raw_df, by = c("decimalLongitude", "decimalLatitude", "year"))
  
  if (nrow(matches) == n_total) {
    message("  [PASS] All ", n_total, " points perfectly match a coordinate/year in the raw data.")
  } else {
    missing <- n_total - nrow(matches)
    message("  [FAIL] ", missing, " points DO NOT exist in the raw data with that specific year!")
  }
}

# 5. RUN THE TESTS
verify_dataset(train, raw_data, "Train Data (2000-2010)")
verify_dataset(test, raw_data, "Test Data - All (2011-2021)")
verify_dataset(test_tn, raw_data, "Test Data - TN Only (2011-2021)")

message("\nValidation complete.")