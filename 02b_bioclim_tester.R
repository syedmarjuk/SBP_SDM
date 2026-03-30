library(terra)
library(here)

# Paths to your shiny new climatology files
train_file <- here("data", "processed", "bioclim", "bioclim_train_2000_2010.tif")
test_file  <- here("data", "processed", "bioclim", "bioclim_test_2011_2021.tif")

files_to_check <- list("TRAINING (2000-2010)" = train_file, 
                       "TESTING (2011-2021)" = test_file)

for (name in names(files_to_check)) {
  f <- files_to_check[[name]]
  message("\n========================================")
  message("CHECKING: ", name)
  message("========================================")
  
  if (!file.exists(f)) {
    message("  [!] File not found. Are you sure it finished running?")
    next
  }
  
  r <- terra::rast(f)
  
  # 1. Band Check
  bands <- terra::nlyr(r)
  message("  Layers (Bands): ", bands, if(bands == 19) "  [PASS]" else "  [FAIL - Expected 19]")
  
  # 2. Temperature Check (BIO1)
  t_range <- round(terra::minmax(r[["BIO1"]]), 2)
  message("  BIO1 (Mean Temp)   : ", t_range[1], " to ", t_range[2], " °C")
  
  # 3. Precipitation Check (BIO12)
  p_range <- round(terra::minmax(r[["BIO12"]]), 1)
  message("  BIO12 (Annual Rain): ", p_range[1], " to ", p_range[2], " mm")
}