# ================================================================
# 01_preprocess_newdata.R
#
# Notes:
#   - This is NOT yet the final PVIR-formatting script.
#   - This script demonstrates common preprocessing cases:
#       1. ALS / LiDAR reference data
#       2. NFI / plot-level or polygon-level data
#       3. NFI / source-specific preprocessing before RawPlots()
#
#   - Plot2Map::RawPlots() is interactive. Users are asked to identify
#     the relevant columns themselves.
#   - Final harmonization into the PVIR structure will be handled later in:
#       02_PVIR_format.R
#   - New PVIR 7 datasets produced here will later be combined with the
#     PVIR 6 reference dataset.
# ================================================================


# ================================================================
# 0. Setup
# ================================================================

# rm(list = ls())

# Prefer opening the project .Rproj instead of using setwd().
# If needed, uncomment and edit:
# setwd("C:/Users/araza001/OneDrive - Wageningen University & Research/AGB4/data/AGBref_2026")

suppressPackageStartupMessages({
  library(Plot2Map)
  library(sf)
  library(dplyr)
  library(data.table)
})

dir.create("outputs/preprocessed_demo", recursive = TRUE, showWarnings = FALSE)

MAP_YEAR <- 2015
MAP_RESOLUTION <- 100


# ================================================================
# 1. Helper: raw-data inspection before RawPlots()
# ================================================================
# Use this BEFORE RawPlots().
#
# Raw datasets can have different source-specific column names.
# Since RawPlots() is interactive, this step only inspects the table
# and supports light source-specific preparation where needed.
# ================================================================

inspect_raw_dataset <- function(dat, dataset_name = "unknown dataset") {
  
  message("\n================================================")
  message("Raw dataset inspection: ", dataset_name)
  message("================================================")
  message("Rows: ", nrow(dat))
  message("Columns: ", ncol(dat))
  message("\nColumn names:")
  print(names(dat))
  
  message("\nFirst rows:")
  print(utils::head(dat))
  
  invisible(TRUE)
}


# ================================================================
# 2. Helper: check standardized output AFTER RawPlots()
# ================================================================
# Use this AFTER RawPlots(), not before.
#
# RawPlots() should return/prepare a dataset with the minimum fields
# needed by the next Plot2Map functions.
# ================================================================

check_after_rawplots <- function(dat, dataset_name = "unknown dataset") {
  
  required_cols <- c("PLOT_ID", "POINT_X", "POINT_Y", "AGB_T_HA", "SIZE_HA")
  
  missing_cols <- setdiff(required_cols, names(dat))
  
  if (length(missing_cols) > 0) {
    stop(
      dataset_name,
      ": RawPlots() output is still missing required fields: ",
      paste(missing_cols, collapse = ", "),
      "\nCheck the interactive column selection or add missing fields manually.",
      call. = FALSE
    )
  }
  
  if (any(is.na(dat$longitude)) || any(is.na(dat$latitude))) {
    stop(dataset_name, ": missing coordinates after RawPlots().", call. = FALSE)
  }
  
  if (any(dat$longitude < -180 | dat$longitude > 180, na.rm = TRUE)) {
    stop(dataset_name, ": invalid longitude values after RawPlots().", call. = FALSE)
  }
  
  if (any(dat$latitude < -90 | dat$latitude > 90, na.rm = TRUE)) {
    stop(dataset_name, ": invalid latitude values after RawPlots().", call. = FALSE)
  }
  
  if (any(dat$SIZE_HA <= 0, na.rm = TRUE)) {
    stop(dataset_name, ": invalid SIZE_HA values after RawPlots().", call. = FALSE)
  }
  
  if (any(dat$AGB_T_HA < 0, na.rm = TRUE)) {
    stop(dataset_name, ": negative AGB_T_HA values after RawPlots().", call. = FALSE)
  }
  
  message(dataset_name, ": RawPlots() output passed minimum checks.")
  
  invisible(TRUE)
}


# ================================================================
# 3. Helper: standard post-RawPlots Plot2Map chain
# ================================================================
# Assumption:
#   dat is already the output from RawPlots(), or already formatted in
#   a Plot2Map-compatible structure.
# ================================================================

run_plot2map_preprocessing <- function(dat,
                                       dataset_name,
                                       map_year = 2015,
                                       map_resolution = 100,
                                       apply_deforestation_filter = TRUE) {
  
  check_after_rawplots(dat, dataset_name)
  
  if (apply_deforestation_filter) {
    dat_def <- Deforested(dat, map_year = map_year)
    dat_use <- dat_def$non_deforested_plots
  } else {
    dat_use <- dat
  }
  
  dat_biome <- BiomePair(dat_use)
  
  dat_uncertainty <- calculateTotalUncertainty(
    dat_biome,
    map_year = map_year,
    map_resolution = map_resolution
  )
  
  return(dat_uncertainty$data)
}


# ================================================================
# 4. Example 1: ALS / LiDAR reference data
# ================================================================
# Case:
#   Data are already accessible using Plot2Map::RefLidar().
#   Minimal preprocessing is needed.
#
# Example:
#   ALS_NL
#
# Notes:
#   - This represents a post-reference product type.
#   - We only standardize fields needed for uncertainty calculation.
# ================================================================

message("\nRunning Example 1: ALS / LiDAR reference data")

als_demo <- RefLidar("ALS_NL") #
als_demo <- BiomePair(als_demo)

# For 100 m resolution, 100 m x 100 m = 1 ha
als_demo$SIZE_HA <- round( (terra::res(terra::rast(list.files("ALS_NL", pattern = "\\.tif$", full.names = TRUE)[1]))) / 100, 0)

# Rename AGB to AGB_T_HA if needed
if ("AGB" %in% names(als_demo) && !"AGB_T_HA" %in% names(als_demo)) {
  names(als_demo)[names(als_demo) == "AGB"] <- "AGB_T_HA"
}

als_out <- calculateTotalUncertainty(
  als_demo,
  map_year = MAP_YEAR,
  map_resolution = MAP_RESOLUTION
)

als_code <- "EU_NLD"
als_type <- "LIDAR"

write.csv(
  als_out$data,
  file.path(
    "outputs/preprocessed_demo",
    paste0("ValidationData_DEMO_", als_code, "_", als_type, "_", MAP_YEAR, ".csv")
  ),
  row.names = FALSE
)

message("Example 1 completed: outputs/preprocessed_demo/ValidationData_DEMO_ALS_NL_2015.csv")


# ================================================================
# 5. Example 2: NFI / plot-level or polygon-level data
# ================================================================
# Case:
#   Dataset is plot-level or polygon-level.
#   Coordinates and biomass are mostly ready, but RawPlots() is used
#   to standardize the data interactively.
#
# Example:
#   Mongolia_NFI_SSU.gpkg
#
# Notes:
#   - If geometry is polygon, centroid is used as plot location.
#   - RawPlots() will ask the user to identify the relevant columns.
# ================================================================

message("\nRunning Example 2: NFI / plot-level or polygon-level data")

nfi_file <- "Mongolia_NFI_SSU.gpkg"

nfi_sf <- sf::st_read(nfi_file)

# If the file has polygon geometry, use centroid as plot location.
# If it has point geometry, centroid still returns the point location.
nfi_sf$longitude <- sf::st_coordinates(sf::st_centroid(nfi_sf))[, 1]
nfi_sf$latitude  <- sf::st_coordinates(sf::st_centroid(nfi_sf))[, 2]

nfi_df <- as.data.frame(nfi_sf)

inspect_raw_dataset(nfi_df, "DEMO_MONGOLIA_NFI") 

# Interactive Plot2Map formatting
nfi_plots <- RawPlots(nfi_df)
  nfi_plots <- nfi_plots[1:100,] #for demonstration! removing deforested plots take a while

# Standard Plot2Map chain after RawPlots()
nfi_out <- run_plot2map_preprocessing(
  dat = nfi_plots,
  dataset_name = "DEMO_MONGOLIA_NFI",
  map_year = MAP_YEAR,
  map_resolution = MAP_RESOLUTION,
  apply_deforestation_filter = TRUE
)
nfi_code <- "ASI_MONG"
nfi_type <- "NFI"

write.csv(
  nfi_out,
  file.path(
    "outputs/preprocessed_demo",
    paste0("ValidationData_DEMO_", nfi_code, "_", nfi_type, "_", MAP_YEAR, ".csv")
  ),
  row.names = FALSE
)
message("Example 2 completed: outputs/preprocessed_demo/ValidationData_DEMO_MONGOLIA_NFI_2015.csv")


# ================================================================
# 6. Example 3: Laos NFI with source-specific preprocessing
# ================================================================
# Case:
#   Dataset requires additional preparation before RawPlots().
#
# Example source-specific tasks:
#   - plot ID is absent or needs to be created
#   - coordinate fields may have source-specific names
#   - plot size may be known from report/publication, not from table
#   - measurement year may need to be added
#   - AGB is already available but under a different column name
#
# Important:
#   - This example assumes AGB is already available.
#   - We are NOT deriving AGB from DBH here.
#   - RawPlots() still handles the interactive column selection.
# ================================================================

message("\nRunning Example 3: Laos NFI with source-specific preprocessing")

laos_file <- "laos2019_nfi.csv"

laos_raw <- sf::st_read(laos_file)
laos_df <- as.data.frame(laos_raw)

inspect_raw_dataset(laos_df, "DEMO_LAOS_NFI_RAW")

# ------------------------------------------------
# Source-specific preparation only.
# Do not force all final Plot2Map names here.
# RawPlots() will ask the user to identify columns.
# ------------------------------------------------

laos_prepared <- laos_df %>%
  mutate(
    # Create temporary plot ID if absent.
    # During RawPlots(), choose this as plot ID if needed.
    plot_id_temp = if ("PLOT_ID" %in% names(.)) {
      as.character(PLOT_ID)
    } else {
      paste0("LAOS_", dplyr::row_number())
    },
    
    # Add known plot size from source documentation if absent.
    # Replace 0.1 with correct source-based plot area if needed.
    plot_size_temp_ha = if ("SIZE_HA" %in% names(.)) {
      as.numeric(SIZE_HA)
    } else {
      0.1
    },
    
    # Add source measurement year if useful.
    year_temp = if ("YEAR" %in% names(.)) {
      as.integer(YEAR)
    } else {
      2019
    }
  )

# Optional:
# Add coordinate helper columns only if obvious raw coordinate fields exist.
# Edit this part depending on the actual dataset.
#
# Example:
# if ("x" %in% names(laos_prepared) && !"longitude_temp" %in% names(laos_prepared)) {
#   laos_prepared$longitude_temp <- laos_prepared$x
# }
#
# if ("y" %in% names(laos_prepared) && !"latitude_temp" %in% names(laos_prepared)) {
#   laos_prepared$latitude_temp <- laos_prepared$y
# }

inspect_raw_dataset(laos_prepared, "DEMO_LAOS_NFI_PREPARED")

# Interactive Plot2Map formatting
laos_plots <- RawPlots(laos_prepared)

# Standard Plot2Map chain after RawPlots()
laos_out <- run_plot2map_preprocessing(
  dat = laos_plots,
  dataset_name = "DEMO_LAOS_NFI",
  map_year = MAP_YEAR,
  map_resolution = MAP_RESOLUTION,
  apply_deforestation_filter = TRUE
)
laos_code <- "ASI_LAO"
laos_type <- "NFI"

write.csv(
  laos_out,
  file.path(
    "outputs/preprocessed_demo",
    paste0("ValidationData_DEMO_", laos_code, "_", laos_type, "_", MAP_YEAR, ".csv")
  ),
  row.names = FALSE
)
message("Example 3 completed: outputs/preprocessed_demo/ValidationData_DEMO_LAOS_NFI_2015.csv")


# ================================================================
# 7. End
# ================================================================

message("\n================================================")
message("01_preprocess_newdata.R demo completed.")
message("Outputs written to: outputs/preprocessed_demo/")
message("Next script later: 02_PVIR_format.R")
message("================================================")