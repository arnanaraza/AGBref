# ================================================================
# 03_validation.R
# Run PVIR validation from combined PVIR6 + newly processed PVIR7 data
#
# Fully independent from local D:/treecover2010_v3.
# Uses gfcanalysis to download raw Hansen GFC treecover2000 + lossyear.
# No epoch-specific 30m treecover tiles are written.
#
# Epoch-specific forest fraction is computed on the fly using:
#   treecover2000, with lossyear <= map_year - 2000 set to non-forest.
# ================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  dplyr,
  stringr,
  plotrix,
  terra,
  raster,
  sp,
  sf,
  plyr,
  foreach,
  doParallel,
  parallel,
  ggplot2,
  Metrics,
  gfcanalysis
)

# ----------------------------
# 0. Project paths and settings
# ----------------------------

projectDir <- getwd()
dataDir    <- file.path(projectDir, "data")
rDir       <- file.path(projectDir, "R")
readyDir   <- file.path(projectDir, "outputs", "pvir_ready_demo")
resultsDir <- file.path(projectDir, "outputs", "pvir_validation")

dir.create(dataDir, recursive = TRUE, showWarnings = FALSE)
dir.create(resultsDir, recursive = TRUE, showWarnings = FALSE)

# Large external CCI dependencies
CCImaps <- "E:/CCI_v7/AGB_GTiff_100m"

agbTilesFolderBase <- Sys.getenv(
  "PVIR_AGB_TILES_100M",
  unset = CCImaps
)

# Raw GFC download folder only. No epoch treecover writing.
gfcDir <- file.path(dataDir, "GFC")
dir.create(gfcDir, recursive = TRUE, showWarnings = FALSE)

# Important:
# treeCoverFolder is kept for compatibility with old callAID/AGBInv arguments.
# It now points to raw GFC folder, not treecover2010_v3.
treeCoverFolder <- gfcDir

sdTilesFolderBase <- Sys.getenv(
  "PVIR_SD_TILES_100M",
  unset = NA_character_
)

pvirProductVersion <- "fv7.0"
pvirMapResolution  <- "100m"
forestTHs <- c(10)

test_mode <- TRUE
test_n <- 10000

SRS <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")

message("Project root: ", projectDir)
message("Data dir:     ", dataDir)
message("R dir:        ", rDir)
message("Results dir:  ", resultsDir)
message("AGB tiles:    ", agbTilesFolderBase)
message("Raw GFC dir:  ", gfcDir)

if (!dir.exists(agbTilesFolderBase)) {
  stop("Missing AGB tiles folder: ", agbTilesFolderBase, call. = FALSE)
}

# ----------------------------
# 1. Inputs
# ----------------------------

combined_reference_rds <- file.path(
  readyDir,
  "PVIR7_combined_reference_by_year_demo.rds"
)

if (!file.exists(combined_reference_rds)) {
  stop("Missing combined reference data: ", combined_reference_rds, call. = FALSE)
}

callaid_file <- file.path(rDir, "callAIDFlex.R")

if (!file.exists(callaid_file)) {
  stop("Missing callAIDFlex.R: ", callaid_file, call. = FALSE)
}

gfc_helper_file <- file.path(rDir, "GFC_on_the_fly.R")

if (!file.exists(gfc_helper_file)) {
  stop("Missing GFC_on_the_fly.R: ", gfc_helper_file, call. = FALSE)
}

# Source callAID first.
source(callaid_file)

# Source GFC helper after callAIDFlex so the on-the-fly functions are available globally.
source(gfc_helper_file)

reference_by_year <- readRDS(combined_reference_rds)

if (test_mode) {
  reference_by_year <- lapply(reference_by_year, function(x) {
    if (is.null(x)) return(NULL)
    x[seq_len(min(test_n, nrow(x))), , drop = FALSE]
  })
}

available_years <- names(reference_by_year)
available_years <- available_years[!is.na(as.integer(available_years))]
available_years <- sort(as.integer(available_years))

message("Available years: ", paste(available_years, collapse = ", "))

# ----------------------------
# 2. gfcanalysis download helpers
# ----------------------------

get_all_reference_points <- function(reference_by_year) {
  
  pts <- dplyr::bind_rows(
    lapply(names(reference_by_year), function(yr) {
      
      x <- reference_by_year[[yr]]
      
      if (is.null(x) || nrow(x) == 0) return(NULL)
      
      if (!all(c("POINT_X", "POINT_Y") %in% names(x))) {
        return(NULL)
      }
      
      data.frame(
        POINT_X = suppressWarnings(as.numeric(x$POINT_X)),
        POINT_Y = suppressWarnings(as.numeric(x$POINT_Y)),
        year = yr
      )
    })
  )
  
  pts <- pts %>%
    filter(
      is.finite(POINT_X),
      is.finite(POINT_Y),
      POINT_X >= -180,
      POINT_X <= 180,
      POINT_Y >= -90,
      POINT_Y <= 90
    )
  
  pts
}

make_aoi_from_points <- function(pts) {
  
  if (nrow(pts) == 0) {
    stop("No valid POINT_X/POINT_Y available to build GFC AOI.", call. = FALSE)
  }
  
  pts_sf <- sf::st_as_sf(
    pts,
    coords = c("POINT_X", "POINT_Y"),
    crs = 4326,
    remove = FALSE
  )
  
  aoi <- sf::st_as_sfc(sf::st_bbox(pts_sf))
  aoi <- sf::st_as_sf(data.frame(id = 1, geometry = aoi))
  
  aoi
}

download_gfc_tiles_for_reference <- function(reference_by_year,
                                             gfcDir,
                                             dataset = "GFC-2024-v1.12") {
  
  pts <- get_all_reference_points(reference_by_year)
  
  message("Valid points for GFC AOI: ", nrow(pts))
  
  aoi <- make_aoi_from_points(pts)
  
  message("Calculating required GFC tiles with gfcanalysis...")
  tiles <- gfcanalysis::calc_gfc_tiles(aoi)
  
  message("Required GFC tiles:")
  print(tiles)
  
  message("Downloading/checking Hansen GFC treecover2000 + lossyear...")
  
  ok <- tryCatch({
    
    gfcanalysis::download_tiles(
      tiles = tiles,
      output_folder = gfcDir,
      images = c("treecover2000", "lossyear"),
      dataset = dataset
    )
    
    TRUE
    
  }, error = function(e) {
    
    message("Download with dataset='", dataset, "' failed. Trying package default.")
    message("Original error: ", conditionMessage(e))
    
    gfcanalysis::download_tiles(
      tiles = tiles,
      output_folder = gfcDir,
      images = c("treecover2000", "lossyear")
    )
    
    TRUE
  })
  
  invisible(ok)
}

# ----------------------------
# 3. Download/check raw Hansen GFC once
# ----------------------------

download_gfc_tiles_for_reference(
  reference_by_year = reference_by_year,
  gfcDir = gfcDir,
  dataset = "GFC-2024-v1.12"
)

tc_downloaded <- list.files(
  gfcDir,
  pattern = "treecover2000.*\\.tif$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

ly_downloaded <- list.files(
  gfcDir,
  pattern = "lossyear.*\\.tif$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

message("Downloaded/available treecover2000 tiles: ", length(tc_downloaded))
message("Downloaded/available lossyear tiles:      ", length(ly_downloaded))

if (length(tc_downloaded) == 0) {
  stop("No treecover2000 files available in: ", gfcDir, call. = FALSE)
}

# ----------------------------
# 4. Filters
# ----------------------------

filter_non_lidar <- function(dat) {
  
  dat %>%
    filter(
      !CODE %in% c("LIDAR", "AFR_COF", "EMAP"),
      !grepl("_LDR$|_LIDAR$|LIDAR|COF|COFOR|EMAP", CODE, ignore.case = TRUE),
      !INVENTORY %in% c("LIDAR"),
      if ("DATA_TYPE" %in% names(.)) {
        !DATA_TYPE %in% c("LIDAR", "COF", "COFOR", "EMAP")
      } else {
        TRUE
      }
    )
}

filter_lidar <- function(dat) {
  
  dat %>%
    filter(
      CODE %in% c("LIDAR", "AFR_COF", "EMAP") |
        grepl("_LDR$|_LIDAR$|LIDAR|COF|COFOR|EMAP", CODE, ignore.case = TRUE) |
        INVENTORY %in% c("LIDAR") |
        if ("DATA_TYPE" %in% names(.)) {
          DATA_TYPE %in% c("LIDAR", "COF", "COFOR", "EMAP")
        } else {
          FALSE
        }
    )
}

# ----------------------------
# 5. Run validation
# ----------------------------

for (yr in available_years) {
  
  yr_char <- as.character(yr)
  
  if (!yr_char %in% names(reference_by_year)) {
    warning("Year not found in reference list: ", yr)
    next
  }
  
  agbTilesFolder <- file.path(agbTilesFolderBase, yr_char)
  
  if (!dir.exists(agbTilesFolder)) {
    warning("Missing AGB tile folder for year ", yr, ": ", agbTilesFolder)
    next
  }
  
  if (!is.na(sdTilesFolderBase)) {
    sdTilesFolder <- file.path(sdTilesFolderBase, yr_char)
  } else {
    sdTilesFolder <- agbTilesFolder
  }
  
  resultsFolder <- file.path(resultsDir, paste0("PVIR7_", yr_char))
  dir.create(resultsFolder, recursive = TRUE, showWarnings = FALSE)
  
  # Globals expected by AGBInvDasymetryVarMerSD.R and GFC_on_the_fly.R
  assign("pvirMapYear", as.integer(yr), envir = .GlobalEnv)
  assign("pvirProductVersion", pvirProductVersion, envir = .GlobalEnv)
  assign("pvirMapResolution", pvirMapResolution, envir = .GlobalEnv)
  assign("sdTilesFolder", sdTilesFolder, envir = .GlobalEnv)
  assign("treeCoverFolder", treeCoverFolder, envir = .GlobalEnv)
  assign("gfcDir", gfcDir, envir = .GlobalEnv)
  assign("forestTHs", forestTHs, envir = .GlobalEnv)
  assign("SRS", SRS, envir = .GlobalEnv)
  
  # Re-source helper inside loop in case AGBInv/callAID source overwrites old names.
  source(gfc_helper_file)
  
  plots_all <- reference_by_year[[yr_char]]
  
  if (is.null(plots_all) || nrow(plots_all) == 0) {
    warning("No plots for year: ", yr_char)
    next
  }
  
  plots_all$CODE <- as.character(plots_all$CODE)
  
  message("\n================================================")
  message("Running validation for map year: ", yr_char)
  message("Rows available: ", nrow(plots_all))
  message("AGB tiles:   ", agbTilesFolder)
  message("SD tiles:    ", sdTilesFolder)
  message("Raw GFC:     ", gfcDir)
  message("Results:     ", resultsFolder)
  message("================================================")
  
  # ----------------------------
  # 5.1 Plot/NFI validation
  # ----------------------------
  
  plots <- filter_non_lidar(plots_all)
  
  if (nrow(plots) > 0) {
    
    source(gfc_helper_file)
    
    callAID(
      df = plots,
      group_var = "all",
      scale = "Agg",
      minPlt = 5,
      blockRes = 0.1,
      resultsFolder = resultsFolder,
      rDir = rDir,
      agbTilesFolder = agbTilesFolder,
      treeCoverFolder = treeCoverFolder,
      forestTHs = forestTHs
    )
    
    source(gfc_helper_file)
    
    callAID(
      df = plots,
      group_var = "biome",
      scale = "Agg",
      minPlt = 5,
      blockRes = 0.1,
      resultsFolder = resultsFolder,
      rDir = rDir,
      agbTilesFolder = agbTilesFolder,
      treeCoverFolder = treeCoverFolder,
      forestTHs = forestTHs
    )
    
    mangrove_plots <- plots %>%
      filter(BIO == "Mangroves", AGB_T_HA < 100)
    
    if (nrow(mangrove_plots) > 0) {
      
      source(gfc_helper_file)
      
      callAID(
        df = mangrove_plots,
        group_var = "biome",
        scale = "Agg",
        minPlt = 2,
        blockRes = 0.1,
        resultsFolder = resultsFolder,
        rDir = rDir,
        agbTilesFolder = agbTilesFolder,
        treeCoverFolder = treeCoverFolder,
        forestTHs = forestTHs,
        result_prefix = "mangrove_lowAGB"
      )
    }
    
  } else {
    warning("No non-LIDAR plots for year: ", yr_char)
  }
  
  # ----------------------------
  # 5.2 LIDAR / special reference validation
  # ----------------------------
  
  plots_lidar <- filter_lidar(plots_all)
  
  if (nrow(plots_lidar) > 0) {
    
    plots_lidar$CODE <- "LIDAR"
    plots_lidar$varTot <- 1
    plots_lidar$TIER <- "t0"
    
    source(gfc_helper_file)
    
    callAID(
      df = plots_lidar,
      group_var = "all",
      scale = "Agg",
      minPlt = 5,
      blockRes = 0.1,
      resultsFolder = resultsFolder,
      rDir = rDir,
      agbTilesFolder = agbTilesFolder,
      treeCoverFolder = treeCoverFolder,
      forestTHs = forestTHs,
      result_prefix = "t0_lidar"
    )
    
  } else {
    message("No LIDAR/special reference plots for year: ", yr_char)
  }
  
  gc(verbose = FALSE)
}

message("\nValidation finished.")
message("Results written to: ", resultsDir)