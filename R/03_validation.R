# ================================================================
# 03_validation.R
# Run PVIR validation from combined PVIR6 + newly processed PVIR7 data
# ================================================================
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(dplyr,stringr,plotrix, terra, raster,sp,sf,plyr,foreach,doParallel,parallel,ggplot2,Metrics)



# ----------------------------
# 0. Project paths and settings
# ----------------------------

projectDir <- getwd()
dataDir    <- file.path(projectDir, "data")
rDir       <- file.path(projectDir, "R")
readyDir   <- file.path(projectDir, "outputs", "pvir_ready_demo")
resultsDir <- file.path(projectDir, "outputs", "pvir_validation")

dir.create(resultsDir, recursive = TRUE, showWarnings = FALSE)

# Large external CCI dependencies
CCImaps <- "E:/CCI_v7/AGB_GTiff_100m"

agbTilesFolderBase <- Sys.getenv(
  "PVIR_AGB_TILES_100M",
  unset = CCImaps
)

treeCoverFolder <- Sys.getenv(
  "PVIR_TREECOVER_2010",
  unset = "D:/treecover2010_v3"
)

# Use this only if SD tiles are separate from AGB tiles.
# If SD tiles are inside each yearly AGB folder, this can remain NA.
sdTilesFolderBase <- Sys.getenv(
  "PVIR_SD_TILES_100M",
  unset = NA_character_
)

pvirProductVersion <- "fv7.0"
pvirMapResolution  <- "100m"
forestTHs <- c(10)

test_mode <- T
test_n <- 10000

message("Project root: ", projectDir)
message("Data dir:     ", dataDir)
message("R dir:        ", rDir)
message("Results dir:  ", resultsDir)
message("AGB tiles:    ", agbTilesFolderBase)
message("Tree cover:   ", treeCoverFolder)

if (!dir.exists(agbTilesFolderBase)) {
  stop("Missing AGB tiles folder: ", agbTilesFolderBase, call. = FALSE)
}

if (!dir.exists(treeCoverFolder)) {
  stop("Missing tree-cover folder: ", treeCoverFolder, call. = FALSE)
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

source(file.path(rDir, "callAIDFlex.R"))

reference_by_year <- readRDS(combined_reference_rds)

if (test_mode) {
  reference_by_year <- lapply(reference_by_year, function(x) {
    if (is.null(x)) return(NULL)
    x[seq_len(min(test_n, nrow(x))), ]
  })
}

available_years <- names(reference_by_year)
available_years <- available_years[!is.na(as.integer(available_years))]
available_years <- sort(as.integer(available_years))

# ----------------------------
# 2. Filters
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
# 3. Run validation
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
  
  # Globals expected by AGBInvDasymetryVarMerSD.R
  assign("pvirMapYear", as.integer(yr), envir = .GlobalEnv)
  assign("pvirProductVersion", pvirProductVersion, envir = .GlobalEnv)
  assign("pvirMapResolution", pvirMapResolution, envir = .GlobalEnv)
  assign("sdTilesFolder", sdTilesFolder, envir = .GlobalEnv)
  
  plots_all <- reference_by_year[[yr_char]]
  plots_all$CODE <- as.character(plots_all$CODE)
  
  message("\n================================================")
  message("Running validation for map year: ", yr_char)
  message("Rows available: ", nrow(plots_all))
  message("AGB tiles: ", agbTilesFolder)
  message("SD tiles:  ", sdTilesFolder)
  message("Results:   ", resultsFolder)
  message("================================================")
  
  # ----------------------------
  # 3.1 Plot/NFI validation
  # ----------------------------
  
  plots <- filter_non_lidar(plots_all)
  
  if (nrow(plots) > 0) {
    
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
  # 3.2 LIDAR / special reference validation
  # ----------------------------
  
  plots_lidar <- filter_lidar(plots_all)
  
  if (nrow(plots_lidar) > 0) {
    
    plots_lidar$CODE <- "LIDAR"
    plots_lidar$varTot <- 1
    plots_lidar$TIER <- "t0"
    
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
}

message("\nValidation finished.")
message("Results written to: ", resultsDir)
