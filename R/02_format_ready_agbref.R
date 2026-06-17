# ================================================================
# 03_make_AGBref_ready_from_valrm_with_GFC_tc.R
# Road:
#   AGBref/PVIR by-year RDS -> val.rm format -> add GFC tc -> AGBref-ready CSV
# ================================================================

suppressPackageStartupMessages({
  library(terra)
})

# ----------------------------
# 0. Settings
# ----------------------------

projectDir  <- "C:/PVIR_Reproducible"
target_year <- 2010

# Existing old output folder from previous script
oldReadyDir <- file.path(projectDir, "outputs", "pvir_ready_demo")

# Cleaner AGBref output folder
readyDir <- file.path(projectDir, "outputs", "agbref_ready_demo")
dir.create(readyDir, recursive = TRUE, showWarnings = FALSE)

gfcDir <- file.path(projectDir, "data", "GFC")

# ----------------------------
# 1. Load by-year AGBref/PVIR object
# ----------------------------

# Prefer clean AGBref name if already created; otherwise use your existing PVIR7 RDS
rds_file <- file.path(readyDir, "AGBref7_new_reference_by_year_demo.rds")

if (!file.exists(rds_file)) {
  rds_file <- file.path(oldReadyDir, "PVIR7_new_reference_by_year_demo.rds")
}

if (!file.exists(rds_file)) {
  stop("Cannot find by-year RDS file: ", rds_file, call. = FALSE)
}

agbref_by_year <- readRDS(rds_file)

# ----------------------------
# 2. Create val.rm from target MapYear
# ----------------------------

yr_chr <- as.character(target_year)

if (!yr_chr %in% names(agbref_by_year)) {
  stop(
    "MapYear ", target_year, " not found in RDS. Available years: ",
    paste(names(agbref_by_year), collapse = ", "),
    call. = FALSE
  )
}

val.rm <- agbref_by_year[[yr_chr]]
val.rm$MapYear <- target_year

# Original val.rm / AGBref columns only
# DATA_TYPE is intentionally excluded
valrm_cols <- c(
  "CODE", "AGB_T_HA", "SIZE_HA", "GEZ", "AVG_YEAR", "ZONE",
  "sdTree", "sdSE", "AGB_T_HA_ORIG", "sdGrowth", "varTot",
  "MapName", "VER", "sdMap", "BIO", "REALM", "OPEN", "INVENTORY",
  "TIER", "MapYear", "FAO.ecozone", "SD", "POINT_X", "POINT_Y"
)

# Add missing val.rm columns as NA
missing_cols <- setdiff(valrm_cols, names(val.rm))
for (nm in missing_cols) val.rm[[nm]] <- NA

# Keep only original val.rm columns and reorder
val.rm <- val.rm[, valrm_cols, drop = FALSE]

# Make sure coordinates are numeric
val.rm$POINT_X <- as.numeric(val.rm$POINT_X)
val.rm$POINT_Y <- as.numeric(val.rm$POINT_Y)

message("val.rm rows: ", nrow(val.rm))
message("val.rm columns:")
print(names(val.rm))

# ----------------------------
# 3. GFC treecover extraction
# ----------------------------

extract_gfc_tc <- function(dat,
                           gfc_dir,
                           target_year = 2010,
                           out_col = "tc") {
  
  tc_files <- list.files(
    gfc_dir,
    pattern = "treecover2000.*\\.tif$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  ly_files <- list.files(
    gfc_dir,
    pattern = "lossyear.*\\.tif$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  if (length(tc_files) == 0) {
    stop("No treecover2000 tif files found in: ", gfc_dir, call. = FALSE)
  }
  
  message("Treecover files found: ", length(tc_files))
  message("Lossyear files found:  ", length(ly_files))
  
  pts <- terra::vect(
    dat,
    geom = c("POINT_X", "POINT_Y"),
    crs = "EPSG:4326"
  )
  
  # Tree cover 2000
  tc_vrt <- file.path(gfc_dir, "treecover2000.vrt")
  terra::vrt(tc_files, filename = tc_vrt, overwrite = TRUE)
  tc_r <- terra::rast(tc_vrt)
  
  tc <- terra::extract(tc_r, pts)[, 2]
  
  # Loss-year correction up to target year
  if (length(ly_files) > 0 && target_year > 2000) {
    
    ly_vrt <- file.path(gfc_dir, "lossyear.vrt")
    terra::vrt(ly_files, filename = ly_vrt, overwrite = TRUE)
    ly_r <- terra::rast(ly_vrt)
    
    ly <- terra::extract(ly_r, pts)[, 2]
    
    # GFC lossyear: 1 = 2001, 2 = 2002, ..., 10 = 2010
    tc <- ifelse(
      !is.na(ly) & ly >= 1 & ly <= target_year - 2000,
      0,
      tc
    )
  }
  
  dat[[out_col]] <- tc
  dat
}

val.rm <- extract_gfc_tc(
  dat = val.rm,
  gfc_dir = gfcDir,
  target_year = target_year,
  out_col = "tc"
)

# Keep BIO compatible with old workflow
val.rm$BIO <- ifelse(is.na(val.rm$BIO), "NA", val.rm$BIO)

# ----------------------------
# 4. Save AGBref-ready output
# ----------------------------

out_csv <- file.path(
  readyDir,
  paste0("agbref_addv7_", target_year, "_with_tc.csv")
)

write.csv(
  val.rm,
  out_csv,
  row.names = FALSE
)

message("Saved AGBref-ready val.rm with tc:")
message(out_csv)

message("tc summary:")
print(summary(val.rm$tc))

message("Non-NA tc count: ", sum(!is.na(val.rm$tc)))

# final object in memory:
# val.rm