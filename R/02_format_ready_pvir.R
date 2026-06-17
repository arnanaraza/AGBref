# ================================================================
# 02_format_harmonize_merge_PVIR7_demo.R
# Format post-Plot2Map outputs, harmonize by map year,
# and merge newly processed PVIR7 additions with PVIR6 reference data
# ================================================================

suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(stringr)
  library(terra)
})

# ----------------------------
# 0. Project paths and settings
# ----------------------------

projectDir <- getwd()
dataDir    <- file.path(projectDir, "data")
rDir       <- file.path(projectDir, "R")
inDir      <- file.path(projectDir, "outputs", "preprocessed_demo")
readyDir   <- file.path(projectDir, "outputs", "pvir_ready_demo")

dir.create(readyDir, recursive = TRUE, showWarnings = FALSE)

pvir_version <- 7

# IMPORTANT:
# legacy_position_recycle reproduces the old mapply-style merge:
#   new list years define output names
#   old PVIR6 list is merged by list position and recycled if needed
#
# by_year is stricter and only merges when names match exactly.
merge_mode <- "legacy_position_recycle"
# merge_mode <- "by_year"

pvir6_reference_rds <- file.path(dataDir, "AGBref_PVIR6.rds")

sd_map_file <- file.path(dataDir, "ESACCI-BIOMASS-L4-AGB_SD-MERGED-1000m-fv7.0.tif")
biome_file  <- file.path(dataDir, "Ecoregions2017_biome.tif")
realm_file  <- file.path(dataDir, "Ecoregions2017_realm.tif")
bio_id_file <- file.path(dataDir, "biome_id.csv")
rea_id_file <- file.path(dataDir, "realm_id.csv")

needed_files <- c(
  pvir6_reference_rds,
  sd_map_file,
  biome_file,
  realm_file,
  bio_id_file,
  rea_id_file,
  file.path(rDir, "BiomePair.R"),
  file.path(rDir, "TempFix_NoFilter.R"),
  file.path(rDir, "TempVis.R")
)

missing_files <- needed_files[!file.exists(needed_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required file(s):\n",
    paste(missing_files, collapse = "\n"),
    call. = FALSE
  )
}

message("Project root: ", projectDir)
message("Data dir:     ", dataDir)
message("R dir:        ", rDir)
message("Input dir:    ", inDir)
message("Output dir:   ", readyDir)
message("Merge mode:   ", merge_mode)

source(file.path(rDir, "BiomePair.R"))
source(file.path(rDir, "TempFix_NoFilter.R"))
source(file.path(rDir, "TempVis.R"))

setwd(dataDir)

# ----------------------------
# 1. Read filename-based metadata
# Expected filename:
# ValidationData_DEMO_<CODE>_<DATA_TYPE>_<YEAR>.csv
# Example:
# ValidationData_DEMO_EU_NLD_LIDAR_2015.csv
# ValidationData_DEMO_ASI_MONG_NFI_2015.csv
# ValidationData_DEMO_ASI_LAO_NFI_2015.csv
# ----------------------------

derive_input_registry <- function(inDir) {
  
  files <- list.files(
    inDir,
    pattern = "^ValidationData_DEMO_[A-Z]+_[A-Z0-9]+_(LIDAR|NFI)_(19|20)\\d{2}\\.csv$",
    full.names = FALSE
  )
  
  if (length(files) == 0) {
    stop("No correctly named preprocessed files found in: ", inDir, call. = FALSE)
  }
  
  out <- lapply(files, function(f) {
    
    m <- stringr::str_match(
      f,
      "^ValidationData_DEMO_([A-Z]+_[A-Z0-9]+)_(LIDAR|NFI)_((19|20)\\d{2})\\.csv$"
    )
    
    data.frame(
      file = f,
      CODE = as.character(m[, 2]),
      DATA_TYPE = as.character(m[, 3]),
      YEAR = as.integer(m[, 4]),
      INVENTORY = ifelse(m[, 3] == "LIDAR", "local", "national"),
      stringsAsFactors = FALSE
    )
  })
  
  dplyr::bind_rows(out)
}

input_registry <- derive_input_registry(inDir)

input_files <- file.path(inDir, input_registry$file)
missing_inputs <- input_files[!file.exists(input_files)]

if (length(missing_inputs) > 0) {
  stop(
    "Missing preprocessed input file(s):\n",
    paste(missing_inputs, collapse = "\n"),
    call. = FALSE
  )
}

message("Detected input files:")
print(input_registry)

get_mapname <- function(year) {
  paste0(substr(year, 3, 4), "_CCIBiomass")
}

read_one <- function(file, CODE, DATA_TYPE, YEAR, INVENTORY) {
  
  x <- read.csv(file.path(inDir, file), stringsAsFactors = FALSE)
  names(x)[1] <- "CODE"
  x$CODE <- as.character(x$CODE)
  
  yr <- as.integer(YEAR)
  
  x %>%
    mutate(
      CODE = as.character(CODE),
      INVENTORY = as.character(INVENTORY),
      DATA_TYPE = as.character(DATA_TYPE),
      MapName = get_mapname(yr),
      VER = pvir_version,
      OPEN = 0,
      AVG_YEAR = if ("AVG_YEAR" %in% names(.)) AVG_YEAR else yr,
      AGB_T_HA_ORIG = AGB_T_HA,
      TIER = case_when(
        DATA_TYPE == "LIDAR" ~ "t0",
        SIZE_HA < 0.6 ~ "t1",
        SIZE_HA >= 0.6 & SIZE_HA < 3 ~ "t2",
        SIZE_HA >= 3 ~ "t3",
        TRUE ~ NA_character_
      )
    )
}

new_reference_raw <- dplyr::bind_rows(
  lapply(seq_len(nrow(input_registry)), function(i) {
    read_one(
      file = input_registry$file[i],
      CODE = input_registry$CODE[i],
      DATA_TYPE = input_registry$DATA_TYPE[i],
      YEAR = input_registry$YEAR[i],
      INVENTORY = input_registry$INVENTORY[i]
    )
  })
)

remove_cols <- c("FEZ", "FAO.ecozone", "RS_HA", "ratio")
new_reference_raw <- new_reference_raw[, !(names(new_reference_raw) %in% remove_cols)]

# ----------------------------
# 2. Add biome and realm
# ----------------------------

new_reference_raw <- new_reference_raw %>%
  dplyr::select(-any_of(c("BIO", "REALM", "BIO_ID", "REALM_ID")))

pts <- terra::vect(
  new_reference_raw,
  geom = c("POINT_X", "POINT_Y"),
  crs = "EPSG:4326"
)

new_reference_raw$BIO_ID <- terra::extract(terra::rast(biome_file), pts)[, 2]
new_reference_raw$REALM_ID <- terra::extract(terra::rast(realm_file), pts)[, 2]

bioID <- read.csv(bio_id_file, stringsAsFactors = FALSE)
realmID <- read.csv(rea_id_file, stringsAsFactors = FALSE)

bio_name_col <- setdiff(names(bioID), "ID")[1]
realm_name_col <- setdiff(names(realmID), "ID")[1]

names(bioID)[names(bioID) == bio_name_col] <- "BIO"
names(realmID)[names(realmID) == realm_name_col] <- "REALM"

new_reference_raw <- new_reference_raw %>%
  left_join(bioID, by = c("BIO_ID" = "ID")) %>%
  left_join(realmID, by = c("REALM_ID" = "ID")) %>%
  dplyr::select(-BIO_ID, -REALM_ID)

for (v in c("sdTree", "sdSE", "sdGrowth", "varPlot")) {
  if (!v %in% names(new_reference_raw)) new_reference_raw[[v]] <- NA_real_
}

for (v in c("ZONE", "GEZ")) {
  if (!v %in% names(new_reference_raw)) new_reference_raw[[v]] <- NA_character_
}

new_reference_pvir <- new_reference_raw %>%
  transmute(
    CODE,
    AGB_T_HA,
    SIZE_HA,
    GEZ,
    AVG_YEAR,
    ZONE,
    POINT_X,
    POINT_Y,
    sdTree,
    sdSE,
    AGB_T_HA_ORIG,
    sdGrowth,
    varTot = varPlot,
    MapName,
    VER,
    sdMap = NA_real_,
    BIO,
    REALM,
    OPEN,
    INVENTORY,
    TIER,
    DATA_TYPE
  )

write.csv(
  new_reference_pvir,
  file.path(readyDir, "PVIR7_new_reference_formatted_demo.csv"),
  row.names = FALSE
)

save(
  new_reference_pvir,
  file = file.path(readyDir, "PVIR7_new_reference_formatted_demo.RData")
)

message("New formatted reference rows: ", nrow(new_reference_pvir))

# ----------------------------
# 3. Align new data to PVIR6 structure
# Keep DATA_TYPE even if not in PVIR6 template
# ----------------------------

pvir6_reference <- readRDS(pvir6_reference_rds)
template <- pvir6_reference[[1]]

template_cols <- names(template)
extra_cols <- setdiff(names(new_reference_pvir), template_cols)
final_cols <- c(template_cols, extra_cols)

new_reference_template <- new_reference_pvir
new_reference_template$CODE <- as.character(new_reference_template$CODE)

missing_cols <- setdiff(final_cols, names(new_reference_template))
for (nm in missing_cols) new_reference_template[[nm]] <- NA

new_reference_template <- new_reference_template[, final_cols]

new_reference_template[] <- lapply(
  new_reference_template,
  function(x) if (is.list(x)) sapply(x, toString) else x
)

# ----------------------------
# 4. Load SD map stack and detect map years
# ----------------------------

sd_stack <- terra::rast(sd_map_file)

sd_years <- as.integer(stringr::str_extract(names(sd_stack), "(19|20)\\d{2}"))
sd_years <- sd_years[!is.na(sd_years)]

if (length(sd_years) == 0) {
  stop("No years detected from SD raster layer names.", call. = FALSE)
}

map_years <- sort(unique(sd_years))

message("Map years detected from SD stack:")
print(map_years)

# ----------------------------
# 5. Temporal harmonization and SD extraction for new additions
# ----------------------------

process_year <- function(dat, yr) {
  
  plt <- BiomePair(dat)
  
  gez <- sort(unique(plt$GEZ))
  gez <- gez[!is.na(gez)]
  
  if (length(gez) == 0) {
    warning("No GEZ classes found for year ", yr)
    return(plt[0, ])
  }
  
  out <- plyr::ldply(
    lapply(gez, function(z) TempApply(plt, z, yr)),
    data.frame
  )
  
  if (nrow(out) == 0) {
    warning("No rows after TempApply for year ", yr)
    return(out)
  }
  
  out <- plyr::ldply(
    lapply(gez, function(z) TempVar(out, z, yr)),
    data.frame
  )
  
  if (!"sdGrowth" %in% names(out)) out$sdGrowth <- NA_real_
  if (!"varTot" %in% names(out)) out$varTot <- out$varPlot
  
  out$sdGrowth <- ifelse(
    is.nan(out$sdGrowth),
    mean(out$sdGrowth, na.rm = TRUE),
    out$sdGrowth
  )
  
  out$varTot <- out$varTot + out$sdGrowth^2
  out$SD <- sqrt(out$varTot)
  out$MapYear <- as.integer(yr)
  
  out
}

extract_sdmap <- function(dat, yr) {
  
  if (is.null(dat) || nrow(dat) == 0) return(dat)
  
  yr <- as.integer(yr)
  
  layer_idx <- if (yr %in% sd_years) {
    which(sd_years == yr)[1]
  } else {
    which.min(abs(sd_years - yr))
  }
  
  pts <- terra::vect(
    dat,
    geom = c("POINT_X", "POINT_Y"),
    crs = "EPSG:4326"
  )
  
  dat$sdMap <- terra::extract(sd_stack[[layer_idx]], pts)[, 2]
  dat$MapYear <- yr
  
  dat
}

new_reference_by_year <- vector("list", length(map_years))
names(new_reference_by_year) <- as.character(map_years)

for (yr in map_years) {
  
  message("Processing new reference data for map year: ", yr)
  
  tmp <- process_year(new_reference_template, yr)
  tmp <- extract_sdmap(tmp, yr)
  
  new_reference_by_year[[as.character(yr)]] <- tmp
  
  message("New rows for ", yr, ": ", nrow(tmp))
}

saveRDS(
  new_reference_by_year,
  file = file.path(readyDir, "PVIR7_new_reference_by_year_demo.rds")
)

message("New rows by year:")
print(vapply(new_reference_by_year, nrow, integer(1)))

# ----------------------------
# 6. Merge PVIR6 with newly processed PVIR7 additions
# ----------------------------

harmonize_columns <- function(old_df, new_df) {
  
  if (!is.null(old_df)) old_df$CODE <- as.character(old_df$CODE)
  if (!is.null(new_df)) new_df$CODE <- as.character(new_df$CODE)
  
  if (is.null(old_df)) return(new_df)
  if (is.null(new_df)) return(old_df)
  
  missing_in_new <- setdiff(names(old_df), names(new_df))
  for (nm in missing_in_new) new_df[[nm]] <- NA
  
  missing_in_old <- setdiff(names(new_df), names(old_df))
  for (nm in missing_in_old) old_df[[nm]] <- NA
  
  new_df <- new_df[, names(old_df)]
  
  rbind(old_df, new_df)
}

if (merge_mode == "legacy_position_recycle") {
  
  output_years <- names(new_reference_by_year)
  combined_reference_by_year <- vector("list", length(output_years))
  names(combined_reference_by_year) <- output_years
  
  merge_log <- data.frame(
    output_year = character(),
    old_source_index = integer(),
    old_source_name = character(),
    old_rows = integer(),
    new_rows = integer(),
    combined_rows = integer(),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(output_years)) {
    
    yr <- output_years[i]
    
    old_i <- ((i - 1) %% length(pvir6_reference)) + 1
    
    old_df <- pvir6_reference[[old_i]]
    new_df <- new_reference_by_year[[yr]]
    
    old_source_name <- names(pvir6_reference)[old_i]
    
    if (!is.null(old_df)) {
      old_df$MapYear <- as.integer(yr)
    }
    
    if (!is.null(new_df)) {
      new_df$MapYear <- as.integer(yr)
    }
    
    out <- harmonize_columns(old_df, new_df)
    
    combined_reference_by_year[[yr]] <- out
    
    merge_log <- rbind(
      merge_log,
      data.frame(
        output_year = yr,
        old_source_index = old_i,
        old_source_name = old_source_name,
        old_rows = ifelse(is.null(old_df), 0, nrow(old_df)),
        new_rows = ifelse(is.null(new_df), 0, nrow(new_df)),
        combined_rows = ifelse(is.null(out), 0, nrow(out)),
        stringsAsFactors = FALSE
      )
    )
    
    message(
      "Merging output year ", yr,
      " | old source = ", old_source_name,
      " | old rows = ", ifelse(is.null(old_df), 0, nrow(old_df)),
      " | new rows = ", ifelse(is.null(new_df), 0, nrow(new_df)),
      " | combined = ", ifelse(is.null(out), 0, nrow(out))
    )
  }
  
} else if (merge_mode == "by_year") {
  
  all_years <- sort(unique(c(
    names(pvir6_reference),
    names(new_reference_by_year)
  )))
  
  combined_reference_by_year <- vector("list", length(all_years))
  names(combined_reference_by_year) <- all_years
  
  merge_log <- data.frame(
    output_year = character(),
    old_source_index = integer(),
    old_source_name = character(),
    old_rows = integer(),
    new_rows = integer(),
    combined_rows = integer(),
    stringsAsFactors = FALSE
  )
  
  for (yr in all_years) {
    
    old_df <- pvir6_reference[[yr]]
    new_df <- new_reference_by_year[[yr]]
    
    if (!is.null(old_df)) old_df$MapYear <- as.integer(yr)
    if (!is.null(new_df)) new_df$MapYear <- as.integer(yr)
    
    out <- harmonize_columns(old_df, new_df)
    
    combined_reference_by_year[[yr]] <- out
    
    merge_log <- rbind(
      merge_log,
      data.frame(
        output_year = yr,
        old_source_index = NA_integer_,
        old_source_name = yr,
        old_rows = ifelse(is.null(old_df), 0, nrow(old_df)),
        new_rows = ifelse(is.null(new_df), 0, nrow(new_df)),
        combined_rows = ifelse(is.null(out), 0, nrow(out)),
        stringsAsFactors = FALSE
      )
    )
    
    message(
      "Merging year ", yr,
      " | old rows = ", ifelse(is.null(old_df), 0, nrow(old_df)),
      " | new rows = ", ifelse(is.null(new_df), 0, nrow(new_df)),
      " | combined = ", ifelse(is.null(out), 0, nrow(out))
    )
  }
  
} else {
  
  stop("Unknown merge_mode: ", merge_mode, call. = FALSE)
}

# ----------------------------
# 7. Save outputs and diagnostics
# ----------------------------

saveRDS(
  combined_reference_by_year,
  file = file.path(readyDir, "PVIR7_combined_reference_by_year_demo.rds")
)

write.csv(
  merge_log,
  file.path(readyDir, "PVIR7_merge_log_demo.csv"),
  row.names = FALSE
)
setwd(projectDir)

message("Combined rows by year:")
print(vapply(combined_reference_by_year, nrow, integer(1)))

message("Merge log:")
print(merge_log)

message("Saved:")
message(file.path(readyDir, "PVIR7_new_reference_formatted_demo.csv"))
message(file.path(readyDir, "PVIR7_new_reference_by_year_demo.rds"))
message(file.path(readyDir, "PVIR7_combined_reference_by_year_demo.rds"))
message(file.path(readyDir, "PVIR7_merge_log_demo.csv"))
