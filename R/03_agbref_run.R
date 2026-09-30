# ================================================================
# 02_run_AGBref_epoch_window_temporal_adjusted.R
#
# AGBref multiresolution / multiepoch run aligned to the Sept-2026 workflow:
#   1. AVG_YEAR epoch-window selection
#   2. Target epochs 2005-2025
#   3. Temporal adjustment using BiomePair(), TempApply(), TempVar()
#   4. 500m / 1km / 10km / 25km only (no 100m)
# ================================================================

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")

pacman::p_load(
  sp,
  raster,
  plyr,
  dplyr,
  foreach,
  doParallel,
  parallel
)

# ----------------------------
# 0. Settings
# ----------------------------

projectDir <- "C:/PVIR_Reproducible"

input_csv <- "C:/PVIR_Reproducible/outputs/agbref_ready_demo/agbref_addv7_2015_with_tc.csv"

scriptsDir <- file.path(projectDir, "R")
dataDir    <- file.path(projectDir, "data")
gfcDir     <- file.path(projectDir, "data", "GFC")
outDir     <- file.path(projectDir, "outputs", "agbref_epoch_window_temporal_adjusted")

dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# Multi-epoch run used by the latest combined AGBref
target_epochs <- c(2005, 2010, 2015, 2020, 2025)

# With strict > and <, epoch_window = 11 means integer years ±10.
epoch_window <- 10

# Early-epoch exception:
# these codes are included for 1995 and 2000 even if AVG_YEAR window excludes them.
early_epoch_exception_codes <- c(
  "AUS1", "NAM2", "AFR3", "AFR11",
  "ASI_PHL", "ASI_PH",
  "CAM1", "AFR_LW"
)

# Not used in the current 2005-2025 production run.
early_exception_epochs <- integer(0)

# Production resolutions. These labels match the current combined AGBref.
# Current QC rule: n > 4 for 500m/1km/10km, n > 5 for 25km.
scales <- data.frame(
  label = c("500m", "1km", "10km", "25km"),
  aggr  = c(0.005, 0.01, 0.1, 0.25),
  minPlots = c(5, 5, 5, 6),
  stringsAsFactors = FALSE
)

forestTH <- 10

# IMPORTANT:
# The latest working AGBref validation uses the temporally harmonized,
# grid-aggregated AGB as the primary reference. Keep GFC FF as a diagnostic.
# Set TRUE only when intentionally validating wall-to-wall cell means.
apply_grid_forest_fraction <- FALSE

# If FF scaling is enabled and GFC is missing, keep the unscaled value.
fallback_to_unscaled_when_gfc_missing <- TRUE

# Latest targeted QC correction:
# AUS1, only 10km/25km, only AGB > 400 Mg/ha; use 80% of TC_GRID_MEAN as
# the extra scaling factor, without altering the stored TC_GRID_MEAN.
apply_aus1_high_agb_fix <- TRUE
aus1_tc_discount <- 0.80

# Set TRUE ONLY if input_csv is raw and these legacy unit conversions have
# not already been applied. Leave FALSE to avoid double conversion.
apply_legacy_unit_conversions <- FALSE

# Use neighboring Hansen tiles but extract only from overlapping raster extents.
use_3x3_gfc_tile_window <- TRUE

ncores <- max(1, min(4, parallel::detectCores() - 1))

SRS <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")

# ----------------------------
# 1. Source temporal adjustment scripts
# ----------------------------

source(file.path(scriptsDir, "BiomePair.R"))
source(file.path(scriptsDir, "TempFix.R"))
source(file.path(scriptsDir, "TempVis.R"))

# ----------------------------
# 2. Utilities
# ----------------------------

num_clean <- function(x) {
  if (is.numeric(x)) return(x)
  x <- as.character(x)
  x <- gsub(",", "", x)
  x <- gsub(" ", "", x)
  x <- gsub("[^0-9eE.+-]", "", x)
  suppressWarnings(as.numeric(x))
}

modalClass <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  names(sort(table(x), decreasing = TRUE))[1]
}

safe_mean <- function(x) {
  x <- num_clean(x)
  if (all(!is.finite(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

safe_sd <- function(x) {
  x <- num_clean(x)
  if (sum(is.finite(x)) <= 1) return(NA_real_)
  sd(x, na.rm = TRUE)
}

safe_wmean <- function(x, w) {
  x <- num_clean(x)
  w <- num_clean(w)
  ok <- is.finite(x) & is.finite(w) & w > 0
  
  if (sum(ok) == 0) return(mean(x, na.rm = TRUE))
  
  weighted.mean(x[ok], w[ok], na.rm = TRUE)
}

safe_inv_var <- function(x) {
  x <- num_clean(x)
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) return(NA_real_)
  1 / sum(1 / x)
}

# ----------------------------
# 3. Load input
# ----------------------------

if (exists("val.rm", envir = .GlobalEnv)) {
  rm(val.rm, envir = .GlobalEnv)
}

if (!file.exists(input_csv)) {
  stop("Input CSV does not exist: ", input_csv, call. = FALSE)
}

val.rm <- read.csv(input_csv, stringsAsFactors = FALSE)

message("Loaded input: ", input_csv)
message("Rows loaded: ", nrow(val.rm))

# Required/fallback columns
if (!"ZONE" %in% names(val.rm)) val.rm$ZONE <- "All"
if (!"BIO" %in% names(val.rm)) val.rm$BIO <- NA
if (!"GEZ" %in% names(val.rm)) val.rm$GEZ <- NA
if (!"AGB_T_HA_ORIG" %in% names(val.rm)) val.rm$AGB_T_HA_ORIG <- NA_real_
if (!"SIZE_HA" %in% names(val.rm)) val.rm$SIZE_HA <- NA_real_
if (!"varTot" %in% names(val.rm)) val.rm$varTot <- 1
if (!"varPlot" %in% names(val.rm)) val.rm$varPlot <- val.rm$varTot
if (!"OPEN" %in% names(val.rm)) val.rm$OPEN <- NA
if (!"VER" %in% names(val.rm)) val.rm$VER <- NA
if (!"INVENTORY" %in% names(val.rm)) val.rm$INVENTORY <- NA
if (!"TIER" %in% names(val.rm)) val.rm$TIER <- NA
if (!"CODE" %in% names(val.rm)) val.rm$CODE <- NA
if (!"AVG_YEAR" %in% names(val.rm)) stop("AVG_YEAR is missing.", call. = FALSE)
if (!"AGB_T_HA" %in% names(val.rm)) stop("AGB_T_HA is missing.", call. = FALSE)

val.rm$POINT_X <- num_clean(val.rm$POINT_X)
val.rm$POINT_Y <- num_clean(val.rm$POINT_Y)
val.rm$AGB_T_HA <- num_clean(val.rm$AGB_T_HA)
val.rm$AGB_T_HA_ORIG <- num_clean(val.rm$AGB_T_HA_ORIG)
val.rm$AVG_YEAR <- num_clean(val.rm$AVG_YEAR)
val.rm$SIZE_HA <- num_clean(val.rm$SIZE_HA)
val.rm$varTot <- num_clean(val.rm$varTot)
val.rm$varPlot <- num_clean(val.rm$varPlot)

# Row-wise fallback
bad_orig <- !is.finite(val.rm$AGB_T_HA_ORIG) & is.finite(val.rm$AGB_T_HA)
val.rm$AGB_T_HA_ORIG[bad_orig] <- val.rm$AGB_T_HA[bad_orig]

# If varPlot missing, use varTot
bad_varplot <- !is.finite(val.rm$varPlot) & is.finite(val.rm$varTot)
val.rm$varPlot[bad_varplot] <- val.rm$varTot[bad_varplot]

val.rm$BIO <- ifelse(is.na(val.rm$BIO), "NA", val.rm$BIO)

# Optional legacy unit conversions. Run exactly once on raw inputs.
if (isTRUE(apply_legacy_unit_conversions)) {
  idx049 <- val.rm$CODE %in% c("ASI_PH", "SAM_guy", "SAM_ECU")
  val.rm$AGB_T_HA[idx049] <- val.rm$AGB_T_HA[idx049] / 0.49
  val.rm$AGB_T_HA_ORIG[idx049] <- val.rm$AGB_T_HA_ORIG[idx049] / 0.49

  idxjap <- val.rm$CODE == "ASI_JAP"
  val.rm$AGB_T_HA[idxjap] <- val.rm$AGB_T_HA[idxjap] / 0.2
  val.rm$AGB_T_HA_ORIG[idxjap] <- val.rm$AGB_T_HA_ORIG[idxjap] / 0.2
}

# No legacy filtering here.
# No ASI_JAP removal.
# No PH/Guyana/Ecuador correction here unless you intentionally add it.

# ----------------------------
# 4. GFC tile helpers
# ----------------------------

MakeBlockPolygon <- function(x, y, size) {
  xll <- size * (x %/% size)
  yll <- size * (y %/% size)
  
  pol0 <- sp::Polygon(
    cbind(
      c(xll, xll + size, xll + size, xll, xll),
      c(yll, yll, yll + size, yll + size, yll)
    )
  )
  
  sp::SpatialPolygons(
    list(sp::Polygons(list(pol0), "pol")),
    proj4string = SRS
  )
}
gfc_code_from_xy <- function(x, y) {
  lon <- 10 * (x %/% 10)
  lat <- 10 * (y %/% 10) + 10
  
  LtX <- ifelse(lon < 0, "W", "E")
  LtY <- ifelse(lat < 0, "S", "N")
  
  paste0(
    sprintf("%02d", abs(lat)), LtY,
    "_",
    sprintf("%03d", abs(lon)), LtX
  )
}
parse_gfc_code <- function(code) {
  lat_part <- sub("_.*$", "", code)
  lon_part <- sub("^.*_", "", code)
  
  lat_val <- as.numeric(substr(lat_part, 1, nchar(lat_part) - 1))
  lat_hemi <- substr(lat_part, nchar(lat_part), nchar(lat_part))
  
  lon_val <- as.numeric(substr(lon_part, 1, nchar(lon_part) - 1))
  lon_hemi <- substr(lon_part, nchar(lon_part), nchar(lon_part))
  
  data.frame(
    lat_edge = ifelse(lat_hemi == "S", -lat_val, lat_val),
    lon_base = ifelse(lon_hemi == "W", -lon_val, lon_val)
  )
}
format_gfc_code <- function(lat_edge, lon_base) {
  LtX <- ifelse(lon_base < 0, "W", "E")
  LtY <- ifelse(lat_edge < 0, "S", "N")
  
  paste0(
    sprintf("%02d", abs(lat_edge)), LtY,
    "_",
    sprintf("%03d", abs(lon_base)), LtX
  )
}
expand_gfc_codes_3x3 <- function(codes) {
  out <- character(0)
  
  for (cd in unique(codes)) {
    p <- parse_gfc_code(cd)
    
    grid <- expand.grid(
      lat_edge = p$lat_edge + c(-10, 0, 10),
      lon_base = p$lon_base + c(-10, 0, 10)
    )
    
    out <- c(out, mapply(format_gfc_code, grid$lat_edge, grid$lon_base))
  }
  
  unique(out)
}
GFCtileCodes_core <- function(pol) {
  bb <- unname(sp::bbox(pol))
  
  xs <- unique(c(bb[1, 1], bb[1, 2], mean(bb[1, ])))
  ys <- unique(c(bb[2, 1], bb[2, 2], mean(bb[2, ])))
  
  xy <- expand.grid(x = xs, y = ys)
  
  unique(mapply(gfc_code_from_xy, xy$x, xy$y))
}
GFCtileCodes <- function(pol) {
  codes <- GFCtileCodes_core(pol)
  
  if (isTRUE(use_3x3_gfc_tile_window)) {
    codes <- expand_gfc_codes_3x3(codes)
  }
  
  unique(codes)
}
find_gfc_tile <- function(tile_code, layer = c("treecover2000", "lossyear")) {
  layer <- match.arg(layer)
  
  f <- list.files(
    gfcDir,
    pattern = paste0(layer, "_", tile_code, "\\.tif$"),
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  if (length(f) == 0) return(NA_character_)
  f[1]
}
extent_overlaps <- function(r_ext, p_ext) {
  !(r_ext@xmax <= p_ext@xmin ||
      r_ext@xmin >= p_ext@xmax ||
      r_ext@ymax <= p_ext@ymin ||
      r_ext@ymin >= p_ext@ymax)
}
get_tc_tiles <- function(pol) {
  codes <- GFCtileCodes(pol)
  
  f <- vapply(codes, function(z) find_gfc_tile(z, "treecover2000"), character(1))
  f <- unique(f[!is.na(f) & file.exists(f)])
  
  if (length(f) == 0) return(character(0))
  
  p_ext <- raster::extent(pol)
  
  keep <- vapply(
    f,
    function(ff) {
      tryCatch({
        extent_overlaps(raster::extent(raster::raster(ff)), p_ext)
      }, error = function(e) FALSE)
    },
    logical(1)
  )
  
  f[keep]
}
# ----------------------------
# 5. GFC extraction
# ----------------------------

apply_epoch_loss_to_tc <- function(tc, ly, target_year) {
  target_year <- as.integer(target_year)
  
  if (target_year > 2000 && length(tc) == length(ly)) {
    loss_cutoff <- target_year - 2000
    
    tc <- ifelse(
      !is.na(ly) & ly >= 1 & ly <= loss_cutoff,
      0,
      tc
    )
  }
  
  tc
}

extract_gfc_polygon <- function(pol, target_year) {
  tc_tiles <- get_tc_tiles(pol)
  
  if (length(tc_tiles) == 0) {
    return(list(
      ff = NA_real_,
      tc_mean = NA_real_,
      tc_sd = NA_real_,
      pixel_n = 0,
      status = "NO_TC_TILE"
    ))
  }
  
  all_tc <- numeric(0)
  all_ly <- numeric(0)
  status <- "OK"
  p_ext <- raster::extent(pol)
  
  for (tc_file in tc_tiles) {
    tile_code <- sub("^.*treecover2000_", "", basename(tc_file), ignore.case = TRUE)
    tile_code <- sub("\\.tif$", "", tile_code, ignore.case = TRUE)
    
    ly_file <- find_gfc_tile(tile_code, "lossyear")
    
    tc_vals <- tryCatch({
      r_tc <- raster::raster(tc_file)
      if (!extent_overlaps(raster::extent(r_tc), p_ext)) return(numeric(0))
      r_tc_crop <- raster::crop(r_tc, p_ext)
      as.numeric(raster::extract(r_tc_crop, pol)[[1]])
    }, error = function(e) {
      status <<- "TC_ERROR"
      numeric(0)
    })
    
    if (length(tc_vals) == 0) next
    
    ly_vals <- rep(NA_real_, length(tc_vals))
    
    if (!is.na(ly_file) && file.exists(ly_file) && target_year > 2000) {
      ly_vals <- tryCatch({
        r_ly <- raster::raster(ly_file)
        if (!extent_overlaps(raster::extent(r_ly), p_ext)) return(rep(NA_real_, length(tc_vals)))
        r_ly_crop <- raster::crop(r_ly, p_ext)
        v <- as.numeric(raster::extract(r_ly_crop, pol)[[1]])
        if (length(v) != length(tc_vals)) rep(NA_real_, length(tc_vals)) else v
      }, error = function(e) {
        status <<- "LY_ERROR"
        rep(NA_real_, length(tc_vals))
      })
    }
    
    all_tc <- c(all_tc, tc_vals)
    all_ly <- c(all_ly, ly_vals)
  }
  
  ok <- !is.na(all_tc)
  all_tc <- all_tc[ok]
  all_ly <- all_ly[ok]
  
  if (length(all_tc) == 0) {
    return(list(
      ff = NA_real_,
      tc_mean = NA_real_,
      tc_sd = NA_real_,
      pixel_n = 0,
      status = "NO_TC_PIXELS"
    ))
  }
  
  tc_epoch <- apply_epoch_loss_to_tc(all_tc, all_ly, target_year)
  tc_epoch <- tc_epoch[!is.na(tc_epoch)]
  
  if (length(tc_epoch) == 0) {
    return(list(
      ff = NA_real_,
      tc_mean = NA_real_,
      tc_sd = NA_real_,
      pixel_n = 0,
      status = "NO_TC_AFTER_LOSS"
    ))
  }
  
  ff <- mean(tc_epoch > forestTH, na.rm = TRUE)
  ff <- pmin(pmax(ff, 0), 1)
  
  list(
    ff = ff,
    tc_mean = mean(tc_epoch, na.rm = TRUE),
    tc_sd = ifelse(length(tc_epoch) > 1, sd(tc_epoch, na.rm = TRUE), NA_real_),
    pixel_n = length(tc_epoch),
    status = status
  )
}

extract_tc_points_epoch <- function(dat, target_year) {
  out <- rep(NA_real_, nrow(dat))
  dat$.row_id <- seq_len(nrow(dat))
  dat$.tile_code <- gfc_code_from_xy(dat$POINT_X, dat$POINT_Y)
  
  for (tile_code in unique(dat$.tile_code)) {
    ids <- which(dat$.tile_code == tile_code)
    
    tc_file <- find_gfc_tile(tile_code, "treecover2000")
    ly_file <- find_gfc_tile(tile_code, "lossyear")
    
    if (is.na(tc_file) || !file.exists(tc_file)) next
    
    pts <- sp::SpatialPoints(dat[ids, c("POINT_X", "POINT_Y")], proj4string = SRS)
    
    tc_vals <- tryCatch(
      raster::extract(raster::raster(tc_file), pts),
      error = function(e) rep(NA_real_, length(ids))
    )
    
    if (!is.na(ly_file) && file.exists(ly_file) && target_year > 2000) {
      ly_vals <- tryCatch(
        raster::extract(raster::raster(ly_file), pts),
        error = function(e) rep(NA_real_, length(ids))
      )
      
      if (length(ly_vals) == length(tc_vals)) {
        tc_vals <- ifelse(
          !is.na(ly_vals) & ly_vals >= 1 & ly_vals <= target_year - 2000,
          0,
          tc_vals
        )
      }
    }
    
    out[dat$.row_id[ids]] <- tc_vals
  }
  
  out
}

# ----------------------------
# 6. Temporal adjustment
# ----------------------------

apply_temporal_adjustment <- function(dat, epoch) {

  if (nrow(dat) == 0) return(dat)

  # Stable row key so TempApply/TempVar output can be checked/restored.
  dat$.ROW_ID <- seq_len(nrow(dat))

  dat <- BiomePair(dat)

  if (!".ROW_ID" %in% names(dat)) {
    stop("BiomePair() dropped .ROW_ID; temporal adjustment cannot be safely aligned.",
         call. = FALSE)
  }

  gez <- sort(unique(dat$GEZ))
  gez <- gez[!is.na(gez)]

  if (length(gez) == 0) {
    warning("No GEZ found for epoch ", epoch, ". Returning unadjusted filtered data.")
    dat$MapYear <- epoch
    return(dat)
  }

  dat2 <- plyr::ldply(
    lapply(gez, function(z) TempApply(dat, z, epoch)),
    data.frame
  )

  if (is.null(dat2) || nrow(dat2) == 0) {
    warning("TempApply returned zero rows for epoch ", epoch)
    return(dat[0, ])
  }

  if (!".ROW_ID" %in% names(dat2)) {
    stop("TempApply() dropped .ROW_ID; fix TempFix.R before continuing.",
         call. = FALSE)
  }

  # Preserve TempApply-adjusted AGB by row ID before TempVar.
  agb_after_tempapply <- setNames(
    num_clean(dat2$AGB_T_HA),
    as.character(dat2$.ROW_ID)
  )

  dat3 <- plyr::ldply(
    lapply(gez, function(z) TempVar(dat2, z, epoch)),
    data.frame
  )

  if (is.null(dat3) || nrow(dat3) == 0) {
    warning("TempVar returned zero rows for epoch ", epoch)
    dat3 <- dat2
  }

  if (!".ROW_ID" %in% names(dat3)) {
    stop("TempVar() dropped .ROW_ID; fix TempFix.R before continuing.",
         call. = FALSE)
  }

  # TempVar should alter uncertainty, not the temporally adjusted biomass.
  # Restore AGB_T_HA explicitly by stable row ID. This avoids the known
  # row-order reassignment problem in older TempVar implementations.
  m <- match(as.character(dat3$.ROW_ID), names(agb_after_tempapply))
  okm <- !is.na(m)
  dat3$AGB_T_HA[okm] <- unname(agb_after_tempapply[m[okm]])

  if (!"sdGrowth" %in% names(dat3)) dat3$sdGrowth <- NA_real_
  if (!"varTot" %in% names(dat3)) dat3$varTot <- dat3$varPlot
  if (!"varPlot" %in% names(dat3)) dat3$varPlot <- dat3$varTot

  dat3$sdGrowth <- num_clean(dat3$sdGrowth)
  fill_growth <- mean(dat3$sdGrowth[is.finite(dat3$sdGrowth)], na.rm = TRUE)
  if (!is.finite(fill_growth)) fill_growth <- 0
  dat3$sdGrowth[!is.finite(dat3$sdGrowth)] <- fill_growth

  dat3$varTot <- num_clean(dat3$varTot)
  dat3$varPlot <- num_clean(dat3$varPlot)
  bad_v <- !is.finite(dat3$varTot)
  dat3$varTot[bad_v] <- dat3$varPlot[bad_v]
  dat3$varTot[!is.finite(dat3$varTot)] <- 1

  dat3$varTot <- dat3$varTot + dat3$sdGrowth^2
  dat3$SD <- sqrt(dat3$varTot)
  dat3$MapYear <- epoch

  dat3
}

# ----------------------------
# 7. Main aggregation
# ----------------------------

make_agbref_epoch_resolution <- function(dat, epoch, res_label, aggr, minPlots = 1) {

  dat <- dat %>%
    dplyr::filter(
      is.finite(POINT_X),
      is.finite(POINT_Y),
      is.finite(AGB_T_HA),
      is.finite(AVG_YEAR)
    ) %>%
    dplyr::filter(
      (AVG_YEAR > epoch - epoch_window &
         AVG_YEAR < epoch + epoch_window) |
        (epoch %in% early_exception_epochs &
           CODE %in% early_epoch_exception_codes)
    )

  if (nrow(dat) == 0) return(data.frame())

  dat$AGB_T_HA_ORIG[!is.finite(dat$AGB_T_HA_ORIG)] <-
    dat$AGB_T_HA[!is.finite(dat$AGB_T_HA_ORIG)]

  dat$varTot[!is.finite(dat$varTot)] <- dat$varPlot[!is.finite(dat$varTot)]
  dat$varTot[!is.finite(dat$varTot)] <- 1

  # Full temporal harmonization to target epoch.
  dat <- apply_temporal_adjustment(dat, epoch)

  if (nrow(dat) == 0) return(data.frame())

  # Clean again after temporal scripts.
  dat$POINT_X <- num_clean(dat$POINT_X)
  dat$POINT_Y <- num_clean(dat$POINT_Y)
  dat$AGB_T_HA <- num_clean(dat$AGB_T_HA)
  dat$AGB_T_HA_ORIG <- num_clean(dat$AGB_T_HA_ORIG)
  dat$AVG_YEAR <- num_clean(dat$AVG_YEAR)
  dat$varTot <- num_clean(dat$varTot)

  if (!"FAO.ecozone" %in% names(dat)) dat$FAO.ecozone <- NA_character_
  if (!"ZONE" %in% names(dat)) dat$ZONE <- NA_character_
  if (!"GEZ" %in% names(dat)) dat$GEZ <- NA_character_
  if (!"sdGrowth" %in% names(dat)) dat$sdGrowth <- 0
  if (!"SD" %in% names(dat)) dat$SD <- sqrt(dat$varTot)

  dat <- dat %>%
    dplyr::filter(
      is.finite(POINT_X),
      is.finite(POINT_Y),
      is.finite(AGB_T_HA)
    )

  if (nrow(dat) == 0) return(data.frame())

  dat$AGB_T_HA_ORIG[!is.finite(dat$AGB_T_HA_ORIG)] <-
    dat$AGB_T_HA[!is.finite(dat$AGB_T_HA_ORIG)]

  # Point-level epoch tree cover is retained for TC_PLT summaries.
  dat$tc <- extract_tc_points_epoch(dat, epoch)

  dat$Xnew <- aggr * (0.5 + dat$POINT_X %/% aggr)
  dat$Ynew <- aggr * (0.5 + dat$POINT_Y %/% aggr)
  dat$inv <- ifelse(is.finite(dat$varTot) & dat$varTot > 0,
                    1 / dat$varTot, NA_real_)

  cells <- dat %>%
    dplyr::group_by(Xnew, Ynew) %>%
    dplyr::summarise(
      POINT_X = dplyr::first(Xnew),
      POINT_Y = dplyr::first(Ynew),
      n = sum(is.finite(AGB_T_HA)),

      # Keep original and temporally harmonized quantities separate.
      AGB_T_HA_ORIG = safe_mean(AGB_T_HA_ORIG),
      AGB_T_HA_PRE_FF = safe_wmean(AGB_T_HA, inv),

      SIZE_HA = safe_mean(SIZE_HA),

      TC_PLT_MEAN = safe_mean(tc),
      TC_PLT_SD = safe_sd(tc),

      AVG_YEAR = round(safe_mean(AVG_YEAR)),
      MapYear = epoch,

      BIO = modalClass(BIO),
      CODE = modalClass(CODE),
      INVENTORY = modalClass(INVENTORY),
      TIER = modalClass(TIER),
      OPEN = modalClass(OPEN),
      VER = modalClass(VER),

      ZONE = modalClass(ZONE),
      FAO.ecozone = modalClass(FAO.ecozone),
      GEZ = modalClass(GEZ),

      sdGrowth = safe_mean(sdGrowth),
      varTot = safe_inv_var(varTot),
      .groups = "drop"
    ) %>%
    dplyr::filter(n >= minPlots)

  if (nrow(cells) == 0) return(data.frame())

  cells$SD <- sqrt(cells$varTot)

  cl <- parallel::makeCluster(max(1, min(ncores, nrow(cells))))
  doParallel::registerDoParallel(cl)

  gfc_df <- foreach(
    i = seq_len(nrow(cells)),
    .combine = "rbind",
    .packages = c("sp", "raster"),
    .export = c(
      "SRS", "gfcDir", "forestTH", "use_3x3_gfc_tile_window",
      "MakeBlockPolygon", "gfc_code_from_xy", "parse_gfc_code",
      "format_gfc_code", "expand_gfc_codes_3x3", "GFCtileCodes_core",
      "GFCtileCodes", "find_gfc_tile", "extent_overlaps", "get_tc_tiles",
      "apply_epoch_loss_to_tc", "extract_gfc_polygon"
    )
  ) %dopar% {
    pol <- MakeBlockPolygon(cells$POINT_X[i], cells$POINT_Y[i], aggr)
    z <- extract_gfc_polygon(pol, epoch)

    data.frame(
      FF_USED = z$ff,
      TC_GRID_MEAN = z$tc_mean,
      TC_GRID_SD = z$tc_sd,
      GFC_PIXEL_N = z$pixel_n,
      GFC_STATUS = z$status,
      stringsAsFactors = FALSE
    )
  }

  parallel::stopCluster(cl)
  foreach::registerDoSEQ()

  out <- dplyr::bind_cols(cells, gfc_df)

  # Primary AGB:
  # latest working AGBref keeps the temporally harmonized aggregated value
  # as the reference and retains FF as a diagnostic.
  if (isTRUE(apply_grid_forest_fraction)) {
    out$AGB_T_HA <- dplyr::case_when(
      is.finite(out$FF_USED) ~ out$AGB_T_HA_PRE_FF * out$FF_USED,
      !is.finite(out$FF_USED) &
        isTRUE(fallback_to_unscaled_when_gfc_missing) ~ out$AGB_T_HA_PRE_FF,
      TRUE ~ NA_real_
    )

    out$CELL_STATUS <- dplyr::case_when(
      is.finite(out$FF_USED) ~ "FOREST_SCALED",
      !is.finite(out$FF_USED) ~ "UNSCALED_GFC_MISSING",
      TRUE ~ "OTHER"
    )
  } else {
    out$AGB_T_HA <- out$AGB_T_HA_PRE_FF
    out$CELL_STATUS <- "UNSCALED_PRIMARY"
  }

  # Targeted current-QC correction. Do not overwrite TC_GRID_MEAN.
  out$AUS1_SPECIAL_CORR <- FALSE
  if (isTRUE(apply_aus1_high_agb_fix) &&
      res_label %in% c("10km", "25km")) {

    idx_aus <- out$CODE == "AUS1" &
      out$AGB_T_HA > 400 &
      is.finite(out$TC_GRID_MEAN)

    out$AGB_T_HA[idx_aus] <-
      out$AGB_T_HA[idx_aus] *
      ((out$TC_GRID_MEAN[idx_aus] * aus1_tc_discount) / 100)

    out$AUS1_SPECIAL_CORR[idx_aus] <- TRUE
  }

  # Output schema aligned to the latest combined AGBref, with extra QC fields.
  out %>%
    dplyr::transmute(
      TC_PLT_SD,
      TC_PLT_MEAN,
      TC_GRID_SD,
      TC_GRID_MEAN,
      n,
      AGB_T_HA_ORIG,
      SIZE_HA,
      OPEN,
      VER,
      varTot,
      AVG_YEAR,
      BIO,
      CODE,
      INVENTORY,
      TIER,
      POINT_X,
      POINT_Y,
      ZONE,
      FAO.ecozone,
      GEZ,
      sdGrowth,
      SD,
      Resolution = res_label,
      Year = epoch,
  #    AGB_T_HA_PRE_FF,
   #   FF_USED,
    #  GFC_PIXEL_N,
     # GFC_STATUS,
      #CELL_STATUS,
      )
}

# ----------------------------
# 8. Run
# ----------------------------
setwd(dataDir)
data_frames <- list()
agbref <- data.frame()

for (ep in target_epochs) {
  for (j in seq_len(nrow(scales))) {
    
    res_label <- scales$label[j]
    aggr <- scales$aggr[j]
    minPlots <- scales$minPlots[j]
    
    message("Running: epoch=", ep, " | resolution=", res_label)
    
    x <- make_agbref_epoch_resolution(
      dat = val.rm,
      epoch = ep,
      res_label = res_label,
      aggr = aggr,
      minPlots = minPlots
    )
    
    nm <- paste(res_label, ep, sep = "_")
    data_frames[[nm]] <- x
    
    if (nrow(x) > 0) {
      agbref <- dplyr::bind_rows(agbref, x)
      names(agbref)[names(agbref) == "AGB_T_HA_ORIG"] <- "AGB_T_HA"
      write.csv(
        x,
        file.path(outDir, paste0("AGBref_", res_label, "_", ep, ".csv")),
        row.names = FALSE
      )
    }
  }
}


