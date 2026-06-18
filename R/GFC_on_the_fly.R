# ================================================================
# R/GFC_on_the_fly.R
# On-the-fly Hansen GFC forest fraction for PVIR validation
#
# No writing of treecover epoch tiles.
# Uses raw gfcanalysis downloads:
#   Hansen_GFC-..._treecover2000_<tile>.tif
#   Hansen_GFC-..._lossyear_<tile>.tif
#
# Epoch logic:
#   tree cover at map year Y =
#   treecover2000 with lossyear 1:(Y - 2000) treated as non-forest
# ================================================================

suppressPackageStartupMessages({
  library(sp)
  library(raster)
})

# ----------------------------
# Utility
# ----------------------------

get_global_safe <- function(name, default = NULL) {
  if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
    get(name, envir = .GlobalEnv)
  } else {
    default
  }
}

# ----------------------------
# Old-compatible block polygon
# ----------------------------

MakeBlockPolygon <- function(x, y, size) {
  
  SRS <- get_global_safe(
    "SRS",
    sp::CRS("+proj=longlat +datum=WGS84 +no_defs")
  )
  
  xll <- size * (x %/% size)
  yll <- size * (y %/% size)
  
  pol0 <- sp::Polygon(
    cbind(
      c(xll, xll + size, xll + size, xll, xll),
      c(yll, yll, yll + size, yll + size, yll)
    )
  )
  
  pol1 <- sp::Polygons(list(pol0), "pol")
  sp::SpatialPolygons(list(pol1), proj4string = SRS)
}

# ----------------------------
# Hansen/GFC tile naming from polygon bbox
# ----------------------------

GFCtileCodes <- function(pol) {
  
  bb <- unname(sp::bbox(pol))
  crds <- expand.grid(x = bb[1, ], y = bb[2, ])
  
  codes <- character(nrow(crds))
  
  for (i in seq_len(nrow(crds))) {
    
    lon <- 10 * (crds[i, 1] %/% 10)
    lat <- 10 * (crds[i, 2] %/% 10) + 10
    
    LtX <- ifelse(lon < 0, "W", "E")
    LtY <- ifelse(lat < 0, "S", "N")
    
    WE <- paste0(sprintf("%03d", abs(lon)), LtX)
    NS <- paste0(sprintf("%02d", abs(lat)), LtY)
    
    # gfcanalysis/Hansen style, e.g. 50N_110E
    codes[i] <- paste0(NS, "_", WE)
  }
  
  unique(codes)
}

find_gfc_tile <- function(tile_code, layer = c("treecover2000", "lossyear")) {
  
  layer <- match.arg(layer)
  
  gfcDir <- get_global_safe("gfcDir", NULL)
  treeCoverFolder <- get_global_safe("treeCoverFolder", NULL)
  
  searchDir <- if (!is.null(gfcDir)) gfcDir else treeCoverFolder
  
  if (is.null(searchDir) || !dir.exists(searchDir)) {
    return(NA_character_)
  }
  
  pat <- paste0(layer, "_", tile_code, "\\.tif$")
  
  f <- list.files(
    searchDir,
    pattern = pat,
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  if (length(f) == 0) return(NA_character_)
  
  f[1]
}

TCtileNames <- function(pol) {
  
  codes <- GFCtileCodes(pol)
  
  f <- vapply(
    codes,
    function(z) find_gfc_tile(z, layer = "treecover2000"),
    character(1)
  )
  
  unique(f[!is.na(f) & file.exists(f)])
}

LYtileNames <- function(pol) {
  
  codes <- GFCtileCodes(pol)
  
  f <- vapply(
    codes,
    function(z) find_gfc_tile(z, layer = "lossyear"),
    character(1)
  )
  
  unique(f[!is.na(f) & file.exists(f)])
}

# ----------------------------
# Core on-the-fly extraction
# ----------------------------

extract_gfc_values_for_polygon <- function(pol, target_year = NULL) {
  
  if (is.null(target_year)) {
    target_year <- get_global_safe("pvirMapYear", 2010)
  }
  
  target_year <- as.integer(target_year)
  
  tc_tiles <- TCtileNames(pol)
  
  if (length(tc_tiles) == 0) {
    return(data.frame(tc = numeric(0), ly = numeric(0)))
  }
  
  all_tc <- numeric(0)
  all_ly <- numeric(0)
  
  for (tc_file in tc_tiles) {
    
    tile_code <- sub("^.*treecover2000_", "", basename(tc_file), ignore.case = TRUE)
    tile_code <- sub("\\.tif$", "", tile_code, ignore.case = TRUE)
    
    ly_file <- find_gfc_tile(tile_code, layer = "lossyear")
    
    tc_vals <- tryCatch({
      r_tc <- raster::raster(tc_file)
      r_tc_crop <- raster::crop(r_tc, raster::extent(pol))
      v <- raster::extract(r_tc_crop, pol)[[1]]
      as.numeric(v)
    }, error = function(e) {
      numeric(0)
    })
    
    if (length(tc_vals) == 0) next
    
    if (!is.na(ly_file) && file.exists(ly_file) && target_year > 2000) {
      
      ly_vals <- tryCatch({
        r_ly <- raster::raster(ly_file)
        r_ly_crop <- raster::crop(r_ly, raster::extent(pol))
        v <- raster::extract(r_ly_crop, pol)[[1]]
        as.numeric(v)
      }, error = function(e) {
        rep(NA_real_, length(tc_vals))
      })
      
      if (length(ly_vals) != length(tc_vals)) {
        ly_vals <- rep(NA_real_, length(tc_vals))
      }
      
    } else {
      ly_vals <- rep(NA_real_, length(tc_vals))
    }
    
    all_tc <- c(all_tc, tc_vals)
    all_ly <- c(all_ly, ly_vals)
  }
  
  ok <- !is.na(all_tc)
  
  data.frame(
    tc = all_tc[ok],
    ly = all_ly[ok]
  )
}

apply_epoch_loss_to_tc <- function(tc, ly, target_year) {
  
  target_year <- as.integer(target_year)
  
  if (target_year > 2000 && length(ly) == length(tc)) {
    
    loss_cutoff <- target_year - 2000
    
    tc <- ifelse(
      !is.na(ly) & ly >= 1 & ly <= loss_cutoff,
      0,
      tc
    )
  }
  
  tc
}

# ----------------------------
# Replacement for old sampleTreeCover()
# ----------------------------

sampleTreeCover <- function(pol, thresholds, wghts = FALSE) {
  
  target_year <- get_global_safe("pvirMapYear", 2010)
  
  vals <- extract_gfc_values_for_polygon(
    pol = pol,
    target_year = target_year
  )
  
  if (nrow(vals) == 0) {
    return(rep(0, length(thresholds)))
  }
  
  tc <- apply_epoch_loss_to_tc(
    tc = vals$tc,
    ly = vals$ly,
    target_year = target_year
  )
  
  tc <- tc[!is.na(tc)]
  
  if (length(tc) == 0) {
    return(rep(0, length(thresholds)))
  }
  
  TCs <- numeric(0)
  
  for (threshold in thresholds) {
    ff <- mean(ifelse(tc > threshold, 1.0, 0.0), na.rm = TRUE)
    if (!is.finite(ff)) ff <- 0
    TCs <- c(TCs, ff)
  }
  
  TCs
}

# ----------------------------
# Replacement for old sampleAGBmap()
# NOTE: despite name, old code used TC tiles here for TC_GRID stats.
# ----------------------------

sampleAGBmap <- function(pol, wghts = FALSE) {
  
  target_year <- get_global_safe("pvirMapYear", 2010)
  
  vals <- extract_gfc_values_for_polygon(
    pol = pol,
    target_year = target_year
  )
  
  if (nrow(vals) == 0) {
    return(list(NA_real_, NA_real_))
  }
  
  tc <- apply_epoch_loss_to_tc(
    tc = vals$tc,
    ly = vals$ly,
    target_year = target_year
  )
  
  tc <- tc[!is.na(tc)]
  
  if (length(tc) == 0) {
    return(list(NA_real_, NA_real_))
  }
  
  SD <- ifelse(length(tc) > 1, sd(tc, na.rm = TRUE), NA_real_)
  MEAN <- mean(tc, na.rm = TRUE)
  
  list(SD, MEAN)
}