# ================================================================
# 02_agbref_ready_multiepoch.R
# Multi-resolution x multi-epoch AGBref-ready aggregation
#
# Uses old-style exact tile sampling:
#   - MakeBlockPolygon()
#   - tile lookup from polygon bbox
#   - raster::crop(raster(tile), extent(pol))
#   - raster::extract(..., pol)
#   - foreach parallel over aggregated cells
#
# No temporal biomass adjustment here.
# Epoch is used only for GFC lossyear correction.
# ================================================================

suppressPackageStartupMessages({
  library(sp)
  library(raster)
  library(dplyr)
  library(plyr)
  library(foreach)
  library(doParallel)
  library(parallel)
})

# ----------------------------
# 0. Paths and settings
# ----------------------------

projectDir <- getwd()
dataDir    <- file.path(projectDir, "data")
readyDir   <- file.path(projectDir, "outputs", "pvir_ready_demo")
outDir     <- file.path(projectDir, "outputs", "agbref_ready_multiepoch")
gfcDir     <- file.path(projectDir, "data", "GFC")

dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

treeCoverFolder <- gfcDir
lossYearFolder  <- gfcDir

forestTHs <- c(10)

target_epochs <- c(2005, 2010, 2015, 2020)

scales <- data.frame(
  label = c("100m", "500m", "1km", "10km", "25km"),
  aggr  = c(0.001, 0.005, 0.01, 0.1, 0.25),
  minPlots = c(1, 1, 1, 1, 1),
  stringsAsFactors = FALSE
)

# Use controlled parallelism; too many workers can duplicate raster reads.
ncores <- max(1, min(4, parallel::detectCores() - 1))

SRS <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")

message("Project root:     ", projectDir)
message("Output dir:       ", outDir)
message("GFC dir:          ", gfcDir)
message("Parallel workers: ", ncores)

# ----------------------------
# 1. Utility functions
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

safe_inv_var <- function(x) {
  x <- num_clean(x)
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) return(NA_real_)
  1 / sum(1 / x)
}

safe_wmean <- function(x, w) {
  x <- num_clean(x)
  w <- num_clean(w)
  ok <- is.finite(x) & is.finite(w) & w > 0
  
  if (sum(ok) == 0) return(mean(x, na.rm = TRUE))
  
  weighted.mean(x[ok], w[ok], na.rm = TRUE)
}

safe_sd <- function(x) {
  x <- num_clean(x)
  if (sum(is.finite(x)) <= 1) return(NA_real_)
  sd(x, na.rm = TRUE)
}

safe_mean <- function(x) {
  x <- num_clean(x)
  if (all(!is.finite(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

# ----------------------------
# 2. Old-style polygon and tile lookup
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
  
  pol1 <- sp::Polygons(list(pol0), "pol")
  sp::SpatialPolygons(list(pol1), proj4string = SRS)
}

gfc_tile_codes_from_pol <- function(pol) {
  
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
    
    # Hansen/gfcanalysis format: 50N_110E
    codes[i] <- paste0(NS, "_", WE)
  }
  
  unique(codes)
}

find_gfc_tile <- function(tile_code, layer = c("treecover2000", "lossyear")) {
  
  layer <- match.arg(layer)
  
  pat <- paste0(layer, "_", tile_code, "\\.tif$")
  
  f <- list.files(
    gfcDir,
    pattern = pat,
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  
  if (length(f) == 0) return(NA_character_)
  
  f[1]
}

TCtileNames_gfc <- function(pol) {
  codes <- gfc_tile_codes_from_pol(pol)
  f <- vapply(codes, find_gfc_tile, character(1), layer = "treecover2000")
  unique(f[!is.na(f) & file.exists(f)])
}

LYtileNames_gfc <- function(pol) {
  codes <- gfc_tile_codes_from_pol(pol)
  f <- vapply(codes, find_gfc_tile, character(1), layer = "lossyear")
  unique(f[!is.na(f) & file.exists(f)])
}

# ----------------------------
# 3. GFC sampling old-style, per exact tile
# ----------------------------

sampleTreeCover_gfc <- function(pol,
                                thresholds,
                                target_year = 2010,
                                wghts = FALSE) {
  
  tc_tiles <- TCtileNames_gfc(pol)
  
  if (length(tc_tiles) == 0) {
    return(rep(0, length(thresholds)))
  }
  
  all_tc <- numeric(0)
  all_ly <- numeric(0)
  
  for (tc_file in tc_tiles) {
    
    tile_code <- sub("^.*treecover2000_", "", basename(tc_file))
    tile_code <- sub("\\.tif$", "", tile_code, ignore.case = TRUE)
    
    ly_file <- find_gfc_tile(tile_code, layer = "lossyear")
    
    tc_vals <- NULL
    ly_vals <- NULL
    
    # Extract treecover
    tc_vals <- tryCatch({
      r_tc <- raster::raster(tc_file)
      r_tc_crop <- raster::crop(r_tc, raster::extent(pol))
      v <- raster::extract(r_tc_crop, pol)[[1]]
      as.numeric(v)
    }, error = function(e) {
      numeric(0)
    })
    
    if (length(tc_vals) == 0) next
    
    # Extract lossyear for same polygon/tile
    if (!is.na(ly_file) && file.exists(ly_file) && target_year > 2000) {
      ly_vals <- tryCatch({
        r_ly <- raster::raster(ly_file)
        r_ly_crop <- raster::crop(r_ly, raster::extent(pol))
        v <- raster::extract(r_ly_crop, pol)[[1]]
        as.numeric(v)
      }, error = function(e) {
        rep(NA_real_, length(tc_vals))
      })
      
      # If mismatch happens, ignore loss adjustment for this tile rather than crash.
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
  all_tc <- all_tc[ok]
  all_ly <- all_ly[ok]
  
  if (length(all_tc) == 0) {
    return(rep(0, length(thresholds)))
  }
  
  if (target_year > 2000 && length(all_ly) == length(all_tc)) {
    all_tc <- ifelse(
      !is.na(all_ly) & all_ly >= 1 & all_ly <= target_year - 2000,
      0,
      all_tc
    )
  }
  
  out <- numeric(0)
  
  for (threshold in thresholds) {
    ff <- mean(ifelse(all_tc > threshold, 1.0, 0.0), na.rm = TRUE)
    if (!is.finite(ff)) ff <- 0
    out <- c(out, ff)
  }
  
  out
}

sampleTCStats_gfc <- function(pol,
                              target_year = 2010) {
  
  tc_tiles <- TCtileNames_gfc(pol)
  
  if (length(tc_tiles) == 0) {
    return(list(SD = NA_real_, MEAN = NA_real_))
  }
  
  all_tc <- numeric(0)
  all_ly <- numeric(0)
  
  for (tc_file in tc_tiles) {
    
    tile_code <- sub("^.*treecover2000_", "", basename(tc_file))
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
  all_tc <- all_tc[ok]
  all_ly <- all_ly[ok]
  
  if (length(all_tc) == 0) {
    return(list(SD = NA_real_, MEAN = NA_real_))
  }
  
  if (target_year > 2000 && length(all_ly) == length(all_tc)) {
    all_tc <- ifelse(
      !is.na(all_ly) & all_ly >= 1 & all_ly <= target_year - 2000,
      0,
      all_tc
    )
  }
  
  list(
    SD = ifelse(length(all_tc) > 1, sd(all_tc, na.rm = TRUE), NA_real_),
    MEAN = mean(all_tc, na.rm = TRUE)
  )
}

extract_tc_points_gfc <- function(dat,
                                  target_year = 2010) {
  
  tc_out <- rep(NA_real_, nrow(dat))
  
  # group points by GFC tile code
  tmp <- dat
  tmp$.row_id <- seq_len(nrow(tmp))
  
  tmp$.tile_lon <- 10 * (tmp$POINT_X %/% 10)
  tmp$.tile_lat <- 10 * (tmp$POINT_Y %/% 10) + 10
  
  tmp$.LtX <- ifelse(tmp$.tile_lon < 0, "W", "E")
  tmp$.LtY <- ifelse(tmp$.tile_lat < 0, "S", "N")
  
  tmp$.tile_code <- paste0(
    sprintf("%02d", abs(tmp$.tile_lat)), tmp$.LtY,
    "_",
    sprintf("%03d", abs(tmp$.tile_lon)), tmp$.LtX
  )
  
  for (tile_code in unique(tmp$.tile_code)) {
    
    ids <- which(tmp$.tile_code == tile_code)
    
    tc_file <- find_gfc_tile(tile_code, layer = "treecover2000")
    ly_file <- find_gfc_tile(tile_code, layer = "lossyear")
    
    if (is.na(tc_file) || !file.exists(tc_file)) next
    
    pts <- sp::SpatialPoints(
      tmp[ids, c("POINT_X", "POINT_Y")],
      proj4string = SRS
    )
    
    tc_vals <- tryCatch({
      raster::extract(raster::raster(tc_file), pts)
    }, error = function(e) {
      rep(NA_real_, length(ids))
    })
    
    if (!is.na(ly_file) && file.exists(ly_file) && target_year > 2000) {
      ly_vals <- tryCatch({
        raster::extract(raster::raster(ly_file), pts)
      }, error = function(e) {
        rep(NA_real_, length(ids))
      })
      
      if (length(ly_vals) == length(tc_vals)) {
        tc_vals <- ifelse(
          !is.na(ly_vals) & ly_vals >= 1 & ly_vals <= target_year - 2000,
          0,
          tc_vals
        )
      }
    }
    
    tc_out[tmp$.row_id[ids]] <- tc_vals
  }
  
  tc_out
}

# ----------------------------
# 4. Load val.rm
# ----------------------------

input_csv <- file.path(readyDir, "agbref_addv7_with_tc.csv")

if (!exists("val.rm")) {
  if (!file.exists(input_csv)) {
    stop(
      "val.rm not found in memory and input CSV missing:\n",
      input_csv,
      call. = FALSE
    )
  }
  val.rm <- read.csv(input_csv, stringsAsFactors = FALSE)
}

if (!"ZONE" %in% names(val.rm)) val.rm$ZONE <- "All"
if (!"BIO" %in% names(val.rm)) stop("val.rm must contain BIO.", call. = FALSE)

val.rm$BIO <- ifelse(is.na(val.rm$BIO), "NA", val.rm$BIO)

# ----------------------------
# 5. Clean required numeric columns
# ----------------------------

message("\nCleaning numeric columns...")

val.rm$POINT_X <- num_clean(val.rm$POINT_X)
val.rm$POINT_Y <- num_clean(val.rm$POINT_Y)

if ("AGB_T_HA" %in% names(val.rm)) {
  val.rm$AGB_T_HA <- num_clean(val.rm$AGB_T_HA)
}

if ("AGB_T_HA_ORIG" %in% names(val.rm)) {
  val.rm$AGB_T_HA_ORIG <- num_clean(val.rm$AGB_T_HA_ORIG)
}

if (!"AGB_T_HA_ORIG" %in% names(val.rm)) {
  if (!"AGB_T_HA" %in% names(val.rm)) {
    stop("Both AGB_T_HA_ORIG and AGB_T_HA are missing.", call. = FALSE)
  }
  val.rm$AGB_T_HA_ORIG <- val.rm$AGB_T_HA
}

if (all(is.na(val.rm$AGB_T_HA_ORIG)) && "AGB_T_HA" %in% names(val.rm)) {
  val.rm$AGB_T_HA_ORIG <- val.rm$AGB_T_HA
}

if (!"AGB_T_HA" %in% names(val.rm) || all(is.na(val.rm$AGB_T_HA))) {
  val.rm$AGB_T_HA <- val.rm$AGB_T_HA_ORIG
}

if ("SIZE_HA" %in% names(val.rm)) {
  val.rm$SIZE_HA <- num_clean(val.rm$SIZE_HA)
} else {
  val.rm$SIZE_HA <- NA_real_
}

if ("varTot" %in% names(val.rm)) {
  val.rm$varTot <- num_clean(val.rm$varTot)
} else {
  val.rm$varTot <- NA_real_
}

if (all(is.na(val.rm$varTot))) {
  val.rm$varTot <- 1
}

if ("tc" %in% names(val.rm)) {
  val.rm$tc <- num_clean(val.rm$tc)
}

# ----------------------------
# 6. Check GFC files
# ----------------------------

tc_files <- list.files(
  gfcDir,
  pattern = "treecover2000.*\\.tif$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

ly_files <- list.files(
  gfcDir,
  pattern = "lossyear.*\\.tif$",
  full.names = TRUE,
  recursive = TRUE,
  ignore.case = TRUE
)

if (length(tc_files) == 0) {
  stop("No treecover2000 tiles found in: ", gfcDir, call. = FALSE)
}

message("Treecover tiles found: ", length(tc_files))
message("Lossyear tiles found:  ", length(ly_files))

# ----------------------------
# 7. Input diagnostics
# ----------------------------

message("\n--- Input diagnostics after cleaning ---")
message("nrow(val.rm): ", nrow(val.rm))

message("POINT_X summary:")
print(summary(val.rm$POINT_X))

message("POINT_Y summary:")
print(summary(val.rm$POINT_Y))

message("AGB_T_HA_ORIG summary:")
print(summary(val.rm$AGB_T_HA_ORIG))

message("AGB_T_HA summary:")
print(summary(val.rm$AGB_T_HA))

message("Valid coord + AGB_ORIG rows:")
print(sum(
  is.finite(val.rm$POINT_X) &
    is.finite(val.rm$POINT_Y) &
    is.finite(val.rm$AGB_T_HA_ORIG)
))

diag_cells <- val.rm %>%
  mutate(
    Xnew_100m = 0.001 * (0.5 + POINT_X %/% 0.001),
    Ynew_100m = 0.001 * (0.5 + POINT_Y %/% 0.001)
  ) %>%
  summarise(
    n_points = n(),
    n_100m_cells = n_distinct(paste(Xnew_100m, Ynew_100m))
  )

print(diag_cells)

# ----------------------------
# 8. Main invDasymetry old-style with exact GFC tiles
# ----------------------------

invDasymetry_gfc_exact <- function(plots,
                                   clmn = "ZONE",
                                   value = "All",
                                   aggr = 0.001,
                                   minPlots = 1,
                                   forestTHs = c(10),
                                   target_year = 2010,
                                   ncores = 1) {
  
  if (value == "All") plots$ZONE <- "All"
  
  if (!clmn %in% names(plots)) {
    stop("Attribute ", clmn, " not found.", call. = FALSE)
  }
  
  plots <- plots[plots[[clmn]] == value, , drop = FALSE]
  
  if (nrow(plots) == 0) {
    stop("There are no records satisfying the selection criterion.", call. = FALSE)
  }
  
  needed <- c(
    "POINT_X", "POINT_Y",
    "AGB_T_HA", "AGB_T_HA_ORIG",
    "SIZE_HA", "BIO", "CODE", "OPEN", "VER",
    "INVENTORY", "TIER", "AVG_YEAR",
    "varTot"
  )
  
  miss <- setdiff(needed, names(plots))
  if (length(miss) > 0) {
    stop("Missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  
  plots$POINT_X <- num_clean(plots$POINT_X)
  plots$POINT_Y <- num_clean(plots$POINT_Y)
  plots$AGB_T_HA <- num_clean(plots$AGB_T_HA)
  plots$AGB_T_HA_ORIG <- num_clean(plots$AGB_T_HA_ORIG)
  plots$SIZE_HA <- num_clean(plots$SIZE_HA)
  plots$varTot <- num_clean(plots$varTot)
  
  if (all(is.na(plots$AGB_T_HA_ORIG)) && any(is.finite(plots$AGB_T_HA))) {
    plots$AGB_T_HA_ORIG <- plots$AGB_T_HA
  }
  
  if (all(is.na(plots$AGB_T_HA)) && any(is.finite(plots$AGB_T_HA_ORIG))) {
    plots$AGB_T_HA <- plots$AGB_T_HA_ORIG
  }
  
  if (all(is.na(plots$varTot))) {
    plots$varTot <- 1
  }
  
  plots <- plots[
    is.finite(plots$POINT_X) &
      is.finite(plots$POINT_Y) &
      is.finite(plots$AGB_T_HA_ORIG),
    ,
    drop = FALSE
  ]
  
  message("Valid plots after coordinate/AGB filter: ", nrow(plots))
  
  if (nrow(plots) == 0) {
    stop("No valid plots after coordinate/AGB filtering.", call. = FALSE)
  }
  
  # point-level tc for this target year
  plots$tc <- extract_tc_points_gfc(
    dat = plots,
    target_year = target_year
  )
  
  plots$tc <- num_clean(plots$tc)
  
  # same old aggregation center formula
  plots$Xnew <- aggr * (0.5 + plots$POINT_X %/% aggr)
  plots$Ynew <- aggr * (0.5 + plots$POINT_Y %/% aggr)
  
  message("Distinct grid cells at aggr ", aggr, ": ",
          dplyr::n_distinct(paste(plots$Xnew, plots$Ynew)))
  
  plots$inv <- ifelse(
    is.finite(plots$varTot) & plots$varTot > 0,
    1 / plots$varTot,
    NA_real_
  )
  
  plots$sdMap <- 1
  
  # ----------------------------
  # old-style aggregation
  # ----------------------------
  
  plotsTMP <- aggregate(
    plots[, c("AGB_T_HA_ORIG", "SIZE_HA")],
    list(plots$Xnew, plots$Ynew),
    mean,
    na.rm = TRUE
  )
  
  modalTMP <- aggregate(
    plots[, c("BIO", "CODE", "OPEN", "VER", "INVENTORY", "TIER", "AVG_YEAR")],
    list(plots$Xnew, plots$Ynew),
    modalClass
  )
  
  plotsTMP <- cbind(plotsTMP, modalTMP[, -c(1, 2), drop = FALSE])
  
  varPlot <- aggregate(
    plots[, "varTot", drop = FALSE],
    list(plots$Xnew, plots$Ynew),
    safe_inv_var
  )
  
  varMap <- aggregate(
    plots[, "sdMap", drop = FALSE],
    list(plots$Xnew, plots$Ynew),
    safe_inv_var
  )
  
  agbW <- plyr::ddply(
    plots,
    .(paste(Ynew, Xnew)),
    function(z) {
      data.frame(
        Xnew = mean(z$Xnew),
        Ynew = mean(z$Ynew),
        AGB_T_HA = safe_wmean(z$AGB_T_HA, z$inv)
      )
    }
  )
  
  tcSD <- aggregate(
    plots$tc,
    list(plots$Xnew, plots$Ynew),
    safe_sd
  )
  
  tcMean <- aggregate(
    plots$tc,
    list(plots$Xnew, plots$Ynew),
    safe_mean
  )
  
  blockCOUNT <- aggregate(
    plots[, "AGB_T_HA", drop = FALSE],
    list(plots$Xnew, plots$Ynew),
    function(x) length(na.omit(x))
  )
  
  names(plotsTMP)[1:2] <- c("POINT_X", "POINT_Y")
  names(varPlot) <- c("POINT_X", "POINT_Y", "varPlot")
  names(varMap) <- c("POINT_X", "POINT_Y", "varMap")
  names(tcSD) <- c("POINT_X", "POINT_Y", "TC_PLT_SD")
  names(tcMean) <- c("POINT_X", "POINT_Y", "TC_PLT_MEAN")
  names(blockCOUNT) <- c("POINT_X", "POINT_Y", "n")
  
  plotsTMP <- plotsTMP %>%
    left_join(varPlot, by = c("POINT_X", "POINT_Y")) %>%
    left_join(varMap, by = c("POINT_X", "POINT_Y")) %>%
    left_join(
      agbW[, c("Xnew", "Ynew", "AGB_T_HA")],
      by = c("POINT_X" = "Xnew", "POINT_Y" = "Ynew")
    ) %>%
    left_join(tcSD, by = c("POINT_X", "POINT_Y")) %>%
    left_join(tcMean, by = c("POINT_X", "POINT_Y")) %>%
    left_join(blockCOUNT, by = c("POINT_X", "POINT_Y")) %>%
    filter(n >= minPlots)
  
  message("Aggregated cells after minPlots: ", nrow(plotsTMP))
  
  if (nrow(plotsTMP) == 0) {
    warning("No cells left after minPlots filtering.")
    return(plotsTMP)
  }
  
  rsl <- aggr
  
  # ----------------------------
  # parallel exact-tile GFC sampling per cell
  # ----------------------------
  
  nworkers <- max(1, min(ncores, nrow(plotsTMP)))
  
  cl <- parallel::makeCluster(nworkers)
  doParallel::registerDoParallel(cl)
  
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
  }, add = TRUE)
  
  FFAGB <- foreach(
    i = seq_len(nrow(plotsTMP)),
    .combine = "rbind",
    .errorhandling = "remove",
    .packages = c("sp", "raster"),
    .export = c(
      "SRS",
      "gfcDir",
      "forestTHs",
      "num_clean",
      "MakeBlockPolygon",
      "gfc_tile_codes_from_pol",
      "find_gfc_tile",
      "TCtileNames_gfc",
      "LYtileNames_gfc",
      "sampleTreeCover_gfc",
      "sampleTCStats_gfc"
    )
  ) %dopar% {
    
    pol <- MakeBlockPolygon(
      plotsTMP$POINT_X[i],
      plotsTMP$POINT_Y[i],
      rsl
    )
    
    gridVals <- sampleTCStats_gfc(
      pol = pol,
      target_year = target_year
    )
    
    treeCovers <- sampleTreeCover_gfc(
      pol = pol,
      thresholds = forestTHs,
      target_year = target_year,
      wghts = FALSE
    )
    
    c(
      plotsTMP$POINT_X[i],
      plotsTMP$POINT_Y[i],
      plotsTMP$TC_PLT_SD[i],
      plotsTMP$TC_PLT_MEAN[i],
      gridVals$SD,
      gridVals$MEAN,
      as.numeric(plotsTMP$n[i]),
      treeCovers[1] * plotsTMP$AGB_T_HA_ORIG[i],
      plotsTMP$SIZE_HA[i],
      plotsTMP$OPEN[i],
      plotsTMP$VER[i],
      plotsTMP$varPlot[i],
      plotsTMP$AVG_YEAR[i],
      plotsTMP$BIO[i],
      plotsTMP$CODE[i],
      plotsTMP$INVENTORY[i],
      plotsTMP$TIER[i]
    )
  }
  
  if (is.null(FFAGB) || nrow(as.data.frame(FFAGB)) == 0) {
    warning("All parallel GFC sampling tasks failed.")
    return(data.frame())
  }
  
  FFAGB <- as.data.frame(FFAGB, stringsAsFactors = FALSE)
  
  names(FFAGB) <- c(
    "POINT_X", "POINT_Y",
    "TC_PLT_SD", "TC_PLT_MEAN",
    "TC_GRID_SD", "TC_GRID_MEAN",
    "n",
    "AGB_T_HA",
    "SIZE_HA",
    "OPEN", "VER",
    "varTot",
    "AVG_YEAR",
    "BIO", "CODE", "INVENTORY", "TIER"
  )
  
  FFAGB[, 1:13] <- lapply(FFAGB[, 1:13], num_clean)
  
  FFAGB
}

# ----------------------------
# 9. Run multi-resolution x multi-epoch
# ----------------------------

data_frames <- list()
run_log <- data.frame()

for (i in seq_len(nrow(scales))) {
  
  label <- scales$label[i]
  aggr <- scales$aggr[i]
  minPlots <- scales$minPlots[i]
  
  for (epoch in target_epochs) {
    
    message("\n================================================")
    message("Building: ", label, " | epoch ", epoch)
    message("Aggregation: ", aggr)
    message("================================================")
    
    mpAGB <- invDasymetry_gfc_exact(
      plots = val.rm,
      clmn = "ZONE",
      value = "All",
      aggr = aggr,
      minPlots = minPlots,
      forestTHs = forestTHs,
      target_year = epoch,
      ncores = ncores
    )
    
    mpAGB[] <- lapply(
      mpAGB,
      function(x) if (is.list(x)) sapply(x, toString) else x
    )
    
    list_name <- paste(label, epoch, sep = "_")
    data_frames[[list_name]] <- mpAGB
    
    save(
      mpAGB,
      file = file.path(
        outDir,
        paste0("AGBref_ready_", label, "_", epoch, "_", Sys.Date(), ".Rdata")
      )
    )
    
    write.csv(
      mpAGB,
      file.path(
        outDir,
        paste0("AGBref_ready_", label, "_", epoch, "_", Sys.Date(), ".csv")
      ),
      row.names = FALSE
    )
    
    run_log <- rbind(
      run_log,
      data.frame(
        resolution = label,
        aggr = aggr,
        epoch = epoch,
        rows = nrow(mpAGB),
        mean_AGB = mean(num_clean(mpAGB$AGB_T_HA), na.rm = TRUE),
        na_AGB = sum(is.na(num_clean(mpAGB$AGB_T_HA))),
        stringsAsFactors = FALSE
      )
    )
    
    message("Rows: ", nrow(mpAGB))
    message("Mean AGB: ", mean(num_clean(mpAGB$AGB_T_HA), na.rm = TRUE))
    message("NA AGB: ", sum(is.na(num_clean(mpAGB$AGB_T_HA))))
    
    gc(verbose = FALSE)
  }
}

# ----------------------------
# 10. Save combined outputs
# ----------------------------

save(
  data_frames,
  file = file.path(
    outDir,
    paste0("AGBrefs_ready_multires_multiepoch_", Sys.Date(), ".Rdata")
  )
)

saveRDS(
  data_frames,
  file.path(
    outDir,
    paste0("AGBrefs_ready_multires_multiepoch_", Sys.Date(), ".rds")
  )
)

write.csv(
  run_log,
  file.path(
    outDir,
    paste0("AGBrefs_ready_multires_multiepoch_log_", Sys.Date(), ".csv")
  ),
  row.names = FALSE
)

message("\nFinished multi-resolution x multi-epoch AGBref-ready aggregation.")
message("Output folder: ", outDir)

print(run_log)
print(vapply(data_frames, nrow, integer(1)))