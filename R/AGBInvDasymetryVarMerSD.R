# ================================================================
# AGBInvDasymetryVarMerSD.R
# Reproducible PVIR map sampling and inverse-variance aggregation
#
# Main function:
#   invDasymetry()
#
# Expected globals passed by callAIDFlex.R / 03_validation.R:
#   resultsFolder
#   agbTilesFolder
#   treeCoverFolder
#   forestTHs
#   pvirMapYear
#   pvirProductVersion
#   pvirMapResolution
#
# Optional globals:
#   sdTilesFolder       # defaults to agbTilesFolder
#   pvirCores           # default detectCores() - 1
# ================================================================

suppressPackageStartupMessages({
  library(sp)
  library(raster)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(dplyr)
  library(plyr)
})

dir.create(resultsFolder, recursive = TRUE, showWarnings = FALSE)

SRS <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs")

# ----------------------------
# 1. Global options
# ----------------------------

get_global <- function(name, default = NULL) {
  if (exists(name, inherits = TRUE)) {
    get(name, inherits = TRUE)
  } else if (exists(name, envir = .GlobalEnv, inherits = FALSE)) {
    get(name, envir = .GlobalEnv)
  } else {
    default
  }
}

get_map_year <- function() {
  
  yr <- get_global("pvirMapYear", NA_integer_)
  
  if (!is.na(yr)) return(as.integer(yr))
  
  txt <- paste(
    get_global("resultsFolder", ""),
    get_global("agbTilesFolder", ""),
    sep = " "
  )
  
  yrs <- regmatches(txt, gregexpr("(19|20)\\d{2}", txt))[[1]]
  
  if (length(yrs) == 0) {
    stop(
      "Cannot infer map year. Set global pvirMapYear or include year in resultsFolder/agbTilesFolder.",
      call. = FALSE
    )
  }
  
  as.integer(tail(yrs, 1))
}

get_product_version <- function() {
  get_global("pvirProductVersion", "fv7.0")
}

get_map_resolution <- function() {
  get_global("pvirMapResolution", "100m")
}

get_sd_tiles_folder <- function() {
  get_global("sdTilesFolder", get_global("agbTilesFolder"))
}

# ----------------------------
# 2. Tile naming helpers
# ----------------------------

tile_origin_10deg <- function(x, y) {
  lon <- 10 * (x %/% 10)
  lat <- 10 * (y %/% 10) + 10
  
  list(
    lon = lon,
    lat = lat,
    LtX = ifelse(lon < 0, "W", "E"),
    LtY = ifelse(lat < 0, "S", "N")
  )
}

agb_tile_code <- function(x, y) {
  z <- tile_origin_10deg(x, y)
  
  paste0(
    z$LtY, sprintf("%02d", abs(z$lat)),
    z$LtX, sprintf("%03d", abs(z$lon))
  )
}

treecover_tile_code <- function(x, y) {
  z <- tile_origin_10deg(x, y)
  
  paste0(
    sprintf("%02d", abs(z$lat)), z$LtY,
    "_",
    sprintf("%03d", abs(z$lon)), z$LtX
  )
}

polygon_tile_points <- function(pol) {
  bb <- unname(sp::bbox(pol))
  expand.grid(x = bb[1, ], y = bb[2, ])
}

# ----------------------------
# 3. Polygon maker
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

# ----------------------------
# 4. Tile name functions
# ----------------------------

TCtileNames <- function(pol) {
  
  crds <- polygon_tile_points(pol)
  
  fnms <- vapply(seq_len(nrow(crds)), function(i) {
    
    tile <- treecover_tile_code(crds$x[i], crds$y[i])
    
    file.path(
      treeCoverFolder,
      paste0(tile, "_treecover2010_v3.tif")
    )
    
  }, character(1))
  
  unique(fnms)
}

AGBtileNames <- function(pol) {
  
  crds <- polygon_tile_points(pol)
  
  yr  <- get_map_year()
  ver <- get_product_version()
  res <- get_map_resolution()
  
  fnms <- vapply(seq_len(nrow(crds)), function(i) {
    
    tile <- agb_tile_code(crds$x[i], crds$y[i])
    
    file.path(
      agbTilesFolder,
      paste0(
        tile,
        "_ESACCI-BIOMASS-L4-AGB-MERGED-",
        res,
        "-",
        yr,
        "-",
        ver,
        ".tif"
      )
    )
    
  }, character(1))
  
  unique(fnms)
}

SDtileNames <- function(pol) {
  
  crds <- polygon_tile_points(pol)
  
  yr  <- get_map_year()
  ver <- get_product_version()
  res <- get_map_resolution()
  sdir <- get_sd_tiles_folder()
  
  fnms <- vapply(seq_len(nrow(crds)), function(i) {
    
    tile <- agb_tile_code(crds$x[i], crds$y[i])
    
    file.path(
      sdir,
      paste0(
        tile,
        "_ESACCI-BIOMASS-L4-AGB_SD-MERGED-",
        res,
        "-",
        yr,
        "-",
        ver,
        ".tif"
      )
    )
    
  }, character(1))
  
  unique(fnms)
}

# ----------------------------
# 5. Raster sampling helpers
# ----------------------------

weighted_extract_mean <- function(files, pol) {
  
  files <- files[file.exists(files)]
  
  if (length(files) == 0) return(0)
  
  vls <- matrix(ncol = 2, nrow = 0)
  
  for (f in files) {
    
    r <- raster::raster(f)
    
    ext <- try(
      raster::extract(
        r,
        pol,
        weights = TRUE,
        normalizeWeights = FALSE
      )[[1]],
      silent = TRUE
    )
    
    if (!inherits(ext, "try-error") && !is.null(ext)) {
      vls <- rbind(vls, ext)
    }
  }
  
  if (nrow(vls) == 0) return(0)
  
  vls <- vls[!is.na(vls[, 1]), , drop = FALSE]
  
  if (nrow(vls) == 0) return(0)
  
  sum(vls[, 1] * vls[, 2], na.rm = TRUE) / sum(vls[, 2], na.rm = TRUE)
}

unweighted_extract_mean <- function(files, pol) {
  
  files <- files[file.exists(files)]
  
  if (length(files) == 0) return(0)
  
  vls <- numeric()
  
  for (f in files) {
    
    r <- raster::raster(f)
    
    ext <- try(
      raster::extract(r, pol)[[1]],
      silent = TRUE
    )
    
    if (!inherits(ext, "try-error") && !is.null(ext)) {
      vls <- c(vls, ext)
    }
  }
  
  if (length(stats::na.omit(vls)) == 0) return(0)
  
  mean(vls, na.rm = TRUE)
}

sampleTreeCover <- function(pol, thresholds, wghts = FALSE) {
  
  files <- TCtileNames(pol)
  files <- files[file.exists(files)]
  
  if (length(files) == 0) {
    return(rep(0, length(thresholds)))
  }
  
  out <- numeric()
  
  if (wghts) {
    
    vls <- matrix(ncol = 2, nrow = 0)
    
    for (f in files) {
      
      r <- raster::raster(f)
      
      ext <- try(
        raster::extract(
          r,
          pol,
          weights = TRUE,
          normalizeWeights = FALSE
        )[[1]],
        silent = TRUE
      )
      
      if (!inherits(ext, "try-error") && !is.null(ext)) {
        vls <- rbind(vls, ext)
      }
    }
    
    if (nrow(vls) == 0) {
      return(rep(0, length(thresholds)))
    }
    
    vls <- vls[!is.na(vls[, 1]), , drop = FALSE]
    
    if (nrow(vls) == 0) {
      return(rep(0, length(thresholds)))
    }
    
    for (threshold in thresholds) {
      tmp <- vls
      tmp[, 1] <- ifelse(tmp[, 1] > threshold, 1, 0)
      out <- c(
        out,
        sum(tmp[, 1] * tmp[, 2], na.rm = TRUE) /
          sum(tmp[, 2], na.rm = TRUE)
      )
    }
    
  } else {
    
    vls <- numeric()
    
    for (f in files) {
      r <- raster::raster(f)
      ext <- try(raster::extract(r, pol)[[1]], silent = TRUE)
      
      if (!inherits(ext, "try-error") && !is.null(ext)) {
        vls <- c(vls, ext)
      }
    }
    
    for (threshold in thresholds) {
      out <- c(out, mean(ifelse(vls > threshold, 1, 0), na.rm = TRUE))
    }
  }
  
  out
}

sampleAGBmap <- function(pol, wghts = FALSE) {
  
  files <- AGBtileNames(pol)
  
  if (wghts) {
    weighted_extract_mean(files, pol)
  } else {
    unweighted_extract_mean(files, pol)
  }
}

sampleSDmap <- function(pol, wghts = FALSE) {
  
  files <- SDtileNames(pol)
  
  if (wghts) {
    weighted_extract_mean(files, pol)
  } else {
    unweighted_extract_mean(files, pol)
  }
}

# ----------------------------
# 6. Aggregation helpers
# ----------------------------

modalClass <- function(x) {
  
  x <- x[!is.na(x)]
  
  if (length(x) == 0) return(NA)
  
  y <- table(x)
  names(y)[which.max(y)]
}

first_existing_raster_resolution <- function() {
  
  files <- list.files(
    agbTilesFolder,
    pattern = "\\.tif$",
    full.names = TRUE,
    recursive = FALSE
  )
  
  files <- files[!grepl("AGB_SD|1000m|aux", files, ignore.case = TRUE)]
  
  if (length(files) == 0) {
    stop("No AGB raster files found in agbTilesFolder: ", agbTilesFolder, call. = FALSE)
  }
  
  raster::xres(raster::raster(files[1]))
}

ensure_required_columns <- function(plots) {
  
  required <- c(
    "POINT_X", "POINT_Y",
    "AGB_T_HA", "AGB_T_HA_ORIG",
    "SIZE_HA", "BIO", "CODE",
    "OPEN", "VER", "INVENTORY", "TIER",
    "varTot", "sdMap"
  )
  
  missing <- setdiff(required, names(plots))
  
  if (length(missing) > 0) {
    stop(
      "plotFile is missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

safe_inverse_variance <- function(x) {
  
  x <- as.numeric(x)
  x[is.na(x) | x <= 0] <- NA_real_
  
  if (all(is.na(x))) return(NA_real_)
  
  1 / sum(1 / x, na.rm = TRUE)
}

# ----------------------------
# 7. Main function
# ----------------------------

invDasymetry <- function(clmn = "ZONE",
                         value = "Australia",
                         aggr = NULL,
                         plotFile = plots,
                         minPlots = 1,
                         wghts = FALSE) {
  
  plots <- plotFile
  
  ensure_required_columns(plots)
  
  if (is.null(aggr)) minPlots <- 1
  
  if (!clmn %in% names(plots)) {
    stop("Attribute not found: ", clmn, call. = FALSE)
  }
  
  plots <- plots[plots[[clmn]] == value & !is.na(plots[[clmn]]), , drop = FALSE]
  
  if (nrow(plots) == 0) {
    stop("There are no records satisfying the selection criterion.", call. = FALSE)
  }
  
  plots$CODE <- as.character(plots$CODE)
  
  # ----------------------------
  # 7.1 Aggregate reference data if requested
  # ----------------------------
  
  if (!is.null(aggr)) {
    
    plots$Xnew <- aggr * (0.5 + plots$POINT_X %/% aggr)
    plots$Ynew <- aggr * (0.5 + plots$POINT_Y %/% aggr)
    
    plots$varTot <- as.numeric(plots$varTot)
    plots$varTot[is.na(plots$varTot) | plots$varTot <= 0] <- NA_real_
    plots$inv <- 1 / plots$varTot
    
    block_count <- aggregate(
      plots$AGB_T_HA,
      by = list(Xnew = plots$Xnew, Ynew = plots$Ynew),
      FUN = function(x) length(stats::na.omit(x))
    )
    
    names(block_count)[3] <- "n"
    
    weighted_agb <- plyr::ddply(
      plots,
      .variables = c("Xnew", "Ynew"),
      .fun = function(z) {
        data.frame(
          AGB_T_HA = stats::weighted.mean(z$AGB_T_HA, z$inv, na.rm = TRUE)
        )
      }
    )
    
    plotsTMP <- plots %>%
      dplyr::group_by(.data$Xnew, .data$Ynew) %>%
      dplyr::summarise(
        AGB_T_HA_ORIG = mean(.data$AGB_T_HA_ORIG, na.rm = TRUE),
        SIZE_HA = mean(.data$SIZE_HA, na.rm = TRUE),
        BIO = modalClass(.data$BIO),
        CODE = modalClass(.data$CODE),
        OPEN = modalClass(.data$OPEN),
        VER = modalClass(.data$VER),
        INVENTORY = modalClass(.data$INVENTORY),
        TIER = modalClass(.data$TIER),
        varPlot = safe_inverse_variance(.data$varTot),
        varMap = safe_inverse_variance(.data$sdMap^2),
        .groups = "drop"
      ) %>%
      dplyr::left_join(weighted_agb, by = c("Xnew", "Ynew")) %>%
      dplyr::left_join(block_count, by = c("Xnew", "Ynew")) %>%
      dplyr::filter(.data$n >= minPlots)
    
    if (nrow(plotsTMP) < 2) {
      
      warning(
        "Fewer than 2 aggregated cells meet minPlots = ",
        minPlots,
        ". Keeping first available aggregated cells for compatibility."
      )
      
      plotsTMP <- plots %>%
        dplyr::group_by(.data$Xnew, .data$Ynew) %>%
        dplyr::summarise(
          AGB_T_HA_ORIG = mean(.data$AGB_T_HA_ORIG, na.rm = TRUE),
          SIZE_HA = mean(.data$SIZE_HA, na.rm = TRUE),
          BIO = modalClass(.data$BIO),
          CODE = modalClass(.data$CODE),
          OPEN = modalClass(.data$OPEN),
          VER = modalClass(.data$VER),
          INVENTORY = modalClass(.data$INVENTORY),
          TIER = modalClass(.data$TIER),
          varPlot = safe_inverse_variance(.data$varTot),
          varMap = safe_inverse_variance(.data$sdMap^2),
          n = length(stats::na.omit(.data$AGB_T_HA)),
          .groups = "drop"
        ) %>%
        dplyr::left_join(weighted_agb, by = c("Xnew", "Ynew")) %>%
        utils::head(2)
    }
    
    plots <- plotsTMP %>%
      dplyr::transmute(
        POINT_X = .data$Xnew,
        POINT_Y = .data$Ynew,
        AGB_T_HA_ORIG,
        SIZE_HA,
        BIO,
        CODE,
        OPEN,
        VER,
        INVENTORY,
        TIER,
        varPlot,
        varMap,
        AGB_T_HA,
        n
      )
    
    rsl <- aggr
    
  } else {
    
    rsl <- first_existing_raster_resolution()
    plots$n <- 1
    plots$varPlot <- plots$varTot
    plots$varMap <- plots$sdMap^2
  }
  
  # ----------------------------
  # 8. Sample forest fraction, map AGB, and map SD
  # ----------------------------
  
  ncores <- get_global("pvirCores", max(1, parallel::detectCores() - 1))
  ncores <- min(ncores, max(1, nrow(plots)))
  
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  
  on.exit({
    try(parallel::stopCluster(cl), silent = TRUE)
  }, add = TRUE)
  
  FFAGB <- foreach::foreach(
    i = seq_len(nrow(plots)),
    .combine = "rbind",
    .packages = c("sp", "raster"),
    .export = c(
      "MakeBlockPolygon",
      "sampleTreeCover",
      "TCtileNames",
      "AGBtileNames",
      "SDtileNames",
      "sampleAGBmap",
      "sampleSDmap",
      "weighted_extract_mean",
      "unweighted_extract_mean",
      "tile_origin_10deg",
      "agb_tile_code",
      "treecover_tile_code",
      "polygon_tile_points",
      "get_global",
      "get_map_year",
      "get_product_version",
      "get_map_resolution",
      "get_sd_tiles_folder",
      "SRS",
      "agbTilesFolder",
      "treeCoverFolder",
      "forestTHs",
      "resultsFolder",
      "pvirMapYear",
      "pvirProductVersion",
      "pvirMapResolution",
      "sdTilesFolder"
    )
  ) %dopar% {
    
    pol <- MakeBlockPolygon(
      plots$POINT_X[i],
      plots$POINT_Y[i],
      rsl
    )
    
    if (is.null(aggr)) {
      
      if (is.na(plots$SIZE_HA[i])) {
        treeCovers <- sampleTreeCover(pol, forestTHs, wghts)
      } else if (plots$SIZE_HA[i] >= 1) {
        treeCovers <- rep(1, length(forestTHs))
      } else {
        treeCovers <- sampleTreeCover(pol, forestTHs, wghts)
      }
      
    } else {
      
      treeCovers <- sampleTreeCover(pol, forestTHs, wghts)
    }
    
    wghts2 <- ifelse(is.null(aggr), FALSE, wghts)
    
    c(
      treeCovers * plots$AGB_T_HA[i],
      treeCovers * plots$AGB_T_HA_ORIG[i],
      sampleAGBmap(pol, wghts2),
      plots$POINT_X[i],
      plots$POINT_Y[i],
      sqrt(plots$varPlot[i]),
      sampleSDmap(pol, wghts2),
      plots$SIZE_HA[i],
      plots$n[i],
      as.character(plots$BIO[i]),
      as.character(plots$CODE[i]),
      as.character(plots$OPEN[i]),
      as.character(plots$VER[i]),
      as.character(plots$INVENTORY[i]),
      as.character(plots$TIER[i])
    )
  }
  
  FFAGB <- as.data.frame(FFAGB, stringsAsFactors = FALSE)
  
  n_num <- ncol(FFAGB) - 6
  
  FFAGB_numeric <- data.frame(
    lapply(FFAGB[, seq_len(n_num), drop = FALSE], as.numeric)
  )
  
  FFAGB_char <- FFAGB[, (n_num + 1):ncol(FFAGB), drop = FALSE]
  
  FFAGB <- cbind(FFAGB_numeric, FFAGB_char)
  
  names(FFAGB) <- c(
    "plotAGB_10",
    "orgPlotAGB",
    "mapAGB",
    "x",
    "y",
    "sdPlot",
    "sdMap",
    "SIZE_HA",
    "n",
    "BIO",
    "CODE",
    "OPEN",
    "VER",
    "INVENTORY",
    "TIER"
  )
  
  FFAGB
}