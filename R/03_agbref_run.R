# ================================================================
# 02_agbref_ready.R
# Build AGBref-ready data across resolutions and epochs
# WITHOUT temporal biomass adjustment
#
# Inputs expected:
#   - val.rm in memory OR:
#     outputs/pvir_ready_demo/agbref_addv7_with_tc.csv
#
# GFC expected:
#   data/GFC/
#     Hansen_GFC-...treecover2000_*.tif
#     Hansen_GFC-...lossyear_*.tif
#
# Outputs:
#   outputs/agbref_ready/
#     AGBref_ready_<resolution>_<epoch>_noTempAdj_<date>.Rdata/.csv
#     AGBrefs_ready_all_resolutions_epochs_noTempAdj_<date>.Rdata/.rds
# ================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(foreach)
  library(doParallel)
})

# ----------------------------
# 0. Paths and settings
# ----------------------------

projectDir <- getwd()
dataDir    <- file.path(projectDir, "data")
readyDir   <- file.path(projectDir, "outputs", "pvir_ready_demo")
outDir     <- file.path(projectDir, "outputs", "agbref_ready")
gfcDir     <- file.path(projectDir, "data", "GFC")

dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

forestTHs <- c(10)

# Safer for terra/GFC VRT extraction. Increase only after successful test.
ncores <- 1
# ncores <- max(1, parallel::detectCores() - 2)

target_epochs <- c(2005, 2010, 2015, 2020)

# "epoch" = GFC tree cover adjusted by lossyear for each epoch
# "baseline_2010" = always use 2010 tree-cover approximation
tc_year_mode <- "epoch"
baseline_tc_year <- 2010

message("Project root: ", projectDir)
message("Output dir:   ", outDir)
message("GFC dir:      ", gfcDir)

# ----------------------------
# 1. Load val.rm if needed
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

if (!"BIO" %in% names(val.rm)) {
  stop("val.rm must already contain BIO before running this script.", call. = FALSE)
}

val.rm$BIO <- ifelse(is.na(val.rm$BIO), "NA", val.rm$BIO)

# ----------------------------
# 2. Build GFC VRTs
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

tc_vrt <- file.path(gfcDir, "treecover2000.vrt")
terra::vrt(tc_files, filename = tc_vrt, overwrite = TRUE)

if (length(ly_files) > 0) {
  ly_vrt <- file.path(gfcDir, "lossyear.vrt")
  terra::vrt(ly_files, filename = ly_vrt, overwrite = TRUE)
} else {
  ly_vrt <- NA_character_
}

message("Treecover VRT: ", tc_vrt)
message("Lossyear VRT:  ", ly_vrt)

# ----------------------------
# 3. Helper functions
# ----------------------------

modalClass <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  names(sort(table(x), decreasing = TRUE))[1]
}

safe_inv_var <- function(x) {
  x <- as.numeric(x)
  x[is.na(x) | x <= 0] <- NA_real_
  if (all(is.na(x))) return(NA_real_)
  1 / sum(1 / x, na.rm = TRUE)
}

safe_wmean <- function(x, w) {
  x <- as.numeric(x)
  w <- as.numeric(w)
  ok <- is.finite(x) & is.finite(w) & w > 0
  
  if (sum(ok) == 0) return(mean(x, na.rm = TRUE))
  
  weighted.mean(x[ok], w[ok], na.rm = TRUE)
}

make_block_polygon <- function(x, y, size) {
  
  xmin <- x - size / 2
  xmax <- x + size / 2
  ymin <- y - size / 2
  ymax <- y + size / 2
  
  coords <- matrix(
    c(
      xmin, ymin,
      xmax, ymin,
      xmax, ymax,
      xmin, ymax,
      xmin, ymin
    ),
    ncol = 2,
    byrow = TRUE
  )
  
  terra::vect(coords, type = "polygons", crs = "EPSG:4326")
}

extract_gfc_values <- function(pol,
                               tc_vrt,
                               ly_vrt = NA_character_,
                               target_year = 2010) {
  
  tc_r <- terra::rast(tc_vrt)
  tc <- terra::extract(tc_r, pol)[, 2]
  
  if (!is.na(ly_vrt) && file.exists(ly_vrt) && target_year > 2000) {
    
    ly_r <- terra::rast(ly_vrt)
    ly <- terra::extract(ly_r, pol)[, 2]
    
    tc <- ifelse(
      !is.na(ly) & ly >= 1 & ly <= target_year - 2000,
      0,
      tc
    )
  }
  
  tc
}

sample_tc_block <- function(pol,
                            tc_vrt,
                            ly_vrt = NA_character_,
                            forestTHs = c(10),
                            target_year = 2010) {
  
  v <- extract_gfc_values(
    pol = pol,
    tc_vrt = tc_vrt,
    ly_vrt = ly_vrt,
    target_year = target_year
  )
  
  v <- v[!is.na(v)]
  
  if (length(v) == 0) {
    return(list(
      TC_GRID_SD = NA_real_,
      TC_GRID_MEAN = NA_real_,
      FF = rep(NA_real_, length(forestTHs))
    ))
  }
  
  list(
    TC_GRID_SD = sd(v, na.rm = TRUE),
    TC_GRID_MEAN = mean(v, na.rm = TRUE),
    FF = sapply(forestTHs, function(th) mean(v > th, na.rm = TRUE))
  )
}

ensure_tc_points <- function(dat,
                             tc_vrt,
                             ly_vrt = NA_character_,
                             target_year = 2010) {
  
  if ("tc" %in% names(dat) && any(!is.na(dat$tc))) {
    return(dat)
  }
  
  message("Extracting point-level GFC tree cover to val.rm$tc ...")
  
  pts <- terra::vect(
    dat,
    geom = c("POINT_X", "POINT_Y"),
    crs = "EPSG:4326"
  )
  
  tc_r <- terra::rast(tc_vrt)
  tc <- terra::extract(tc_r, pts)[, 2]
  
  if (!is.na(ly_vrt) && file.exists(ly_vrt) && target_year > 2000) {
    
    ly_r <- terra::rast(ly_vrt)
    ly <- terra::extract(ly_r, pts)[, 2]
    
    tc <- ifelse(
      !is.na(ly) & ly >= 1 & ly <= target_year - 2000,
      0,
      tc
    )
  }
  
  dat$tc <- tc
  dat
}

# ----------------------------
# 4. Main aggregation function
# ----------------------------

make_agbref_ready <- function(dat,
                              aggr,
                              minPlots = 1,
                              forestTHs = c(10),
                              target_year = 2010,
                              tc_vrt,
                              ly_vrt = NA_character_,
                              ncores = 1) {
  
  dat <- ensure_tc_points(
    dat = dat,
    tc_vrt = tc_vrt,
    ly_vrt = ly_vrt,
    target_year = target_year
  )
  
  if (!"ZONE" %in% names(dat)) dat$ZONE <- "All"
  dat$ZONE <- "All"
  
  needed <- c(
    "POINT_X", "POINT_Y",
    "AGB_T_HA", "AGB_T_HA_ORIG",
    "SIZE_HA", "BIO", "CODE",
    "OPEN", "VER", "INVENTORY", "TIER",
    "AVG_YEAR", "varTot", "tc"
  )
  
  missing <- setdiff(needed, names(dat))
  
  if (length(missing) > 0) {
    stop(
      "Missing required column(s): ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  
  dat <- dat %>%
    filter(
      is.finite(POINT_X),
      is.finite(POINT_Y),
      is.finite(AGB_T_HA)
    ) %>%
    mutate(
      Xnew = aggr * (0.5 + POINT_X %/% aggr),
      Ynew = aggr * (0.5 + POINT_Y %/% aggr),
      varTot = as.numeric(varTot),
      inv = ifelse(is.finite(varTot) & varTot > 0, 1 / varTot, NA_real_)
    )
  
  block_count <- dat %>%
    count(Xnew, Ynew, name = "n")
  
  agb_w <- dat %>%
    group_by(Xnew, Ynew) %>%
    summarise(
      AGB_T_HA = safe_wmean(AGB_T_HA, inv),
      .groups = "drop"
    )
  
  agg_df <- dat %>%
    group_by(Xnew, Ynew) %>%
    summarise(
      AGB_T_HA_ORIG = mean(AGB_T_HA_ORIG, na.rm = TRUE),
      SIZE_HA = mean(SIZE_HA, na.rm = TRUE),
      BIO = modalClass(BIO),
      CODE = modalClass(CODE),
      OPEN = modalClass(OPEN),
      VER = modalClass(VER),
      INVENTORY = modalClass(INVENTORY),
      TIER = modalClass(TIER),
      AVG_YEAR = modalClass(AVG_YEAR),
      varTot = safe_inv_var(varTot),
      TC_PLT_SD = sd(tc, na.rm = TRUE),
      TC_PLT_MEAN = mean(tc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(agb_w, by = c("Xnew", "Ynew")) %>%
    left_join(block_count, by = c("Xnew", "Ynew")) %>%
    filter(n >= minPlots)
  
  if (nrow(agg_df) == 0) {
    warning("No cells left after minPlots filtering for aggr = ", aggr)
    return(agg_df)
  }
  
  if (ncores <= 1) {
    
    sampled <- lapply(seq_len(nrow(agg_df)), function(i) {
      
      pol <- make_block_polygon(agg_df$Xnew[i], agg_df$Ynew[i], aggr)
      
      s <- sample_tc_block(
        pol = pol,
        tc_vrt = tc_vrt,
        ly_vrt = ly_vrt,
        forestTHs = forestTHs,
        target_year = target_year
      )
      
      data.frame(
        TC_GRID_SD = s$TC_GRID_SD,
        TC_GRID_MEAN = s$TC_GRID_MEAN,
        FF = s$FF[1]
      )
    })
    
    sampled <- dplyr::bind_rows(sampled)
    
  } else {
    
    cl <- parallel::makeCluster(min(ncores, nrow(agg_df)))
    doParallel::registerDoParallel(cl)
    
    on.exit({
      try(parallel::stopCluster(cl), silent = TRUE)
    }, add = TRUE)
    
    sampled <- foreach(
      i = seq_len(nrow(agg_df)),
      .combine = "rbind",
      .packages = c("terra"),
      .export = c(
        "make_block_polygon",
        "sample_tc_block",
        "extract_gfc_values"
      )
    ) %dopar% {
      
      pol <- make_block_polygon(agg_df$Xnew[i], agg_df$Ynew[i], aggr)
      
      s <- sample_tc_block(
        pol = pol,
        tc_vrt = tc_vrt,
        ly_vrt = ly_vrt,
        forestTHs = forestTHs,
        target_year = target_year
      )
      
      data.frame(
        TC_GRID_SD = s$TC_GRID_SD,
        TC_GRID_MEAN = s$TC_GRID_MEAN,
        FF = s$FF[1]
      )
    }
    
    sampled <- as.data.frame(sampled)
  }
  
  out <- bind_cols(agg_df, sampled) %>%
    transmute(
      POINT_X = Xnew,
      POINT_Y = Ynew,
      TC_PLT_SD,
      TC_PLT_MEAN,
      TC_GRID_SD,
      TC_GRID_MEAN,
      n,
      AGB_T_HA = FF * AGB_T_HA_ORIG,
      SIZE_HA,
      OPEN,
      VER,
      varTot,
      AVG_YEAR,
      BIO,
      CODE,
      INVENTORY,
      TIER,
      MapYear = as.integer(target_year),
      TC_YEAR = as.integer(target_year),
      TEMPORAL_ADJUSTED = FALSE
    )
  
  out[] <- lapply(
    out,
    function(x) if (is.list(x)) sapply(x, toString) else x
  )
  
  out
}

# ----------------------------
# 5. Run resolutions × epochs WITHOUT temporal adjustment
# ----------------------------

scales <- data.frame(
  label = c("100m", "500m", "1km", "10km", "25km"),
  aggr  = c(0.001, 0.005, 0.01, 0.1, 0.25),
  minPlots = c(1, 1, 1, 1, 1),
  stringsAsFactors = FALSE
)

data_frames <- list()
run_log <- data.frame()

for (i in seq_len(nrow(scales))) {
  
  label <- scales$label[i]
  aggr <- scales$aggr[i]
  minPlots <- scales$minPlots[i]
  
  for (epoch in target_epochs) {
    
    tc_year_used <- ifelse(tc_year_mode == "epoch", epoch, baseline_tc_year)
    
    message("\n================================================")
    message("Building AGBref-ready output: ", label, " | epoch ", epoch)
    message("Aggregation: ", aggr, " degrees")
    message("GFC tree-cover target year: ", tc_year_used)
    message("NO temporal biomass adjustment applied")
    message("================================================")
    
    mpAGB <- make_agbref_ready(
      dat = val.rm,
      aggr = aggr,
      minPlots = minPlots,
      forestTHs = forestTHs,
      target_year = tc_year_used,
      tc_vrt = tc_vrt,
      ly_vrt = ly_vrt,
      ncores = ncores
    )
    
    if (nrow(mpAGB) > 0) {
      mpAGB$MapYear <- as.integer(epoch)
      mpAGB$TC_YEAR <- as.integer(tc_year_used)
      mpAGB$TEMPORAL_ADJUSTED <- FALSE
    }
    
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
        paste0("AGBref_ready_", label, "_", epoch, "_noTempAdj_", Sys.Date(), ".Rdata")
      )
    )
    
    write.csv(
      mpAGB,
      file.path(
        outDir,
        paste0("AGBref_ready_", label, "_", epoch, "_noTempAdj_", Sys.Date(), ".csv")
      ),
      row.names = FALSE
    )
    
    run_log <- rbind(
      run_log,
      data.frame(
        resolution = label,
        aggr = aggr,
        minPlots = minPlots,
        epoch = epoch,
        tc_year_mode = tc_year_mode,
        tc_year_used = tc_year_used,
        temporal_adjusted = FALSE,
        rows = nrow(mpAGB),
        stringsAsFactors = FALSE
      )
    )
    
    message("Saved ", label, " | ", epoch, " rows: ", nrow(mpAGB))
  }
}

# ----------------------------
# 6. Save combined outputs
# ----------------------------

save(
  data_frames,
  file = file.path(
    outDir,
    paste0("AGBrefs_ready_all_resolutions_epochs_noTempAdj_", Sys.Date(), ".Rdata")
  )
)

saveRDS(
  data_frames,
  file.path(
    outDir,
    paste0("AGBrefs_ready_all_resolutions_epochs_noTempAdj_", Sys.Date(), ".rds")
  )
)

write.csv(
  run_log,
  file.path(outDir, paste0("AGBrefs_ready_run_log_noTempAdj_", Sys.Date(), ".csv")),
  row.names = FALSE
)

message("Output folder: ", outDir)

print(run_log)
print(vapply(data_frames, nrow, integer(1)))