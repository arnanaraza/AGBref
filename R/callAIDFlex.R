# ================================================================
# callAIDFlex.R
# Flexible wrapper for PVIR validation runs
#
# Main dependency:
#   R/AGBInvDasymetryVarMerSD.R
#
# Other dependencies:
#   R/Acc.R
#   R/OnePlot.R
#
# Required from 03_validation.R:
#   resultsFolder
#   agbTilesFolder
#   treeCoverFolder
#   forestTHs
#
# Notes:
#   - t0 is ONLY for LIDAR / CoFor / EMAP-type reference products
#   - t1/t2/t3 are plot-size based
#   - skips bins with insufficient data instead of crashing
#   - output filenames and plot captions include epoch, tier, group, scale
# ================================================================

callAID <- function(df,
                    group_var = c("all", "biome", "slope"),
                    scale = c("Agg", "nonAgg"),
                    minPlt = 5,
                    blockRes = 0.1,
                    resultsFolder,
                    rDir,
                    agbTilesFolder,
                    treeCoverFolder,
                    forestTHs = c(10),
                    exclude_bins = NULL,
                    result_prefix = NULL,
                    accuracy_digits = 8) {
  
  group_var <- match.arg(group_var)
  scale <- match.arg(scale)
  
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required.", call. = FALSE)
  }
  
  if (!requireNamespace("stringr", quietly = TRUE)) {
    stop("Package 'stringr' is required.", call. = FALSE)
  }
  
  if (!dir.exists(resultsFolder)) {
    dir.create(resultsFolder, recursive = TRUE, showWarnings = FALSE)
  }
  
  required_scripts <- c(
    "AGBInvDasymetryVarMerSD.R",
    "Acc.R",
    "OnePlot.R"
  )
  
  missing_scripts <- file.path(rDir, required_scripts)[
    !file.exists(file.path(rDir, required_scripts))
  ]
  
  if (length(missing_scripts) > 0) {
    stop(
      "Missing required script(s):\n",
      paste(missing_scripts, collapse = "\n"),
      call. = FALSE
    )
  }
  
  if (!dir.exists(agbTilesFolder)) {
    stop("Missing AGB tiles folder: ", agbTilesFolder, call. = FALSE)
  }
  
  if (!dir.exists(treeCoverFolder)) {
    stop("Missing tree-cover folder: ", treeCoverFolder, call. = FALSE)
  }
  
  # ----------------------------------------------------------------
  # Legacy dependency support
  # AGBInvDasymetryVarMerSD.R expects these globals.
  # Keep them explicit here for reproducibility.
  # ----------------------------------------------------------------
  
  assign("resultsFolder", resultsFolder, envir = .GlobalEnv)
  assign("agbTilesFolder", agbTilesFolder, envir = .GlobalEnv)
  assign("treeCoverFolder", treeCoverFolder, envir = .GlobalEnv)
  assign("forestTHs", forestTHs, envir = .GlobalEnv)
  assign("plotsFolder", dirname(rDir), envir = .GlobalEnv)
  
  source(file.path(rDir, "AGBInvDasymetryVarMerSD.R"), local = .GlobalEnv)
  source(file.path(rDir, "Acc.R"), local = .GlobalEnv)
  source(file.path(rDir, "OnePlot.R"), local = .GlobalEnv)
  
  if (!exists("invDasymetry", mode = "function")) {
    stop("Function invDasymetry() not found after sourcing AGBInvDasymetryVarMerSD.R.", call. = FALSE)
  }
  
  if (!exists("Accuracy", mode = "function")) {
    stop("Function Accuracy() not found after sourcing Acc.R.", call. = FALSE)
  }
  
  if (!exists("OnePlot", mode = "function")) {
    stop("Function OnePlot() not found after sourcing OnePlot.R.", call. = FALSE)
  }
  
  # ----------------------------------------------------------------
  # Helper functions
  # ----------------------------------------------------------------
  
  clean_name <- function(x) {
    x <- as.character(x)
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    x <- gsub("_+", "_", x)
    x <- gsub("^_|_$", "", x)
    x
  }
  
  normalize_tier <- function(x) {
    x <- as.character(x)
    x <- tolower(x)
    x <- gsub("tier", "t", x)
    x <- ifelse(x %in% c("0", "t0"), "t0", x)
    x <- ifelse(x %in% c("1", "t1"), "t1", x)
    x <- ifelse(x %in% c("2", "t2"), "t2", x)
    x <- ifelse(x %in% c("3", "t3"), "t3", x)
    x
  }
  
  is_special_reference <- function(dat) {
    
    code_vals <- if ("CODE" %in% names(dat)) as.character(dat$CODE) else NA_character_
    inv_vals  <- if ("INVENTORY" %in% names(dat)) as.character(dat$INVENTORY) else NA_character_
    type_vals <- if ("DATA_TYPE" %in% names(dat)) as.character(dat$DATA_TYPE) else NA_character_
    
    any(
      code_vals %in% c("LIDAR", "AFR_COF", "EMAP"),
      grepl("_LDR$|_LIDAR$|LIDAR|COF|COFOR|EMAP", code_vals, ignore.case = TRUE),
      inv_vals %in% c("LIDAR"),
      type_vals %in% c("LIDAR", "COF", "COFOR", "EMAP"),
      na.rm = TRUE
    )
  }
  
  get_tier_label <- function(dat) {
    
    if (is_special_reference(dat)) return("t0")
    
    if ("TIER" %in% names(dat)) {
      vals <- unique(normalize_tier(dat$TIER))
      vals <- vals[!is.na(vals) & vals != ""]
      vals <- vals[vals != "t0"]
      
      if (length(vals) == 1) return(vals)
      if (length(vals) > 1) return("tMixed")
    }
    
    if ("SIZE_HA" %in% names(dat)) {
      vals <- dplyr::case_when(
        dat$SIZE_HA < 0.6 ~ "t1",
        dat$SIZE_HA >= 0.6 & dat$SIZE_HA < 3 ~ "t2",
        dat$SIZE_HA >= 3 ~ "t3",
        TRUE ~ NA_character_
      )
      
      vals <- unique(vals[!is.na(vals)])
      
      if (length(vals) == 1) return(vals)
      if (length(vals) > 1) return("tMixed")
    }
    
    "tUnknown"
  }
  
  get_bins <- function(dat, group_var) {
    
    if (group_var == "all") {
      dat$all <- "all"
      return(list(dat = dat, feature_col = "all", bins = "all"))
    }
    
    if (group_var == "biome") {
      if (!"BIO" %in% names(dat)) {
        stop("BIO column is required for group_var = 'biome'.", call. = FALSE)
      }
      
      bins <- sort(unique(dat$BIO))
      bins <- bins[!is.na(bins) & bins != ""]
      return(list(dat = dat, feature_col = "BIO", bins = bins))
    }
    
    if (group_var == "slope") {
      if (!"slp_grp" %in% names(dat)) {
        stop("slp_grp column is required for group_var = 'slope'.", call. = FALSE)
      }
      
      bins <- sort(unique(dat$slp_grp))
      bins <- bins[!is.na(bins) & bins != ""]
      return(list(dat = dat, feature_col = "slp_grp", bins = bins))
    }
  }
  
  get_epoch_label <- function() {
    
    if (exists("pvirMapYear", envir = .GlobalEnv, inherits = FALSE)) {
      return(as.character(get("pvirMapYear", envir = .GlobalEnv)))
    }
    
    yr <- stringr::str_extract(resultsFolder, "(19|20)\\d{2}")
    
    if (is.na(yr)) {
      yr <- "unknownYear"
    }
    
    yr
  }
  
  make_file_stub <- function(prefix, scale, blockRes, minPlt, group_var, bin) {
    
    epoch <- get_epoch_label()
    bin_clean <- clean_name(bin)
    
    paste(
      "PVIR7",
      epoch,
      prefix,
      scale,
      paste0("res", blockRes),
      paste0("min", minPlt),
      group_var,
      bin_clean,
      sep = "_"
    )
  }
  
  make_caption <- function(prefix, scale, blockRes, minPlt, group_var, bin) {
    
    epoch <- get_epoch_label()
    
    group_label <- if (group_var == "all") {
      "All reference data"
    } else {
      paste0(stringr::str_to_title(group_var), ": ", bin)
    }
    
    paste(
      "PVIR7",
      epoch,
      "|", prefix,
      "|", group_label,
      "|", scale,
      "|", paste0(blockRes, "°"),
      "|", paste0("min plots = ", minPlt)
    )
  }
  
  has_sufficient_data <- function(dat, feature_col, bin, scale, blockRes, minPlt) {
    
    z <- dat[dat[[feature_col]] == bin & !is.na(dat[[feature_col]]), , drop = FALSE]
    
    if (nrow(z) == 0) {
      return(list(
        ok = FALSE,
        reason = "zero rows for bin",
        n_rows = 0,
        n_cells = 0,
        n_cells_ok = 0
      ))
    }
    
    needed <- c("POINT_X", "POINT_Y", "AGB_T_HA")
    missing_needed <- setdiff(needed, names(z))
    
    if (length(missing_needed) > 0) {
      return(list(
        ok = FALSE,
        reason = paste("missing required column(s):", paste(missing_needed, collapse = ", ")),
        n_rows = nrow(z),
        n_cells = 0,
        n_cells_ok = 0
      ))
    }
    
    z <- z[
      is.finite(z$POINT_X) &
        is.finite(z$POINT_Y) &
        is.finite(z$AGB_T_HA),
      ,
      drop = FALSE
    ]
    
    if (nrow(z) < 2) {
      return(list(
        ok = FALSE,
        reason = "fewer than 2 valid rows",
        n_rows = nrow(z),
        n_cells = 0,
        n_cells_ok = 0
      ))
    }
    
    if (scale == "nonAgg") {
      return(list(
        ok = TRUE,
        reason = "ok",
        n_rows = nrow(z),
        n_cells = NA_integer_,
        n_cells_ok = NA_integer_
      ))
    }
    
    z$Xnew <- blockRes * (0.5 + z$POINT_X %/% blockRes)
    z$Ynew <- blockRes * (0.5 + z$POINT_Y %/% blockRes)
    
    cell_counts <- z %>%
      dplyr::count(.data$Xnew, .data$Ynew, name = "n")
    
    n_cells <- nrow(cell_counts)
    n_cells_ok <- sum(cell_counts$n >= minPlt)
    
    if (n_cells_ok < 2) {
      return(list(
        ok = FALSE,
        reason = paste0("fewer than 2 aggregated cells with minPlt = ", minPlt),
        n_rows = nrow(z),
        n_cells = n_cells,
        n_cells_ok = n_cells_ok
      ))
    }
    
    list(
      ok = TRUE,
      reason = "ok",
      n_rows = nrow(z),
      n_cells = n_cells,
      n_cells_ok = n_cells_ok
    )
  }
  
  write_skip_log <- function(file_stub,
                             caption,
                             reason,
                             n_rows,
                             n_cells,
                             n_cells_ok) {
    
    write.csv(
      data.frame(
        file_stub = file_stub,
        caption = caption,
        status = "skipped",
        reason = reason,
        n_rows = n_rows,
        n_cells = n_cells,
        n_cells_ok = n_cells_ok,
        stringsAsFactors = FALSE
      ),
      file.path(resultsFolder, paste0(file_stub, "_SKIPPED.csv")),
      row.names = FALSE
    )
  }
  
  safe_save_rdata <- function(AGBdata, file_stub) {
    save(
      AGBdata,
      file = file.path(resultsFolder, paste0(file_stub, ".Rdata"))
    )
  }
  
  run_one_bin <- function(dat, feature_col, bin, file_stub) {
    
    caption <- make_caption(
      prefix = result_prefix,
      scale = scale,
      blockRes = blockRes,
      minPlt = minPlt,
      group_var = group_var,
      bin = bin
    )
    
    message("Processing: ", caption)
    
    suff <- has_sufficient_data(
      dat = dat,
      feature_col = feature_col,
      bin = bin,
      scale = scale,
      blockRes = blockRes,
      minPlt = minPlt
    )
    
    if (!suff$ok) {
      
      warning(
        "Skipping ", file_stub, ": ", suff$reason,
        " | rows = ", suff$n_rows,
        " | cells = ", suff$n_cells,
        " | cells_ok = ", suff$n_cells_ok
      )
      
      write_skip_log(
        file_stub = file_stub,
        caption = caption,
        reason = suff$reason,
        n_rows = suff$n_rows,
        n_cells = suff$n_cells,
        n_cells_ok = suff$n_cells_ok
      )
      
      return(list(status = "skipped", data = NULL))
    }
    
    AGBdata <- tryCatch(
      {
        if (scale == "nonAgg") {
          invDasymetry(
            clmn = feature_col,
            value = bin,
            wghts = TRUE,
            plotFile = dat
          )
        } else {
          invDasymetry(
            clmn = feature_col,
            value = bin,
            aggr = blockRes,
            minPlots = minPlt,
            plotFile = dat
          )
        }
      },
      error = function(e) {
        
        warning("Skipping ", file_stub, ": invDasymetry failed: ", conditionMessage(e))
        
        write_skip_log(
          file_stub = file_stub,
          caption = caption,
          reason = paste("invDasymetry failed:", conditionMessage(e)),
          n_rows = suff$n_rows,
          n_cells = suff$n_cells,
          n_cells_ok = suff$n_cells_ok
        )
        
        return(NULL)
      }
    )
    
    if (is.null(AGBdata)) {
      return(list(status = "failed", data = NULL))
    }
    
    if (!is.data.frame(AGBdata)) {
      
      warning("Skipping ", file_stub, ": invDasymetry did not return a data.frame.")
      
      write_skip_log(
        file_stub = file_stub,
        caption = caption,
        reason = "invDasymetry did not return a data.frame",
        n_rows = suff$n_rows,
        n_cells = suff$n_cells,
        n_cells_ok = suff$n_cells_ok
      )
      
      return(list(status = "failed", data = NULL))
    }
    
    if (nrow(AGBdata) < 2) {
      
      warning("Skipping ", file_stub, ": invDasymetry returned fewer than 2 rows.")
      
      write_skip_log(
        file_stub = file_stub,
        caption = caption,
        reason = "invDasymetry returned fewer than 2 rows",
        n_rows = suff$n_rows,
        n_cells = suff$n_cells,
        n_cells_ok = suff$n_cells_ok
      )
      
      return(list(status = "skipped", data = AGBdata))
    }
    
    safe_save_rdata(AGBdata, file_stub)
    
    needed_out <- c("plotAGB_10", "mapAGB")
    missing_out <- setdiff(needed_out, names(AGBdata))
    
    if (length(missing_out) > 0) {
      
      warning("Skipping Accuracy/plots for ", file_stub, ": missing output columns.")
      
      write_skip_log(
        file_stub = file_stub,
        caption = caption,
        reason = paste("missing output column(s):", paste(missing_out, collapse = ", ")),
        n_rows = suff$n_rows,
        n_cells = suff$n_cells,
        n_cells_ok = suff$n_cells_ok
      )
      
      return(list(status = "partial", data = AGBdata))
    }
    
    valid_pairs <- is.finite(AGBdata$plotAGB_10) &
      is.finite(AGBdata$mapAGB) &
      !is.na(AGBdata$plotAGB_10) &
      !is.na(AGBdata$mapAGB)
    
    n_valid <- sum(valid_pairs)
    
    if (n_valid < 2) {
      
      warning(
        "Skipping Accuracy/plots for ", file_stub,
        ": fewer than 2 valid plot-map pairs. RData was saved."
      )
      
      write_skip_log(
        file_stub = file_stub,
        caption = caption,
        reason = "fewer than 2 valid plot-map pairs after map extraction",
        n_rows = suff$n_rows,
        n_cells = suff$n_cells,
        n_cells_ok = suff$n_cells_ok
      )
      
      return(list(status = "partial", data = AGBdata))
    }
    
    AGBdata_valid <- AGBdata[valid_pairs, , drop = FALSE]
    
    tryCatch(
      {
        Accuracy(
          AGBdata_valid,
          accuracy_digits,
          resultsFolder,
          paste0(file_stub, ".csv")
        )
      },
      error = function(e) {
        warning("Accuracy failed for ", file_stub, ": ", conditionMessage(e))
      }
    )
    
    tryCatch(
      {
        OnePlot(
          AGBdata_valid$plotAGB_10,
          AGBdata_valid$mapAGB,
          resultsFolder,
          caption,
          paste0(file_stub, ".png")
        )
      },
      error = function(e) {
        warning("OnePlot failed for ", file_stub, ": ", conditionMessage(e))
      }
    )
    
    if (exists("OnePlotScatter", mode = "function")) {
      tryCatch(
        {
          OnePlotScatter(
            AGBdata_valid$plotAGB_10,
            AGBdata_valid$mapAGB,
            resultsFolder,
            caption,
            paste0(file_stub, "_scatter.png")
          )
        },
        error = function(e) {
          warning("OnePlotScatter failed for ", file_stub, ": ", conditionMessage(e))
        }
      )
    }
    
    list(status = "success", data = AGBdata)
  }
  
  # ----------------------------------------------------------------
  # Main run
  # ----------------------------------------------------------------
  
  if (is.null(df) || nrow(df) == 0) {
    warning("Input df has zero rows.")
    return(invisible(NULL))
  }
  
  df$CODE <- if ("CODE" %in% names(df)) as.character(df$CODE) else NA_character_
  
  tier_label <- get_tier_label(df)
  
  if (is.null(result_prefix)) {
    result_prefix <- tier_label
  }
  
  bin_info <- get_bins(df, group_var)
  df <- bin_info$dat
  feature_col <- bin_info$feature_col
  bins <- bin_info$bins
  
  if (!is.null(exclude_bins)) {
    bins <- setdiff(bins, exclude_bins)
  }
  
  if (length(bins) == 0) {
    warning("No valid bins found for group_var = ", group_var)
    return(invisible(NULL))
  }
  
  output_log <- data.frame(
    group_var = character(),
    bin = character(),
    scale = character(),
    status = character(),
    file_stub = character(),
    result_rdata = character(),
    result_csv = character(),
    result_png = character(),
    result_scatter_png = character(),
    stringsAsFactors = FALSE
  )
  
  last_AGBdata <- NULL
  
  for (bin in bins) {
    
    file_stub <- make_file_stub(
      prefix = result_prefix,
      scale = scale,
      blockRes = blockRes,
      minPlt = minPlt,
      group_var = group_var,
      bin = bin
    )
    
    run_result <- run_one_bin(
      dat = df,
      feature_col = feature_col,
      bin = bin,
      file_stub = file_stub
    )
    
    status <- if (is.null(run_result$status)) "unknown" else run_result$status
    
    if (!is.null(run_result$data)) {
      last_AGBdata <- run_result$data
    }
    
    output_log <- rbind(
      output_log,
      data.frame(
        group_var = group_var,
        bin = as.character(bin),
        scale = scale,
        status = status,
        file_stub = file_stub,
        result_rdata = file.path(resultsFolder, paste0(file_stub, ".Rdata")),
        result_csv = file.path(resultsFolder, paste0(file_stub, ".csv")),
        result_png = file.path(resultsFolder, paste0(file_stub, ".png")),
        result_scatter_png = file.path(resultsFolder, paste0(file_stub, "_scatter.png")),
        stringsAsFactors = FALSE
      )
    )
  }
  
  write.csv(
    output_log,
    file.path(resultsFolder, paste0(result_prefix, "_", group_var, "_", scale, "_run_log.csv")),
    row.names = FALSE
  )
  
  invisible(last_AGBdata)
}