library(terra)

projectDir <- "C:/PVIR_Reproducible"
gfcDir <- file.path(projectDir, "data", "GFC")

# Check downloaded files
list.files(gfcDir, pattern = "\\.tif$", recursive = TRUE, full.names = FALSE)

extract_gfc_tc <- function(dat,
                           gfc_dir = file.path(projectDir, "data", "GFC"),
                           target_year = 2010,
                           out_col = "tc") {
  
  pts <- terra::vect(dat, geom = c("POINT_X", "POINT_Y"), crs = "EPSG:4326")
  
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
  
  stopifnot(length(tc_files) > 0)
  
  tc_vrt <- file.path(gfc_dir, "treecover2000.vrt")
  terra::vrt(tc_files, filename = tc_vrt, overwrite = TRUE)
  tc <- terra::extract(terra::rast(tc_vrt), pts)[, 2]
  
  if (length(ly_files) > 0 && target_year > 2000) {
    ly_vrt <- file.path(gfc_dir, "lossyear.vrt")
    terra::vrt(ly_files, filename = ly_vrt, overwrite = TRUE)
    ly <- terra::extract(terra::rast(ly_vrt), pts)[, 2]
    
    # Approximate target-year tree cover:
    # treecover2000, but zero if forest loss happened up to target_year
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
  target_year = 2010,
  out_col = "tc"
)

hist(val.rm$tc)
summary(val.rm$tc)

val.rm$BIO <- ifelse(is.na(val.rm$BIO), "NA", val.rm$BIO)

write.csv(
  val.rm,
  file.path(projectDir, "outputs", "pvir_ready_demo", "agbref_addv7_with_tc.csv"),
  row.names = FALSE
)