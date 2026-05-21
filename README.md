# AGBref

# AGBref data processing workflow

This repository documents the processing workflow used to generate **AGBref**, a multi-epoch and multi-resolution global reference dataset of above-ground forest biomass (AGB). The workflow extends the logic of **Plot2Map**, an R-based workflow for comparing forest plot data with biomass maps, including plot preprocessing, temporal adjustment, uncertainty estimation, map validation, forest-mask handling, and aggregation across spatial supports. :contentReference[oaicite:0]{index=0}

## Purpose

The workflow converts heterogeneous reference data sources, including National Forest Inventories (NFIs), research plots, and airborne LiDAR-derived AGB maps, into standardized gridded reference datasets for global biomass map assessment.

AGBref is generated for multiple target epochs and spatial resolutions:

- **Epochs:** 2005, 2010, 2015, 2020
- **Spatial resolutions:** 500 m, 1 km, 10 km, 25 km

The resulting files are intended for independent assessment of global and regional AGB maps, with accompanying uncertainty estimates and quality flags.

## Workflow overview

The processing workflow includes the following main steps:

1. **Input harmonization**
   - Standardize coordinates, biomass units, plot size, inventory year, biome/ecozone labels, and data-source metadata.
   - Retain only records with valid coordinates, biomass estimates, plot size, and inventory year.

2. **Forest-mask harmonization**
   - Apply a forest definition based on tree-cover thresholds.
   - Harmonize grid-cell biomass estimates to the forest area represented within each grid cell.

3. **Temporal adjustment**
   - Adjust biomass estimates to target epochs using growth-rate assumptions where plot measurement year and target epoch differ.

4. **Uncertainty estimation**
   - Estimate uncertainty components related to tree-level measurement and allometric uncertainty, temporal growth-rate uncertainty, and sampling/support mismatch.
   - Combine uncertainty components at the variance level.

5. **Spatial aggregation**
   - Aggregate plot-level or local-reference data to target grid-cell resolutions.
   - Use inverse-variance weighting where multiple plots occur within the same grid cell.

6. **Quality flag assignment**
   - Assign binary quality attributes to each grid cell based on data-quality criteria.
   - Quality flags are stored as columns in the output CSV files.

7. **Representativeness assessment**
   - Generate under-sampled raster layers showing where AGBref is less representative of biomass-related environmental conditions.
   - These rasters are provided as guidance layers for interpreting spatial coverage.

## Output files

The Zenodo archive contains:

### Master table

`AGBref_Zenodo_ready_master.csv`

This file contains all grid cells across all epochs and spatial resolutions.

### Epoch-resolution tables

Separate CSV files are also provided for each epoch and spatial resolution, for example:

- `AGBref_2005_500m.csv`
- `AGBref_2005_1km.csv`
- `AGBref_2005_10km.csv`
- `AGBref_2005_25km.csv`
- ...
- `AGBref_2020_25km.csv`

Each file retains the same column structure as the master table.

## Main columns

| Column | Description |
|---|---|
| `ID` | Unique grid-cell record ID |
| `POINT_X` | Longitude of grid-cell centroid |
| `POINT_Y` | Latitude of grid-cell centroid |
| `N` | Number of plots or reference observations within the grid cell |
| `AGB_T_HA` | Above-ground biomass estimate in Mg ha⁻¹ |
| `VAR` | Estimated variance of the grid-cell AGB estimate |
| `SD` | Estimated standard deviation of the grid-cell AGB estimate |
| `FEZ` | FAO ecological zone |
| `GEZ` | General ecological zone |
| `TC_GRID_MEAN` | Mean tree-cover percentage within the grid cell |
| `Year` | Target AGB epoch |
| `Resolution` | Spatial resolution / support of the grid-cell estimate |

## Quality flags

Quality flags are binary attributes:

- `1` = criterion satisfied
- `0` = criterion not satisfied

| Column | Description |
|---|---|
| `QUALITY_MIN_PLOTS` | Minimum plot-number criterion |
| `QUALITY_NOT_OUTDATED` | Plot acquisition within ±3 years of the target epoch |
| `QUALITY_LARGE_SIZE` | Relatively large cumulative plot area within the grid cell |
| `QUALITY_LOCALLY_REP` | Plot locations are locally representative of grid-cell tree-cover conditions |
| `QUALITY_STRICT_FILTER` | All four criteria above are satisfied |

Users can subset AGBref using any of these columns. For example, in R:

```r
AGBref <- read.csv("AGBref_Zenodo_ready_master.csv")

AGBref_strict <- subset(AGBref, QUALITY_STRICT_FILTER == 1)
AGBref_local  <- subset(AGBref, QUALITY_LOCALLY_REP == 1)


Under-sampled rasters

The archive also includes under-sampled raster layers for different filtering scenarios. These rasters are not quality flags for individual grid cells. Instead, they show where the filtered AGBref dataset is less representative of biomass-related environmental conditions.

Raster values are interpreted as:

1 = under-sampled
0 = represented / not under-sampled
NA = outside evaluated domain

These layers can be used to flag, mask, or interpret validation results in regions where AGBref has limited environmental coverage.

Visualization

A companion visualization page is available here:

https://rnvllflores.github.io/agb-ref-data-processing/

This page provides exploratory summaries of AGBref processing outputs, including spatial distributions, quality-filter effects, and representativeness layers.
