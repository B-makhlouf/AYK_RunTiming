################################################################################
# config.R  --  SINGLE SOURCE OF TRUTH FOR ALL PATHS AND SHARED PARAMETERS
################################################################################
# Every analysis script in R/ reads its paths and parameters from the objects
# defined here. This is the ONLY file you should need to edit to run the
# pipeline on your own machine.
#
# HOW PATHS ARE RESOLVED
#   PROJECT_ROOT is detected automatically:
#     1. If the 'here' package is installed, here::here() is used (anchored by
#        the AYK_RunTiming_Publication.Rproj file in this folder).
#     2. Otherwise the current working directory is used, so make sure R's
#        working directory is set to this folder before sourcing.
#
# WHAT YOU MUST SUPPLY
#   The strontium-isoscape stream network and basin shapefiles are NOT
#   redistributed here. Place the two Kuskokwim shapefiles (and their sidecar
#   .dbf/.shx/.prj files) under data/spatial/, preserving the SAME folder
#   structure they had in the original analysis:
#       data/spatial/Spatial Data/KuskoUSGS_HUC_joined.shp
#       data/spatial/isoscapes_new/Kusko/Kusko_basin.shp
#   (These originate from French et al. 2020; see README.)
################################################################################

# ---- PROJECT ROOT -----------------------------------------------------------
if (requireNamespace("here", quietly = TRUE)) {
  PROJECT_ROOT <- here::here()
} else {
  PROJECT_ROOT <- getwd()
}
message("config.R: PROJECT_ROOT = ", PROJECT_ROOT)

# ---- INPUT DATA (shipped with this repository) ------------------------------
NATAL_DATA_DIR <- file.path(PROJECT_ROOT, "data", "natal_origins")

# ---- INPUT DATA YOU MUST SUPPLY (spatial layers; see header) ----------------
# Sub-paths under data/spatial/ replicate the ORIGINAL analysis folder structure
# ("Spatial Data/..." for the stream edges, "isoscapes_new/Kusko/..." for the
# basin), so the source folders can be copied in wholesale without renaming.
SPATIAL_DIR <- file.path(PROJECT_ROOT, "data", "spatial")
KUSKO_EDGES <- file.path(SPATIAL_DIR, "Spatial Data", "KuskoUSGS_HUC_joined.shp")
KUSKO_BASIN <- file.path(SPATIAL_DIR, "isoscapes_new", "Kusko", "Kusko_basin.shp")

# ---- OUTPUT LOCATIONS -------------------------------------------------------
OUTPUT_ROOT          <- file.path(PROJECT_ROOT, "outputs")
MAPS_DIR             <- file.path(OUTPUT_ROOT, "Basin_Maps")
FIGURES_DIR          <- file.path(OUTPUT_ROOT, "Figures")
MGMT_RIVER_DIR       <- file.path(OUTPUT_ROOT, "Management_River_Analysis")
CUMULATIVE_DIR       <- file.path(OUTPUT_ROOT, "Cumulative_Distribution")
CPUE_DIR             <- file.path(OUTPUT_ROOT, "RemainingCPUE")
CLUSTER_GROUPS_DIR   <- file.path(OUTPUT_ROOT, "Management_Distinguishability")   # Fig 2B
CLUSTER_ANALYSIS_DIR <- file.path(OUTPUT_ROOT, "Cluster Analysis Results")        # Figs 3-6

# ---- INTERMEDIATE FILES (written upstream, read downstream) ------------------
# These paths MUST be consistent between the script that writes them and the
# cluster-analysis script that reads them. Do not change one without the other.
MGMT_TIDY_CSV  <- file.path(MGMT_RIVER_DIR, "management_river_analysis_tidy.csv")
CUMULATIVE_CSV <- file.path(CUMULATIVE_DIR, "CSV",
                            "cumulative_distribution_data_Kusko_3day_intervals.csv")
DAILY_CPUE_CSV <- file.path(CPUE_DIR, "daily_cpue_with_remaining_proportion.csv")

# ---- SHARED ANALYSIS PARAMETERS ---------------------------------------------
# (mirrors the values used to produce the published figures; unchanged)
CFG <- list(
  years                   = c(2017, 2018, 2019, 2020, 2021, 2022),
  watersheds              = c("Kusko"),
  cumulative_interval_days = 3,
  closure_start_doy       = 153,   # June 1  (used by 06_cluster_analysis.R)
  closure_end_doy         = 162,   # June 11 (used by 06_cluster_analysis.R)
  kmeans_seed             = 42,    # set.seed() value used in clustering
  kmeans_nstart           = 25
)

invisible(TRUE)
