################################################################################
# run_all.R  --  REPRODUCE THE CanJFS ANALYSIS END TO END
################################################################################
# Runs the complete pipeline in dependency order. Each stage writes its outputs
# (CSV + figures) into outputs/ ; later stages read the intermediates written by
# earlier stages.
#
# BEFORE RUNNING
#   1. Install the required packages (see README "Requirements").
#   2. Put the two Kuskokwim shapefiles in data/spatial/ (see config.R / README).
#   3. Open AYK_RunTiming_Publication.Rproj in RStudio, OR setwd() to this
#      folder, so that all relative paths resolve correctly.
#
# Then:  source("run_all.R")
################################################################################

t_start <- Sys.time()

# ---- Configuration (paths + parameters: the single source of truth) ----------
source("config.R")
dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

# ---- Stage 0: core functions -------------------------------------------------
# 00_setup.R  : parameters + the Bayesian natal-origin assignment engine
# 01_plot...  : shared plotting functions (tributary maps, cumulative overview)
source("R/00_setup.R")
source("R/01_plot_functions.R")

# ---- Stage 1: natal origin -> tributary / management-unit production ---------
# Writes: outputs/Management_River_Analysis/management_river_analysis_tidy.csv
#         outputs/Basin_Maps/DOY_Quartile/... (tributary + management maps; Fig 2A material)
# Requires: data/spatial shapefiles + data/natal_origins CSVs
source("R/02_tributary_management.R")
run_tributary_analysis(years = CONFIG$years,
                       watersheds = CONFIG$watersheds,
                       export_csv = TRUE)

# ---- Stage 2: cumulative distribution over the season ------------------------
# Writes: outputs/Cumulative_Distribution/CSV/cumulative_distribution_data_Kusko_3day_intervals.csv
#         outputs/Cumulative_Distribution/Plots/cumulative_overview_Kusko_PAPER.png
source("R/03_cumulative_distribution.R")
run_cumulative_analysis(years = CONFIG$years,
                        watersheds = CONFIG$watersheds,
                        interval_days = CONFIG$cumulative_interval_days,
                        export_csv = TRUE)

# ---- Stage 3: daily CPUE + remaining proportion ------------------------------
# Writes: outputs/RemainingCPUE/daily_cpue_with_remaining_proportion.csv
# (sourcing this script runs it)
source("R/04_cpue_remaining.R")

# ---- Stage 4: k-means regional groups (Figure 2B) ----------------------------
# Writes: outputs/Management_Distinguishability/ (elbow/silhouette metrics,
#         k=3..6 cluster assignments + basin maps, k=6 violin plot)
# (sourcing this script runs run_isotope_clustering())
source("R/05_clustering_regional_groups.R")

# ---- Stage 5: cluster analyses 1-6 (Figures 3, 4, 5, 6) ----------------------
# Reads the three intermediates produced above and the spatial layers.
# Writes: outputs/Cluster Analysis Results/ (figures + CSV_Exports)
# (sourcing this script runs all six analyses)
source("R/06_cluster_analysis.R")

# ---- Done --------------------------------------------------------------------
elapsed <- round(difftime(Sys.time(), t_start, units = "mins"), 2)
cat("\n================================================================\n")
cat("PIPELINE COMPLETE -", elapsed, "minutes\n")
cat("All outputs written under:", OUTPUT_ROOT, "\n")
cat("================================================================\n")
