################################################################################
# run_all.R  --  REPRODUCE THE CanJFS ANALYSIS END TO END
################################################################################
# Runs the complete pipeline in dependency order. Each stage writes its outputs
# into outputs/ ; later stages read the intermediates written by earlier stages.
#
# BEFORE RUNNING
#   1. Install the required packages (see README "Requirements").
#   2. Put the two Kuskokwim shapefiles in data/spatial/Spatial Data/.
#   3. Open AYK_RunTiming_Publication.Rproj (or setwd() to this folder) so the
#      relative paths resolve.
#
# Then:  source("run_all.R")
################################################################################

t_start <- Sys.time()

# ---- Configuration (paths + parameters: the single source of truth) ----------
source("config.R")
dir.create(OUTPUT_ROOT, recursive = TRUE, showWarnings = FALSE)

# ---- Stage 0: core functions -------------------------------------------------
source("R/00_setup.R")
source("R/01_plot_functions.R")

# ---- Stage 1: natal origin -> tributary / management-unit production ---------
# Writes outputs/Management_River_Analysis/management_river_analysis_tidy.csv
source("R/02_tributary_management.R")
run_tributary_analysis(years = CONFIG$years,
                       watersheds = CONFIG$watersheds,
                       export_csv = TRUE)

# ---- Stage 2: cumulative distribution over the season ------------------------
# Writes outputs/Cumulative_Distribution/CSV/cumulative_distribution_data_Kusko_3day_intervals.csv
source("R/03_cumulative_distribution.R")
run_cumulative_analysis(years = CONFIG$years,
                        watersheds = CONFIG$watersheds,
                        interval_days = CONFIG$cumulative_interval_days,
                        export_csv = TRUE)

# ---- Stage 3: daily CPUE + remaining proportion ------------------------------
source("R/04_cpue_remaining.R")

# ---- Stage 4: k-means regional groups + clustering diagnostics ---------------
# Writes the cluster-selection (elbow/silhouette) figure and the per-cluster
# isotope boxplots, plus the cluster-assignment CSVs.
source("R/05_clustering_regional_groups.R")

# ---- Stage 5: cluster analyses 1-6 (Figures 3, 4, 5, 6) ----------------------
source("R/06_cluster_analysis.R")

# ---- Stage 6: collect final figures into Figures/ in paper order -------------
# Figures/Figure3_*.png .. Figure6_*.png  and  Figures/Clustering/*.png
source("R/07_export_manuscript_figures.R")

# ---- Done --------------------------------------------------------------------
elapsed <- round(difftime(Sys.time(), t_start, units = "mins"), 2)
cat("\n================================================================\n")
cat("PIPELINE COMPLETE -", elapsed, "minutes\n")
cat("Working outputs:", OUTPUT_ROOT, "\n")
cat("Final figures:  ", FIGURES_EXPORT_DIR, "\n")
cat("================================================================\n")
