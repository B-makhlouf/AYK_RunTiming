################################################################################
# 07_EXPORT_MANUSCRIPT_FIGURES.R
################################################################################
# Collects the final figures into a clean, paper-ordered Figures/ folder:
#
#   Figures/                         <- manuscript figures, named in paper order
#     Figure3_QuartileCont.png             (R/06 Analysis 1, combined 3a+3b)
#     Figure3a_contribution_by_group_faceted.png          (R/06 Analysis 1)
#     Figure3b_group_contribution_stacked.png             (R/06 Analysis 1)
#     Figure4_cumulative_distribution_by_group.png        (R/06 Analysis 2)
#     Figure5_closure_protection_by_offset.png            (R/06 Analysis 5)
#     Figure6_protection_vs_foregone_harvest.png          (R/06 Analysis 6)
#     Figure_6_ALT.png                                     (R/06 Analysis 6, axes flipped + 1:1 line)
#   Figures/Clustering/              <- clustering-results figures
#     cluster_selection_elbow_silhouette.png   (4-panel incl. silhouette score)
#     isotope_value_boxplots_by_cluster.png    (boxplots of values per cluster)
#
# Pure file copy: each figure is byte-identical to what the analysis scripts
# render (same data, same rendering) - only the filenames change.
# Manuscript Figures 1, 2A, 2B are GIS-made basemaps, not generated here.
################################################################################

if (!exists("PROJECT_ROOT")) source("config.R")

dir.create(FIGURES_EXPORT_DIR,     recursive = TRUE, showWarnings = FALSE)
dir.create(CLUSTERING_FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

year_range <- paste0(min(CFG$years), "-", max(CFG$years))

manuscript_figs <- list(
  list(dest = "Figure3_QuartileCont.png",
       src  = file.path(CLUSTER_ANALYSIS_DIR, paste0("Figure3_QuartileCont_", year_range, ".png"))),
  list(dest = "Figure3a_contribution_by_group_faceted.png",
       src  = file.path(CLUSTER_ANALYSIS_DIR, paste0("Figure3a_contribution_by_group_faceted_", year_range, ".png"))),
  list(dest = "Figure3b_group_contribution_stacked.png",
       src  = file.path(CLUSTER_ANALYSIS_DIR, paste0("Figure3b_group_contribution_stacked_", year_range, ".png"))),
  list(dest = "Figure4_cumulative_distribution_by_group.png",
       src  = file.path(CLUSTER_ANALYSIS_DIR, "cluster_cumulative_distribution_MULTIPANEL.png")),
  list(dest = "Figure5_closure_protection_by_offset.png",
       src  = file.path(CLUSTER_ANALYSIS_DIR, "cluster_front_end_closure_protection_comparison.png")),
  list(dest = "Figure6_protection_vs_foregone_harvest.png",
       src  = file.path(CLUSTER_ANALYSIS_DIR, "cpue_vs_protection_combined.png")),
  list(dest = "Figure_6_ALT.png",
       src  = file.path(CLUSTER_ANALYSIS_DIR, "cpue_vs_protection_combined_ALT.png"))
)

clustering_figs <- list(
  list(dest = "cluster_selection_elbow_silhouette.png",
       src  = file.path(CLUSTER_GROUPS_DIR, "elbow_method_analysis_combined_upper_kusko.png")),
  list(dest = "isotope_value_boxplots_by_cluster.png",
       src  = file.path(CLUSTER_GROUPS_DIR, "cluster_isotope_boxplots_k6_combined_upper_kusko.png"))
)

copy_set <- function(items, out_dir, label) {
  cat("\n=== ", label, " -> ", out_dir, " ===\n", sep = "")
  missing <- character(0)
  for (m in items) {
    dest <- file.path(out_dir, m$dest)
    if (file.exists(m$src)) {
      file.copy(m$src, dest, overwrite = TRUE)
      cat(sprintf("  %-52s <- %s\n", m$dest, basename(m$src)))
    } else {
      missing <- c(missing, m$dest)
      cat(sprintf("  [MISSING] %s\n            (source not found: %s)\n", m$dest, m$src))
    }
  }
  missing
}

miss <- c(
  copy_set(manuscript_figs, FIGURES_EXPORT_DIR,     "Manuscript figures"),
  copy_set(clustering_figs, CLUSTERING_FIGURES_DIR, "Clustering figures")
)

if (length(miss) > 0) {
  cat("\n  WARNING:", length(miss), "figure(s) missing - run the full pipeline first (run_all.R).\n")
}
cat("\nManuscript Figures 1, 2A, 2B are GIS-made maps and are not generated here (see README).\n")
cat("=== Figure export complete ===\n")
