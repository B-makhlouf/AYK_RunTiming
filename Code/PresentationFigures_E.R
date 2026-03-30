################################################################################
# PRESENTATION FIGURES — PRES-E
# NEW STANDALONE SCRIPT — do not modify any existing scripts
#
# Produces one two-panel PNG per year (2017-2022):
#   LEFT  panel: Map colored by regional group contribution — FULL YEAR
#   RIGHT panel: Map colored by regional group contribution — Q1 ONLY
#
# All stream segments in the same regional group share the same color.
# Color intensity (light pink → dark red) encodes that group's proportional
# contribution to total production (full year on the left, Q1 on the right).
#
# Saves to: <OUTPUT_DIR>/Figures/PresFigures/
#   PRES-E_two_maps_fullyear_q1_<year>.png
################################################################################

library(dplyr)
library(readr)
library(ggplot2)
library(sf)
library(scales)
library(patchwork)

################################################################################
# CONFIGURATION  — identical paths/constants to PresentationFigures.R
################################################################################

OUTPUT_DIR   <- "/Users/benjaminmakhlouf/Research_repos/02_AYK_RunTiming"
PRES_DIR     <- file.path(OUTPUT_DIR, "Figures", "PresFigures")

SPATIAL_PATH         <- "/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC_joined.shp"
BASIN_PATH           <- "/Users/benjaminmakhlouf/Desktop/Research/isoscapes_new/Kusko/Kusko_basin.shp"
MGMT_DATA_PATH       <- file.path(OUTPUT_DIR, "Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv")
CUMULATIVE_DATA_PATH <- file.path(OUTPUT_DIR, "Analysis_Results/Cumulative_Distribution/CSV/cumulative_distribution_data_Kusko_3day_intervals.csv")

MIN_STREAM_ORDER <- 3  # Kuskokwim setting from 00_setup.R

dir.create(PRES_DIR, recursive = TRUE, showWarnings = FALSE)

# ----- Cluster assignments (identical to PresentationFigures.R) -----
cluster_assignments <- tibble::tribble(
  ~mgmt_river,               ~cluster,
  "N. Fork Kusko",           "Cluster 1",
  "E. Fork Kuskokwim River", "Cluster 1",
  "E. Fork Kuskokwim",       "Cluster 1",
  "S. Fork Kusko",           "Cluster 2",
  "Upper Kusko Main",        "Cluster 2",
  "Big River",               "Cluster 2",
  "Swift",                   "Cluster 2",
  "Tatlawiksuk",             "Cluster 2",
  "Stony",                   "Cluster 3",
  "Takotna and Nixon Fork",  "Cluster 4",
  "Holitna",                 "Cluster 4",
  "George",                  "Cluster 4",
  "Hoholitna",               "Cluster 4",
  "Middle Kusko Main",       "Cluster 5",
  "Oskakawlik",              "Cluster 5",
  "Holokuk",                 "Cluster 5",
  "Kisaralik",               "Cluster 5",
  "Tuluksak",                "Cluster 5",
  "Aniak",                   "Cluster 6",
  "Kwethluk",                "Cluster 6"
)

cluster_order <- c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5","Cluster 6")

################################################################################
# DATA LOADING
################################################################################

cat("=== PRES-E: loading data ===\n")

# Management river data — provides Q1 within_quartile_prop
mgmt_data <- read_csv(MGMT_DATA_PATH, show_col_types = FALSE) %>%
  filter(!mgmt_river %in% c("Johnson", "Lower Kusko")) %>%
  left_join(cluster_assignments, by = "mgmt_river") %>%
  filter(!is.na(cluster))

available_years <- sort(unique(mgmt_data$year))
cat("  Years:", paste(available_years, collapse = ", "), "\n")

# Cumulative distribution data — provides full-year production totals
cumulative_data <- read_csv(CUMULATIVE_DATA_PATH, show_col_types = FALSE) %>%
  filter(!mgmt_river %in% c("Johnson", "Lower Kusko")) %>%
  left_join(cluster_assignments, by = "mgmt_river") %>%
  filter(!is.na(cluster))

# ------------------------------------------------------------------
# Full-year cluster production proportion per year
# (identical method to PresentationFigures.R lines 97-101)
# ------------------------------------------------------------------
full_year_cluster <- cumulative_data %>%
  group_by(year, cluster, mgmt_river) %>%
  summarise(unit_total = max(cumulative_production, na.rm = TRUE), .groups = "drop") %>%
  group_by(year, cluster) %>%
  summarise(cluster_total = sum(unit_total, na.rm = TRUE), .groups = "drop") %>%
  group_by(year) %>%
  mutate(prop_full = cluster_total / sum(cluster_total, na.rm = TRUE)) %>%
  ungroup()

# ------------------------------------------------------------------
# Q1 cluster production proportion per year
# within_quartile_prop is each river's share of that quartile's total;
# summing within a cluster gives the cluster's share of Q1 production.
# ------------------------------------------------------------------
q1_cluster <- mgmt_data %>%
  filter(toupper(trimws(quartile)) == "Q1") %>%
  group_by(year, cluster) %>%
  summarise(prop_q1 = sum(within_quartile_prop, na.rm = TRUE), .groups = "drop")

cat("  Full-year cluster proportions (sample):\n")
print(full_year_cluster %>% filter(year == available_years[1]) %>% arrange(cluster))
cat("  Q1 cluster proportions (sample):\n")
print(q1_cluster      %>% filter(year == available_years[1]) %>% arrange(cluster))

################################################################################
# SPATIAL DATA
################################################################################

cat("\n  Loading spatial data ...\n")
edges_raw <- st_read(SPATIAL_PATH, quiet = TRUE)
basin     <- st_read(BASIN_PATH,   quiet = TRUE)

edges_raw <- edges_raw[edges_raw$Str_Order >= MIN_STREAM_ORDER, ]
edges_raw <- st_transform(edges_raw, st_crs(basin))

# All edges joined to cluster info (cluster = NA for background streams)
edges_all <- edges_raw %>%
  left_join(cluster_assignments, by = "mgmt_river") %>%
  mutate(lwd_base = pmin(Str_Order, 7) * 0.18)

edges_background <- edges_all %>% filter(is.na(cluster))
edges_managed    <- edges_all %>% filter(!is.na(cluster)) %>%
  mutate(cluster = factor(cluster, levels = cluster_order))

cat("  Managed edges:", nrow(edges_managed), "   Background:", nrow(edges_background), "\n\n")

################################################################################
# HELPER: build one map panel
################################################################################

# Shared color scale limits — fixed across both panels so colours are comparable
# (set after computing the data per year inside the loop)

make_map <- function(edges_managed, edges_background, basin,
                     prop_col,          # name of the proportion column in yr_props
                     yr_props,          # data frame: cluster, <prop_col>
                     panel_title,
                     scale_limits) {

  # Join cluster-level proportion to every managed edge
  edges_yr <- edges_managed %>%
    left_join(yr_props %>% rename(fill_val = !!prop_col), by = "cluster")

  ggplot() +
    geom_sf(data = basin,
            fill = "#eef2f8", color = "grey45", linewidth = 0.55) +
    geom_sf(data = edges_background,
            color = "grey78", linewidth = 0.20, alpha = 0.50) +
    geom_sf(data = edges_yr,
            aes(color = fill_val, linewidth = lwd_base)) +
    scale_color_distiller(
      palette   = "YlOrRd",
      direction = 1,
      name     = "Proportional\ncontribution",
      limits   = scale_limits,
      labels   = function(x) paste0(round(x * 100, 0), "%"),
      na.value = "grey85"
    ) +
    scale_linewidth_identity() +
    labs(title = panel_title) +
    theme_void(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold", size = 13, hjust = 0.5,
                                      margin = margin(b = 5)),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold", size = 10),
      legend.text      = element_text(size = 9),
      legend.key.width = unit(2.0, "cm"),
      plot.background  = element_rect(fill = "transparent", color = NA),
      plot.margin      = margin(8, 8, 8, 8)
    ) +
    guides(color = guide_colorbar(title.position = "top", title.hjust = 0.5,
                                  barwidth = 8, barheight = 0.7))
}

################################################################################
# MAIN LOOP — one two-panel figure per year
################################################################################

cat("=== Generating PRES-E figures ===\n\n")

for (yr in available_years) {
  cat("  Year:", yr, "\n")

  yr_full <- full_year_cluster %>% filter(year == yr) %>% select(cluster, prop_full)
  yr_q1   <- q1_cluster        %>% filter(year == yr) %>% select(cluster, prop_q1)

  # Fix colour scale so both maps share the same range (0 to max across both)
  combined_max <- max(c(yr_full$prop_full, yr_q1$prop_q1), na.rm = TRUE)
  scale_lims   <- c(0, combined_max)

  p_full <- make_map(
    edges_managed, edges_background, basin,
    prop_col     = "prop_full",
    yr_props     = yr_full,
    panel_title  = paste0("Full Run \u2014 ", yr),
    scale_limits = scale_lims
  )

  p_q1 <- make_map(
    edges_managed, edges_background, basin,
    prop_col     = "prop_q1",
    yr_props     = yr_q1,
    panel_title  = paste0("Within Front End Closure \u2014 ", yr),
    scale_limits = scale_lims
  )

  combined <- (p_full | p_q1) + plot_layout(guides = "collect")

  out_file <- file.path(PRES_DIR,
                        sprintf("PRES-E_two_maps_fullyear_q1_%d.png", yr))
  ggsave(out_file, combined,
         width = 18, height = 8, dpi = 300, bg = "transparent")
  cat("    Saved:", basename(out_file), "\n")
}

cat("\n=== PRES-E COMPLETE ===\n")
cat("  Output folder:", PRES_DIR, "\n")
cat("  Files: PRES-E_two_maps_fullyear_q1_<year>.png\n")
cat("  Years produced:", paste(available_years, collapse = ", "), "\n\n")
