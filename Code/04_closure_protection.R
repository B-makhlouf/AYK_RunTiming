################################################################################
# CLUSTER ANALYSIS - KUSKOKWIM WATERSHED (WITH CORRECTED ANALYSIS 4)
################################################################################
# Aggregates management unit data into 6 spatial clusters:
#   Cluster 1: North Fork, East Fork
#   Cluster 2: South Fork, Upper Kusko, Big River, Swift, Tatlawiksuk
#   Cluster 3: Stony River
#   Cluster 4: Takotna/Nixon Fork, Holitna, George, Hoholitna
#   Cluster 5: Middle Kusko, Oskakawlik, Holokuk, Kisaralik, Tuloksuk
#   Cluster 6: Aniak, Kwethluk
################################################################################

library(dplyr)
library(readr)
library(ggplot2)
library(sf)
library(RColorBrewer)
library(scales)
library(cowplot)
library(patchwork)

################################################################################
# SETUP AND CONFIGURATION
################################################################################

# File paths
OUTPUT_DIR <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Cluster Analysis Results"
CSV_DIR <- file.path(OUTPUT_DIR, "CSV_Exports")
SPATIAL_PATH <- "/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC_joined.shp"
BASIN_PATH <- "/Users/benjaminmakhlouf/Desktop/Research/isoscapes_new/Kusko/Kusko_basin.shp"
MGMT_DATA_PATH <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv"
CUMULATIVE_DATA_PATH <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Cumulative_Distribution/CSV/cumulative_distribution_data_Kusko_3day_intervals.csv"

# Closure window dates (ADJUST THESE IF YOUR CLOSURE DATES ARE DIFFERENT)
CLOSURE_START_DOY <- 153  # June 1
CLOSURE_END_DOY <- 162    # June 11

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(CSV_DIR, recursive = TRUE, showWarnings = FALSE)

# Define cluster assignments
cluster_assignments <- tibble::tribble(
  ~mgmt_river,              ~cluster,
  "N. Fork Kusko",          "Cluster 1",
  "E. Fork Kuskokwim River", "Cluster 1",
  "E. Fork Kuskokwim",      "Cluster 1",  # Handle both naming conventions
  "S. Fork Kusko",          "Cluster 2",
  "Upper Kusko Main",       "Cluster 2",
  "Big River",              "Cluster 2",
  "Swift",                  "Cluster 2",
  "Tatlawiksuk",            "Cluster 2",
  "Stony",                  "Cluster 3",
  "Takotna and Nixon Fork", "Cluster 4",
  "Holitna",                "Cluster 4",
  "George",                 "Cluster 4",
  "Hoholitna",              "Cluster 4",
  "Middle Kusko Main",      "Cluster 5",
  "Oskakawlik",             "Cluster 5",
  "Holokuk",                "Cluster 5",
  "Kisaralik",              "Cluster 5",
  "Tuluksak",               "Cluster 5",
  "Aniak",                  "Cluster 6",
  "Kwethluk",               "Cluster 6"
)

cluster_order <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6")

################################################################################
# ANALYSIS 1: SPATIAL MAPS AND PRODUCTION VARIABILITY BOXPLOTS
################################################################################
# Creates:
# - Four-panel map figure with Q1-Q4 showing spatial patterns (consistent scale)
# - Boxplots showing production variability across years for each cluster
################################################################################

cat("\n=== STARTING ANALYSIS 1: SPATIAL MAPS AND BOXPLOTS ===\n")

# Load and prepare management unit data
mgmt_data <- read_csv(MGMT_DATA_PATH) %>%
  filter(mgmt_river != "Johnson", mgmt_river != "Lower Kusko") %>%
  mutate(quartile_clean = case_when(
    quartile == "Q1" ~ "Q1",
    quartile == "Q2" ~ "Q2",
    quartile == "Q3" ~ "Q3",
    quartile == "Q4" ~ "Q4",
    TRUE ~ quartile
  )) %>%
  left_join(cluster_assignments, by = "mgmt_river") %>%
  filter(!is.na(cluster))

year_range <- paste0(min(mgmt_data$year), "-", max(mgmt_data$year))

# Calculate median within-quartile proportion by cluster and quartile
# IMPORTANT: This shows the MEDIAN proportion across management units, NOT the proportion
# of the cluster's total production. Each management unit is weighted equally regardless
# of its actual production level. This shows typical timing patterns for units in each cluster.
cluster_avg_production <- mgmt_data %>%
  group_by(cluster, quartile_clean) %>%
  summarise(
    production_proportion = median(within_quartile_prop, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cluster = factor(cluster, levels = cluster_order))

# --- Load spatial data ---

edges <- st_read(SPATIAL_PATH, quiet = TRUE)
basin <- st_read(BASIN_PATH, quiet = TRUE)

# Add cluster assignments to spatial data
edges_with_clusters <- edges %>%
  left_join(cluster_assignments, by = "mgmt_river")

# --- Create four-panel map figure with consistent scale ---

# Prepare data for all quartiles
all_quartile_data <- list()
for (q in c("Q1", "Q2", "Q3", "Q4")) {
  quartile_data <- cluster_avg_production %>% filter(quartile_clean == q)
  
  edges_with_data <- edges_with_clusters %>%
    left_join(quartile_data, by = "cluster") %>%
    filter(!is.na(cluster)) %>%
    mutate(
      production_proportion = ifelse(is.na(production_proportion), 0, production_proportion),
      stream_order = ifelse(is.na(Str_Order), 3, Str_Order),
      quartile = q
    )
  
  all_quartile_data[[q]] <- edges_with_data
}

# Combine all data
combined_edges <- do.call(rbind, all_quartile_data)

# Transform basin if needed
if (st_crs(basin) != st_crs(combined_edges)) {
  basin <- st_transform(basin, st_crs(combined_edges))
}

# Create four-panel figure
p_four_panel <- ggplot() +
  geom_sf(data = basin, fill = "gray95", color = "gray70", linewidth = 0.5, alpha = 0.3) +
  geom_sf(data = combined_edges, aes(color = production_proportion, linewidth = stream_order), alpha = 0.8) +
  facet_wrap(~quartile, nrow = 1) +
  scale_color_gradientn(
    colors = brewer.pal(9, "YlOrRd"),
    name = "Production\nProportion",
    limits = c(0, 0.15),
    na.value = "grey60",
    labels = percent_format(accuracy = 1),
    guide = guide_colorbar(barwidth = 1, barheight = 15, frame.colour = "grey40", ticks.colour = "grey40")
  ) +
  scale_linewidth_continuous(range = c(0.2, 1.5), guide = "none") +
  coord_sf(datum = NA) +
  theme_void() +
  theme(
    strip.text = element_text(size = 14, face = "bold", color = "grey30"),
    legend.position = "right",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.spacing = unit(0.5, "lines")
  )

ggsave(file.path(OUTPUT_DIR, paste0("Four_Panel_Quartiles_", year_range, ".png")), 
       p_four_panel, width = 20, height = 5, dpi = 300, bg = "white")

# --- Create 8-panel figure: Maps (top row) + Boxplots (bottom row) ---

# Create production variability boxplots data
quartile_boxplot_data <- mgmt_data %>%
  mutate(
    within_quartile_pct = within_quartile_prop * 100,
    cluster = factor(cluster, levels = rev(cluster_order))
  )

# Create individual map plots without titles
map_plots <- list()
for (q in c("Q1", "Q2", "Q3", "Q4")) {
  
  quartile_edges <- combined_edges %>% filter(quartile == q)
  
  p_map <- ggplot() +
    geom_sf(data = basin, fill = "gray95", color = "gray70", linewidth = 0.5, alpha = 0.3) +
    geom_sf(data = quartile_edges, aes(color = production_proportion, linewidth = stream_order), alpha = 0.8) +
    scale_color_gradientn(
      colors = brewer.pal(9, "YlOrRd"),
      name = "Production\nProportion",
      limits = c(0, 0.15),
      na.value = "grey60",
      labels = percent_format(accuracy = 1),
      guide = guide_colorbar(barwidth = 1, barheight = 15, frame.colour = "grey40", ticks.colour = "grey40")
    ) +
    scale_linewidth_continuous(range = c(0.2, 1.5), guide = "none") +
    coord_sf(datum = NA) +
    labs(title = q) +
    theme_void() +
    theme(
      plot.title = element_text(size = 14, face = "bold", color = "grey30", hjust = 0.5),
      legend.position = "none",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  map_plots[[q]] <- p_map
}

# Create individual boxplot plots without titles
box_plots <- list()
for (q in c("Q1", "Q2", "Q3", "Q4")) {
  
  q_data <- quartile_boxplot_data %>% filter(quartile_clean == q)
  
  p_box <- ggplot(q_data, aes(x = cluster, y = within_quartile_pct)) +
    geom_boxplot(
      fill = "grey75", 
      color = "grey75",
      alpha = 0.8, 
      outlier.shape = 21,
      outlier.fill = "grey75",
      outlier.color = "grey60",
      outlier.size = 2, 
      outlier.alpha = 0.7,
      linewidth = 0.6
    ) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      limits = c(0, 20),
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_flip(clip = "off") +
    labs(
      x = NULL,
      y = "Production (%)"
    ) +
    theme_minimal(base_size = 11, base_family = "sans") +
    theme(
      axis.title.x = element_text(size = 10, color = "grey30", face = "bold", margin = margin(t = 5)),
      axis.text.x = element_text(size = 9, color = "grey40"),
      axis.text.y = element_text(size = 9, color = "grey30", face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      axis.line = element_blank(),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  box_plots[[q]] <- p_box
}

# Extract legend from one of the maps
p_legend_map <- ggplot() +
  geom_sf(data = basin, fill = "gray95", color = "gray70", linewidth = 0.5, alpha = 0.3) +
  geom_sf(data = combined_edges %>% filter(quartile == "Q1"), 
          aes(color = production_proportion, linewidth = stream_order), alpha = 0.8) +
  scale_color_gradientn(
    colors = brewer.pal(9, "YlOrRd"),
    name = "Production\nProportion",
    limits = c(0, 0.12),
    na.value = "grey60",
    labels = percent_format(accuracy = 1),
    guide = guide_colorbar(barwidth = 1, barheight = 15, frame.colour = "grey40", ticks.colour = "grey40")
  ) +
  scale_linewidth_continuous(range = c(0.2, 1.5), guide = "none") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold", color = "grey30"),
    legend.text = element_text(size = 9, color = "grey40")
  )

legend <- cowplot::get_legend(p_legend_map)

# Combine into 8-panel layout
library(patchwork)

eight_panel <- (map_plots[["Q1"]] | map_plots[["Q2"]] | map_plots[["Q3"]] | map_plots[["Q4"]]) /
  (box_plots[["Q1"]] | box_plots[["Q2"]] | box_plots[["Q3"]] | box_plots[["Q4"]]) +
  plot_layout(heights = c(1, 1))

# Combine with legend
final_plot <- cowplot::plot_grid(eight_panel, legend, rel_widths = c(1, 0.08))

ggsave(file.path(OUTPUT_DIR, paste0("Eight_Panel_Maps_Boxplots_", year_range, ".png")), 
       final_plot, width = 20, height = 10, dpi = 300, bg = "white")

# --- Create production variability boxplots ---

# Define color palette for quartiles
quartile_colors <- c("Q1" = "#2C7BB6", "Q2" = "#ABD9E9", "Q3" = "#FDAE61", "Q4" = "#D7191C")

for (q in c("Q1", "Q2", "Q3", "Q4")) {
  
  q_data <- quartile_boxplot_data %>% filter(quartile_clean == q)
  
  p <- ggplot(q_data, aes(x = cluster, y = within_quartile_pct)) +
    geom_boxplot(
      fill = "grey75", 
      color = "grey75",
      alpha = 0.8, 
      outlier.shape = 21,
      outlier.fill = "grey75",
      outlier.color = "grey60",
      outlier.size = 2.5, 
      outlier.alpha = 0.7,
      linewidth = 0.7
    ) +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      limits = c(0, 20),
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_flip(clip = "off") +
    labs(
      title = paste0(q, " Production Variability by Cluster"),
      subtitle = paste0("Distribution of within-", q, " production across years (", year_range, ")"),
      x = NULL,
      y = "Production (%)"
    ) +
    theme_minimal(base_size = 13, base_family = "sans") +
    theme(
      # Title styling
      plot.title = element_text(face = "bold", size = 18, hjust = 0, color = "grey20", margin = margin(b = 5)),
      plot.subtitle = element_text(size = 12, hjust = 0, color = "grey50", margin = margin(b = 15)),
      
      # Axis styling
      axis.title.x = element_text(size = 12, color = "grey30", face = "bold", margin = margin(t = 10)),
      axis.text.x = element_text(size = 11, color = "grey40"),
      axis.text.y = element_text(size = 12, color = "grey30", face = "bold"),
      
      # Grid styling
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      
      # Background - remove panel border
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_blank(),
      axis.line = element_blank(),
      
      # Margins
      plot.margin = margin(20, 20, 20, 20)
    )
  
  ggsave(file.path(OUTPUT_DIR, paste0("cluster_", q, "_boxplot.png")), 
         p, width = 12, height = 10, dpi = 300, bg = "white")
}

# Combined faceted boxplot
p_combined <- ggplot(quartile_boxplot_data, aes(x = cluster, y = within_quartile_pct)) +
  geom_boxplot(
    fill = "grey75",
    color = "grey75",
    alpha = 0.8, 
    outlier.shape = 21,
    outlier.fill = "grey75",
    outlier.color = "grey60",
    outlier.size = 2, 
    outlier.alpha = 0.7,
    linewidth = 0.6
  ) +
  facet_wrap(~quartile_clean, scales = "free_y", ncol = 2) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, 20),
    expand = expansion(mult = c(0, 0.02))
  ) +
  coord_flip(clip = "off") +
  labs(
    title = "Production Variability Across All Quartiles",
    subtitle = paste0("Distribution of within-quartile production by cluster (", year_range, ")"),
    x = NULL,
    y = "Production (%)"
  ) +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    # Title styling
    plot.title = element_text(face = "bold", size = 16, hjust = 0, color = "grey20", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, hjust = 0, color = "grey50", margin = margin(b = 15)),
    
    # Facet labels
    strip.text = element_text(face = "bold", size = 13, color = "grey20", margin = margin(b = 10)),
    strip.background = element_rect(fill = "grey95", color = NA),
    
    # Axis styling
    axis.title.x = element_text(size = 11, color = "grey30", face = "bold", margin = margin(t = 8)),
    axis.text.x = element_text(size = 10, color = "grey40"),
    axis.text.y = element_text(size = 10, color = "grey30", face = "bold"),
    
    # Grid styling
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    
    # Background - remove panel border
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    panel.border = element_blank(),
    axis.line = element_blank(),
    panel.spacing = unit(1.5, "lines"),
    
    # Margins
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave(file.path(OUTPUT_DIR, "cluster_all_quartiles_boxplot.png"),
       p_combined, width = 14, height = 11, dpi = 300, bg = "white")

# ==================== EXPORT ANALYSIS 1 DATA ====================
cat("--- Exporting Analysis 1 data to CSV ---\n")

# Export cluster average production by quartile
write_csv(cluster_avg_production, 
          file.path(CSV_DIR, "analysis1_cluster_avg_production_by_quartile.csv"))

# Export boxplot data (raw data with all years)
write_csv(quartile_boxplot_data %>% 
            select(year, mgmt_river, cluster, quartile_clean, within_quartile_pct, within_quartile_prop),
          file.path(CSV_DIR, "analysis1_production_variability_by_year.csv"))

# Export summary statistics for boxplots
boxplot_summary <- quartile_boxplot_data %>%
  group_by(cluster, quartile_clean) %>%
  summarise(
    mean_production_pct = mean(within_quartile_pct, na.rm = TRUE),
    median_production_pct = median(within_quartile_pct, na.rm = TRUE),
    sd_production_pct = sd(within_quartile_pct, na.rm = TRUE),
    min_production_pct = min(within_quartile_pct, na.rm = TRUE),
    max_production_pct = max(within_quartile_pct, na.rm = TRUE),
    q1_production_pct = quantile(within_quartile_pct, 0.25, na.rm = TRUE),
    q3_production_pct = quantile(within_quartile_pct, 0.75, na.rm = TRUE),
    n_years = n(),
    .groups = "drop"
  )

write_csv(boxplot_summary,
          file.path(CSV_DIR, "analysis1_boxplot_summary_statistics.csv"))

cat("=== ANALYSIS 1 COMPLETE ===\n")
cat("Created: 1 four-panel spatial map (Q1-Q4), 4 individual boxplots, 1 combined boxplot\n")
cat("Exported: 3 CSV files\n\n")

################################################################################
# ANALYSIS 2: CUMULATIVE DISTRIBUTION OVER TIME
################################################################################
# Creates:
# - Individual year plots showing cumulative production curves for all 6 clusters
# - Average plot across all years with standard error bands
# - Multi-panel plot with all years in one figure
################################################################################

cat("=== STARTING ANALYSIS 2: CUMULATIVE DISTRIBUTION ===\n")

# Load cumulative distribution data
cumulative_data <- read_csv(CUMULATIVE_DATA_PATH) %>%
  filter(mgmt_river != "Johnson", mgmt_river != "Lower Kusko") %>%
  left_join(cluster_assignments, by = "mgmt_river") %>%
  filter(!is.na(cluster))


################################################################################
# SECTION 2.1: DATA AGGREGATION
################################################################################

cat("--- Aggregating management units to cluster level ---\n")

cluster_annual_totals <- cumulative_data %>%
  group_by(year, cluster, mgmt_river) %>%
  summarise(final_production = max(cumulative_production, na.rm = TRUE), .groups = "drop") %>%
  group_by(year, cluster) %>%
  summarise(cluster_annual_total = sum(final_production, na.rm = TRUE), .groups = "drop")

cluster_aggregated <- cumulative_data %>%
  group_by(year, doy, cluster) %>%
  summarise(cluster_cumulative_production = sum(cumulative_production, na.rm = TRUE), .groups = "drop") %>%
  left_join(cluster_annual_totals, by = c("year", "cluster")) %>%
  mutate(
    cluster_cumulative_percent = (cluster_cumulative_production / cluster_annual_total) * 100,
    cluster = factor(cluster, levels = cluster_order)
  )

# Define color scheme
cluster_colors <- colorRampPalette(c("#8B0000", "#FF0000", "#FF8C00", "#87CEEB", "#1E90FF", "#000080"))(6)
names(cluster_colors) <- cluster_order

years <- sort(unique(cluster_aggregated$year))
cat("Years to plot:", paste(years, collapse = ", "), "\n\n")


################################################################################
# SECTION 2.2: INDIVIDUAL YEAR PLOTS
################################################################################

cat("--- Creating individual year plots ---\n")

for (year in years) {
  
  year_data <- cluster_aggregated %>% filter(year == !!year)
  
  p <- ggplot(year_data, aes(x = doy, y = cluster_cumulative_percent, color = cluster)) +
    geom_hline(yintercept = c(25, 50, 75), color = "gray70", linetype = "dashed", alpha = 0.7) +
    geom_line(linewidth = 2.5, alpha = 0.9) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(values = cluster_colors, breaks = cluster_order,
                       guide = guide_legend(override.aes = list(linewidth = 3))) +
    scale_y_continuous(
      labels = function(x) paste0(round(x, 1), "%"),
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100)
    ) +
    scale_x_continuous(
      breaks = seq(155, 210, by = 10),
      labels = function(x) {
        date_str <- format(as.Date(x - 1, origin = paste0(year, "-01-01")), "%b %d")
        paste0("DOY ", x, "\n", date_str)
      }
    ) +
    labs(
      title = paste("Cluster Run Timing Progress -", year),
      subtitle = "Kusko Watershed | Cumulative production by cluster",
      x = "Day of Year",
      y = "Cumulative Percent of Cluster's Total Production",
      color = "Cluster"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 13, hjust = 0.5),
      legend.position = "right",
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  ggsave(file.path(OUTPUT_DIR, paste0("cluster_cumulative_distribution_", year, ".png")),
         p, width = 14, height = 10, dpi = 300, bg = "white")
}

cat("Created", length(years), "individual year plots\n\n")


################################################################################
# SECTION 2.3: AVERAGE PLOT WITH STANDARD ERROR BANDS
################################################################################

cat("--- Creating average cumulative curve across all years ---\n")

avg_cumulative <- cluster_aggregated %>%
  group_by(cluster, doy) %>%
  summarise(
    mean_cumulative_percent = mean(cluster_cumulative_percent, na.rm = TRUE),
    se_cumulative_percent = sd(cluster_cumulative_percent, na.rm = TRUE) / sqrt(n()),
    sd_cumulative_percent = sd(cluster_cumulative_percent, na.rm = TRUE),
    n_years = n(),
    .groups = "drop"
  ) %>%
  mutate(cluster = factor(cluster, levels = cluster_order))

p_avg <- ggplot(avg_cumulative, aes(x = doy, y = mean_cumulative_percent, color = cluster, fill = cluster)) +
  geom_ribbon(aes(ymin = mean_cumulative_percent - se_cumulative_percent,
                  ymax = mean_cumulative_percent + se_cumulative_percent),
              alpha = 0.2, color = NA) +
  geom_hline(yintercept = c(25, 50, 75), color = "gray70", linetype = "dashed", alpha = 0.7) +
  geom_line(linewidth = 2.5, alpha = 0.9) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = cluster_colors, breaks = cluster_order,
                     guide = guide_legend(override.aes = list(linewidth = 3))) +
  scale_fill_manual(values = cluster_colors, breaks = cluster_order, guide = "none") +
  scale_y_continuous(
    labels = function(x) paste0(round(x, 1), "%"),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  scale_x_continuous(
    breaks = seq(155, 210, by = 10),
    labels = function(x) {
      date_str <- format(as.Date(x - 1, origin = "2020-01-01"), "%b %d")
      paste0("DOY ", x, "\n", date_str)
    }
  ) +
  labs(
    title = paste("Average Cluster Run Timing -", min(years), "to", max(years)),
    subtitle = "Kusko Watershed | Mean cumulative production (± SE)",
    x = "Day of Year",
    y = "Cumulative Percent of Cluster's Total Production",
    color = "Cluster",
    caption = "Shaded bands show ±1 standard error"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    legend.position = "right",
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(OUTPUT_DIR, "cluster_cumulative_distribution_AVERAGE.png"),
       p_avg, width = 14, height = 10, dpi = 300, bg = "white")

cat("Created average plot with standard error bands\n\n")


################################################################################
# SECTION 2.4: MULTI-PANEL PLOT WITH ALL YEARS
################################################################################

cat("--- Creating multi-panel plot with all years ---\n")

# Calculate plot height based on number of years
n_years <- length(years)
plot_height <- max(15, n_years * 2.5)

p_multi <- ggplot(cluster_aggregated, aes(x = doy, y = cluster_cumulative_percent, color = cluster)) +
  geom_hline(yintercept = c(25, 50, 75), color = "gray70", linetype = "dashed", alpha = 0.5) +
  geom_line(linewidth = 1.2, alpha = 0.9) +
  geom_point(size = 0.8, alpha = 0.7) +
  scale_color_manual(
    values = cluster_colors, 
    breaks = cluster_order,
    guide = guide_legend(override.aes = list(linewidth = 2))
  ) +
  scale_y_continuous(
    labels = function(x) paste0(round(x, 0), "%"),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  scale_x_continuous(
    breaks = seq(155, 210, by = 10),
    labels = function(x) {
      date_str <- format(as.Date(x - 1, origin = "2020-01-01"), "%b %d")
      date_str
    }
  ) +
  facet_wrap(~ year, ncol = 1, scales = "free_x", strip.position = "top") +
  labs(
    title = "Cumulative distribution of returning Chinook salmon over time",
    x = "Date",
    y = "Cumulative Percent of Cluster's Total Production",
    color = "Cluster"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "white", color = NA),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    panel.spacing.y = unit(0.5, "lines"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray80", linewidth = 0.3),
    panel.border = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 11, color = "black")
  )

ggsave(
  file.path(OUTPUT_DIR, "cluster_cumulative_distribution_MULTIPANEL.png"),
  p_multi, 
  width = 10, 
  height = plot_height, 
  dpi = 300, 
  bg = "white",
  limitsize = FALSE
)

cat("Created multi-panel plot with all years\n")
cat("  Plot dimensions: 10 x", plot_height, "inches\n")
cat("  Number of panels:", n_years, "\n\n")

# ==================== EXPORT ANALYSIS 2 DATA ====================
cat("--- Exporting Analysis 2 data to CSV ---\n")

# Export cluster-aggregated cumulative data (all years)
write_csv(cluster_aggregated %>%
            select(year, doy, cluster, cluster_cumulative_production, 
                   cluster_annual_total, cluster_cumulative_percent),
          file.path(CSV_DIR, "analysis2_cluster_cumulative_by_year.csv"))

# Export average cumulative data with SE
write_csv(avg_cumulative,
          file.path(CSV_DIR, "analysis2_average_cumulative_with_se.csv"))

################################################################################
# SECTION 2.5: ANALYSIS SUMMARY
################################################################################

cat("=== ANALYSIS 2 COMPLETE ===\n")
cat("Total figures created:\n")
cat("  - Individual year plots:", length(years), "\n")
cat("  - Average plot: 1\n")
cat("  - Multi-panel plot: 1\n")
cat("  TOTAL:", length(years) + 2, "figures\n")
cat("Exported: 2 CSV files\n\n")


################################################################################
# ANALYSIS 3: DAYS TO 50% RETURN BY CLUSTER
################################################################################
# Creates:
# - Horizontal range plot showing min/max/mean days to reach 50% of return
# - Each cluster shows variability across years
################################################################################

cat("=== STARTING ANALYSIS 3: DAYS TO 50% RETURN ===\n")

# --- Calculate days to 50% for each management unit and year ---

days_to_50_mgmt <- cumulative_data %>%
  group_by(year, mgmt_river) %>%
  arrange(year, mgmt_river, doy) %>%
  summarise(
    run_start_doy = min(doy, na.rm = TRUE),
    doy_above_50 = {
      above_50_idx <- which(cumulative_percent > 50)
      if (length(above_50_idx) > 0) doy[min(above_50_idx)] else NA_real_
    },
    .groups = "drop"
  ) %>%
  filter(!is.na(doy_above_50)) %>%
  mutate(days_to_50_percent = doy_above_50 - run_start_doy)

# --- Aggregate to cluster level ---

days_to_50_cluster <- days_to_50_mgmt %>%
  left_join(cluster_assignments, by = "mgmt_river") %>%
  filter(!is.na(cluster)) %>%
  group_by(year, cluster) %>%
  summarise(
    cluster_days_to_50 = mean(days_to_50_percent, na.rm = TRUE),
    .groups = "drop"
  )

# --- Calculate summary statistics (min, max, mean across years) ---

days_to_50_summary <- days_to_50_cluster %>%
  group_by(cluster) %>%
  summarise(
    mean_days = mean(cluster_days_to_50, na.rm = TRUE),
    median_days = median(cluster_days_to_50, na.rm = TRUE),
    sd_days = sd(cluster_days_to_50, na.rm = TRUE),
    min_days = min(cluster_days_to_50, na.rm = TRUE),
    max_days = max(cluster_days_to_50, na.rm = TRUE),
    n_years = n(),
    .groups = "drop"
  ) %>%
  mutate(cluster = factor(cluster, levels = cluster_order))

# --- Create horizontal range plot ---

p_days_to_50 <- ggplot(days_to_50_summary, aes(y = cluster)) +
  # Range lines (min to max)
  geom_segment(aes(x = min_days, xend = max_days, color = cluster), 
               linewidth = 6, alpha = 0.7) +
  
  # Mean markers (black vertical lines)
  geom_point(aes(x = mean_days), color = "black", size = 4, shape = "|", stroke = 2) +
  
  # Colors
  scale_color_manual(values = cluster_colors, breaks = cluster_order, guide = "none") +
  
  scale_x_continuous(
    name = "Days to 50% return",
    breaks = function(x) pretty(x, n = 6),
    expand = expansion(mult = c(0.02, 0.02)),
    labels = function(x) paste0(x, " days")
  ) +
  
  scale_y_discrete(name = NULL, limits = rev(cluster_order)) +
  
  labs(title = "Days to 50% of return by stock") +
  
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 20, hjust = 0, color = "gray15", margin = margin(b = 5)),
    axis.text.y = element_text(size = 11, color = "gray30", hjust = 1, margin = margin(r = 10)),
    axis.text.x = element_text(size = 10, color = "gray30", margin = margin(t = 8)),
    axis.title.x = element_text(size = 12, color = "gray20", face = "bold", margin = margin(t = 15)),
    panel.grid.major.x = element_line(color = "gray92", linewidth = 0.4, linetype = "solid"),
    panel.grid.minor.x = element_line(color = "gray96", linewidth = 0.2),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(20, 20, 20, 20)
  )

ggsave(file.path(OUTPUT_DIR, "cluster_days_to_50_percent.png"),
       p_days_to_50, width = 10, height = 6, dpi = 300, bg = "white")

# ==================== EXPORT ANALYSIS 3 DATA ====================
cat("--- Exporting Analysis 3 data to CSV ---\n")

# Export raw days to 50% by year and cluster
write_csv(days_to_50_cluster,
          file.path(CSV_DIR, "analysis3_days_to_50_by_year.csv"))

# Export summary statistics
write_csv(days_to_50_summary,
          file.path(CSV_DIR, "analysis3_days_to_50_summary.csv"))

cat("=== ANALYSIS 3 COMPLETE ===\n")
cat("Created: 1 range plot showing days to 50% return\n")
cat("Exported: 2 CSV files\n\n")

################################################################################
# ANALYSIS 4: FRONT-END CLOSURE PROTECTION (CORRECTED)
################################################################################
# Uses ACTUAL CLOSURE WINDOW DATES (not Q1) to match management unit analysis
# Default: DOY 153-162 (approximately June 1-11)
################################################################################

cat("=== STARTING ANALYSIS 4: FRONT-END CLOSURE PROTECTION (CORRECTED) ===\n")
cat("Using closure window: DOY", CLOSURE_START_DOY, "to", CLOSURE_END_DOY, "\n")
cat("(Approximately June 1-11 - adjust CLOSURE_START_DOY and CLOSURE_END_DOY if different)\n\n")

################################################################################
# STEP 4.1: CALCULATE CLOSURE WINDOW PRODUCTION FOR EACH UNIT
################################################################################

cat("Calculating closure window production for each management unit...\n")

closure_protection_by_unit <- cumulative_data %>%
  group_by(year, mgmt_river, cluster) %>%
  arrange(doy) %>%
  summarise(
    # Total annual production
    total_annual_production = max(cumulative_production, na.rm = TRUE),
    
    # Production during closure window
    closure_window_production = {
      # Find production at the last day of closure window (or closest day)
      closure_days <- doy[doy >= CLOSURE_START_DOY & doy <= CLOSURE_END_DOY]
      if (length(closure_days) > 0) {
        max_closure_doy <- max(closure_days)
        cumulative_production[doy == max_closure_doy][1]
      } else {
        0  # No fish during closure window
      }
    },
    
    # Check if we have data for closure window
    has_closure_data = any(doy >= CLOSURE_START_DOY & doy <= CLOSURE_END_DOY),
    
    .groups = "drop"
  ) %>%
  mutate(
    # Handle cases where closure window production is NA
    closure_window_production = ifelse(is.na(closure_window_production), 0, closure_window_production),
    
    # Calculate protection percentage for this unit
    unit_protection_pct = (closure_window_production / total_annual_production) * 100
  )

cat("Unit-level protection calculated for", nrow(closure_protection_by_unit), "unit-year combinations\n")

# Check for any issues
issues <- closure_protection_by_unit %>% 
  filter(!has_closure_data | is.na(unit_protection_pct) | unit_protection_pct > 100)

if (nrow(issues) > 0) {
  cat("⚠ Warning: Found", nrow(issues), "unit-years with data issues\n")
}

################################################################################
# STEP 4.2: AGGREGATE TO CLUSTER LEVEL
################################################################################

cat("\nAggregating to cluster level...\n")

cluster_closure_protection <- closure_protection_by_unit %>%
  group_by(year, cluster) %>%
  summarise(
    # Sum production across all units in cluster
    cluster_total_production = sum(total_annual_production, na.rm = TRUE),
    cluster_closure_production = sum(closure_window_production, na.rm = TRUE),
    n_units = n(),
    .groups = "drop"
  ) %>%
  mutate(
    # Calculate protection percentage for the cluster
    closure_protection_pct = (cluster_closure_production / cluster_total_production) * 100,
    cluster = factor(cluster, levels = rev(cluster_order))
  )

# Summary statistics
cat("\n--- Closure Protection by Cluster ---\n")
closure_summary <- cluster_closure_protection %>%
  group_by(cluster) %>%
  summarise(
    mean_pct = round(mean(closure_protection_pct), 1),
    min_pct = round(min(closure_protection_pct), 1),
    max_pct = round(max(closure_protection_pct), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_pct))

print(closure_summary)

################################################################################
# STEP 4.3: CREATE VISUALIZATION
################################################################################

cat("\n--- Creating visualization ---\n")

p_closure <- ggplot(cluster_closure_protection, 
                    aes(y = cluster, x = closure_protection_pct, 
                        color = cluster, fill = cluster)) +
  geom_violin(alpha = 0.15, scale = "width", width = 0.6, color = NA) +
  stat_summary(fun.min = min, fun.max = max, geom = "linerange", 
               linewidth = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "point", size = 4, 
               color = "white", stroke = 2, shape = 21) +
  scale_color_manual(values = rev(cluster_colors), guide = "none") +
  scale_fill_manual(values = rev(cluster_colors), guide = "none") +
  scale_x_continuous(
    name = "Protection Effectiveness (%)",
    labels = function(x) paste0(round(x, 1), "%"),
    expand = expansion(mult = c(0.02, 0.05)),
    limits = c(0, NA)
  ) +
  scale_y_discrete(name = NULL) +
  labs(
    title = "Front-End Closure Protection by Cluster",
    subtitle = paste0("% of cluster production during DOY ", CLOSURE_START_DOY, "-", 
                      CLOSURE_END_DOY, " (", year_range, ")"),
    caption = "White dots = mean | Colored areas show distribution | Range lines show min-max\nCalculated using actual closure window dates (not Q1 proxy)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 0.5, color = "gray50"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 12),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(OUTPUT_DIR, "cluster_front_end_closure_protection.png"),
       p_closure, width = 12, height = 10, dpi = 300, bg = "white")

cat("✓ Figure saved\n")

# ==================== EXPORT ANALYSIS 4 DATA ====================
cat("--- Exporting Analysis 4 data to CSV ---\n")

# Export cluster-level data
write_csv(cluster_closure_protection %>% 
            select(year, cluster, closure_protection_pct, n_units,
                   cluster_total_production, cluster_closure_production),
          file.path(CSV_DIR, "analysis4_cluster_closure_protection_by_year.csv"))

write_csv(closure_summary,
          file.path(CSV_DIR, "analysis4_cluster_closure_protection_summary.csv"))

# Export unit-level data for transparency
write_csv(closure_protection_by_unit %>% 
            select(year, mgmt_river, cluster, total_annual_production, 
                   closure_window_production, unit_protection_pct),
          file.path(CSV_DIR, "analysis4_unit_closure_protection_by_year.csv"))

cat("=== ANALYSIS 4 COMPLETE ===\n")
cat("Created: 1 violin plot showing front-end closure protection\n")
cat("Exported: 3 CSV files\n\n")

################################################################################
# SUMMARY
################################################################################

cat("=== ALL ANALYSES COMPLETE ===\n\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("CSV output directory:", CSV_DIR, "\n\n")
cat("Files created:\n")
cat("  Analysis 1 (Spatial & Boxplots):\n")
cat("    - 1 four-panel map (Four_Panel_Quartiles_*.png)\n")
cat("    - 1 eight-panel map+boxplot (Eight_Panel_Maps_Boxplots_*.png)\n")
cat("    - 4 individual boxplots (cluster_Q*_boxplot.png)\n")
cat("    - 1 combined boxplot (cluster_all_quartiles_boxplot.png)\n")
cat("    - 3 CSV files\n")
cat("  Analysis 2 (Cumulative Distribution):\n")
cat("    -", length(years), "individual year plots (cluster_cumulative_distribution_*.png)\n")
cat("    - 1 average plot (cluster_cumulative_distribution_AVERAGE.png)\n")
cat("    - 1 multi-panel plot (cluster_cumulative_distribution_MULTIPANEL.png)\n")
cat("    - 2 CSV files\n")
cat("  Analysis 3 (Days to 50%):\n")
cat("    - 1 range plot (cluster_days_to_50_percent.png)\n")
cat("    - 2 CSV files\n")
cat("  Analysis 4 (Front-End Closure Protection - CORRECTED):\n")
cat("    - 1 violin plot (cluster_front_end_closure_protection.png)\n")
cat("    - 3 CSV files\n")
cat("\nTotal:", 7 + length(years) + 3, "figures created\n")
cat("Total: 11 CSV files exported\n")
cat("\n✓ SCRIPT COMPLETE ✓\n")