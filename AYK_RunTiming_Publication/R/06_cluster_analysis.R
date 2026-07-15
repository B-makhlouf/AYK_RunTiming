################################################################################
# CLUSTER ANALYSIS - KUSKOKWIM WATERSHED (WITH CORRECTED ANALYSIS 1, 2 & 4)
# MODIFIED: PUBLICATION-READY BOXPLOT LABEL SIZES
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
library(tidyr)
library(grid)  # ADDED for textGrob function

################################################################################
# SETUP AND CONFIGURATION
################################################################################

# File paths come from config.R (single source of truth)
if (!exists("PROJECT_ROOT")) source("config.R")
OUTPUT_DIR <- CLUSTER_ANALYSIS_DIR
CSV_DIR <- file.path(OUTPUT_DIR, "CSV_Exports")
SPATIAL_PATH <- KUSKO_EDGES
BASIN_PATH <- KUSKO_BASIN
MGMT_DATA_PATH <- MGMT_TIDY_CSV              # written by R/02_tributary_management.R
CUMULATIVE_DATA_PATH <- CUMULATIVE_CSV       # written by R/03_cumulative_distribution.R

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

# *** ADD THIS HERE - DEFINE CLUSTER COLORS EARLY ***
# Define cluster colors (used across all analyses)
cluster_colors <- colorRampPalette(c("#8B0000", "#FF0000", "#FF8C00", "#87CEEB", "#1E90FF", "#000080"))(6)
names(cluster_colors) <- cluster_order

################################################################################
# ANALYSIS 1: SPATIAL MAPS AND PRODUCTION VARIABILITY BOXPLOTS (CORRECTED)
################################################################################
# GOAL: For each quartile, calculate what percentage of that quartile's fish
#       come from each cluster (must sum to 100% within each quartile)
#       Do this for EACH YEAR separately, then average across years
#
# MAP FIX: Maps now colored by CLUSTER MEDIAN values (matching boxplots)
#
# MODIFICATIONS:
# - Only exports 8-panel combined figure (maps + boxplots)
# - Increased boxplot contrast (darker outlines, thicker median)
# - Single x-axis label across all panels
# - PUBLICATION-READY: Much larger text throughout (size 22-28)
# - Figure size increased to 24x12 inches
################################################################################

cat("\n=== STARTING ANALYSIS 1: SPATIAL MAPS AND BOXPLOTS (CORRECTED) ===\n")

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

cat("\n--- Calculating cluster contributions within each quartile (by year) ---\n")
cat("Using within_quartile_prop variable\n\n")

# Calculate cluster proportions within each quartile for each year
cluster_yearly_production <- mgmt_data %>%
  group_by(year, cluster, quartile_clean) %>%
  summarise(
    cluster_proportion = sum(within_quartile_prop, na.rm = TRUE),
    .groups = "drop"
  )

# Verify that proportions sum to ~1 (100%) for each year-quartile
cat("\n--- VERIFICATION: Checking that cluster proportions sum to 1.0 within each year-quartile ---\n")
verification <- cluster_yearly_production %>%
  group_by(year, quartile_clean) %>%
  summarise(
    sum_proportions = sum(cluster_proportion, na.rm = TRUE),
    n_clusters = n(),
    .groups = "drop"
  )

cat("Range of sums:", range(verification$sum_proportions), "\n")
cat("All sums close to 1.0?", all(abs(verification$sum_proportions - 1.0) < 0.01), "\n\n")

if (!all(abs(verification$sum_proportions - 1.0) < 0.01)) {
  cat("⚠ WARNING: Cluster proportions do NOT sum to 1.0 for all year-quartiles!\n")
  cat("Showing cases where sum != 1.0:\n")
  print(verification %>% filter(abs(sum_proportions - 1.0) >= 0.01))
  cat("\n")
}

# Average across years to get mean proportions for mapping
cluster_avg_production <- cluster_yearly_production %>%
  group_by(cluster, quartile_clean) %>%
  summarise(
    production_proportion = mean(cluster_proportion, na.rm = TRUE),
    sd_proportion = sd(cluster_proportion, na.rm = TRUE),
    n_years = n(),
    .groups = "drop"
  ) %>%
  mutate(cluster = factor(cluster, levels = cluster_order))

cat("--- Cluster average proportional contribution by quartile (averaged across years) ---\n")
print(cluster_avg_production)

# Final verification: Show total contribution by quartile
cat("\n--- Sum of average proportions by quartile (should be ~1.0) ---\n")
quartile_sums <- cluster_avg_production %>%
  group_by(quartile_clean) %>%
  summarise(total_proportion = sum(production_proportion))
print(quartile_sums)

# Prepare boxplot data showing year-to-year variability
# Convert proportions to percentages
quartile_boxplot_data <- cluster_yearly_production %>%
  mutate(
    within_quartile_pct = cluster_proportion * 100,
    cluster = factor(cluster, levels = rev(cluster_order)),
    # Add display labels with "Group" instead of "Cluster"
    cluster_label = factor(
      gsub("Cluster", "Group", cluster),
      levels = gsub("Cluster", "Group", rev(cluster_order))
    )
  )

# --- CORRECTED: Calculate MEDIAN values for spatial mapping ---
cat("\n--- Calculating MEDIAN cluster values for spatial mapping ---\n")

cluster_median_values <- quartile_boxplot_data %>%
  group_by(cluster, quartile_clean) %>%
  summarise(
    median_production_pct = median(within_quartile_pct, na.rm = TRUE),
    .groups = "drop"
  )

cat("Cluster median values by quartile:\n")
print(cluster_median_values)

# --- Load spatial data ---

edges <- st_read(SPATIAL_PATH, quiet = TRUE)
basin <- st_read(BASIN_PATH, quiet = TRUE)

# Add cluster assignments to spatial data
edges_with_clusters <- edges %>%
  left_join(cluster_assignments, by = "mgmt_river")

# --- Prepare data using CLUSTER MEDIAN values ---
all_quartile_data <- list()
for (q in c("Q1", "Q2", "Q3", "Q4")) {
  
  # Get median value for each cluster in this quartile
  quartile_medians <- cluster_median_values %>% 
    filter(quartile_clean == q)
  
  edges_with_data <- edges_with_clusters %>%
    left_join(quartile_medians, by = "cluster") %>%
    filter(!is.na(cluster)) %>%
    mutate(
      production_proportion = median_production_pct / 100,  # Convert back to proportion for color scale
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

# Determine appropriate color scale limits based on actual data
color_limits <- c(0, max(combined_edges$production_proportion, na.rm = TRUE) * 1.1)
cat("Color scale range:", round(color_limits[1], 3), "to", round(color_limits[2], 3), "\n")

cat("\n--- Creating Figure 3a: contribution by group, faceted (free y) ---\n")

# Group color palette (the same hues used for each group throughout the analysis)
group_levels <- gsub("Cluster", "Group", cluster_order)        # "Group 1" ... "Group 6"
group_colors <- setNames(cluster_colors, group_levels)         # named by Group label

# Data for Figure 3a: x = seasonal quartile, y = within-quartile contribution (%),
# one facet per group, boxplots capturing year-to-year variability.
fig3a_data <- quartile_boxplot_data %>%
  mutate(
    quartile_clean = factor(quartile_clean, levels = c("Q1", "Q2", "Q3", "Q4")),
    group_label    = factor(gsub("Cluster", "Group", cluster), levels = group_levels)
  )

library(patchwork)

# Modern look: slim translucent boxes (no whisker staples, no separate outlier
# glyphs) with the raw yearly points overlaid (un-jittered) in the group color.
# Font sizes are deliberately large so the panel reads in a combined/publication
# figure where each sub-plot is scaled down.
fig3a <- ggplot(fig3a_data,
                aes(x = quartile_clean, y = within_quartile_pct)) +
  geom_boxplot(
    aes(fill = group_label, color = group_label),
    width = 0.55, alpha = 0.28, linewidth = 0.8,
    outlier.shape = NA, staplewidth = 0, fatten = 1.6
  ) +
  geom_point(
    aes(color = group_label),
    size = 2.2, alpha = 0.8, stroke = 0
  ) +
  facet_wrap(~ group_label, scales = "free_y", ncol = 3) +   # 3 columns x 2 rows
  scale_fill_manual(values = group_colors, guide = "none") +
  scale_color_manual(values = group_colors, guide = "none") +
  scale_y_continuous(labels = function(y) paste0(y, "%"),
                     expand = expansion(mult = c(0.05, 0.10))) +
  labs(x = "Seasonal quartile", y = "Proportional contribution within quartile") +
  theme_minimal(base_size = 20, base_family = "sans") +
  theme(
    strip.text         = element_text(size = 24, face = "bold", color = "grey20", hjust = 0),
    strip.background   = element_blank(),
    axis.title         = element_text(size = 24, face = "bold", color = "grey25"),
    axis.title.x       = element_text(margin = margin(t = 14)),
    axis.title.y       = element_text(margin = margin(r = 14)),
    axis.text.x        = element_text(size = 20, color = "grey30"),
    axis.text.y        = element_text(size = 18, color = "grey35"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey92", linewidth = 0.6),
    panel.spacing      = unit(1.8, "lines"),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(16, 18, 16, 16)
  )

ggsave(file.path(OUTPUT_DIR, paste0("Figure3a_contribution_by_group_faceted_", year_range, ".png")),
       fig3a, width = 15, height = 9, dpi = 300, bg = "white")

cat("✓ Figure 3a (faceted by group) saved\n")
cat("  Facets: one per group (Group 1-6), 3 columns x 2 rows, free y-axes\n")
cat("  X-axis: seasonal quartiles (Q1-Q4); boxplots + yearly points (un-jittered)\n")
cat("  Colors: per-group palette\n\n")


# ----------------------------------------------------------------------------
# FIGURE 3b: stacked bar of each group's contribution to the full annual run
# ----------------------------------------------------------------------------
cat("\n--- Creating Figure 3b: stacked group contribution to annual run by year ---\n")

# Each group's share of the full year's run = total_run_prop summed over the
# group's rivers across all quartiles, renormalized to sum to 100% per year
# across the clustered groups (matching the within-quartile group set).
fig3b_data <- mgmt_data %>%
  group_by(year, cluster) %>%
  summarise(group_run_prop = sum(total_run_prop, na.rm = TRUE), .groups = "drop") %>%
  group_by(year) %>%
  mutate(group_pct = 100 * group_run_prop / sum(group_run_prop)) %>%
  ungroup() %>%
  mutate(group_label = factor(gsub("Cluster", "Group", cluster), levels = group_levels))

fig3b <- ggplot(fig3b_data,
                aes(x = factor(year), y = group_pct, fill = group_label)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.4,
           position = position_stack(reverse = FALSE)) +
  scale_fill_manual(values = group_colors, breaks = group_levels, name = NULL) +
  scale_y_continuous(labels = function(y) paste0(y, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Year", y = "Proportional contribution to annual run") +
  theme_minimal(base_size = 20, base_family = "sans") +
  theme(
    axis.title         = element_text(size = 24, face = "bold", color = "grey20"),
    axis.title.x       = element_text(margin = margin(t = 12)),
    axis.title.y       = element_text(margin = margin(r = 12)),
    axis.text          = element_text(size = 20, color = "grey25"),
    legend.text        = element_text(size = 20, color = "grey20"),
    legend.key.size    = unit(1.3, "lines"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.6),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(16, 16, 16, 16)
  )

ggsave(file.path(OUTPUT_DIR, paste0("Figure3b_group_contribution_stacked_", year_range, ".png")),
       fig3b, width = 10, height = 7, dpi = 300, bg = "white")

cat("✓ Figure 3b (stacked group contribution by year) saved\n")
cat("  X-axis: year; stacked bars sum to 100% per year\n")
cat("  Segments colored by per-group palette\n\n")


# ----------------------------------------------------------------------------
# FIGURE 3b_Q1: stacked bar of each group's contribution to the Q1 run only
# ----------------------------------------------------------------------------
# Identical to Figure 3b, but restricted to the first seasonal quartile (Q1).
# Each group's share of the Q1 run = total_run_prop summed over the group's
# rivers within Q1, renormalized to sum to 100% per year across the clustered
# groups.
cat("\n--- Creating Figure 3b_Q1: stacked group contribution to Q1 run by year ---\n")

fig3b_q1_data <- mgmt_data %>%
  filter(quartile_clean == "Q1") %>%
  group_by(year, cluster) %>%
  summarise(group_run_prop = sum(total_run_prop, na.rm = TRUE), .groups = "drop") %>%
  group_by(year) %>%
  mutate(group_pct = 100 * group_run_prop / sum(group_run_prop)) %>%
  ungroup() %>%
  mutate(group_label = factor(gsub("Cluster", "Group", cluster), levels = group_levels))

fig3b_q1 <- ggplot(fig3b_q1_data,
                aes(x = factor(year), y = group_pct, fill = group_label)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.4,
           position = position_stack(reverse = FALSE)) +
  scale_fill_manual(values = group_colors, breaks = group_levels, name = NULL) +
  scale_y_continuous(labels = function(y) paste0(y, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Year", y = "Proportional contribution to Q1 run") +
  theme_minimal(base_size = 20, base_family = "sans") +
  theme(
    axis.title         = element_text(size = 24, face = "bold", color = "grey20"),
    axis.title.x       = element_text(margin = margin(t = 12)),
    axis.title.y       = element_text(margin = margin(r = 12)),
    axis.text          = element_text(size = 20, color = "grey25"),
    legend.text        = element_text(size = 20, color = "grey20"),
    legend.key.size    = unit(1.3, "lines"),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.6),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(16, 16, 16, 16)
  )

ggsave(file.path(OUTPUT_DIR, paste0("Figure3b_Q1_group_contribution_stacked_", year_range, ".png")),
       fig3b_q1, width = 10, height = 7, dpi = 300, bg = "white")

cat("✓ Figure 3b_Q1 (stacked Q1 group contribution by year) saved\n")
cat("  X-axis: year; stacked bars sum to 100% of the Q1 run per year\n")
cat("  Segments colored by per-group palette\n\n")


# ----------------------------------------------------------------------------
# FIGURE 3_QuartileCont: combined one-row figure (3a | 3b) for the manuscript
# ----------------------------------------------------------------------------
cat("\n--- Creating Figure 3_QuartileCont: combined 3a + 3b (one row) ---\n")

fig3_combined <- (fig3a | fig3b_q1) +
  plot_layout(widths = c(1.7, 1)) +              # 3a (6 facets) wider than 3b_q1
  plot_annotation(
    tag_levels = list(c("A", "B")),
    theme = theme(plot.background = element_rect(fill = "white", color = NA))
  ) &
  theme(plot.tag = element_text(size = 30, face = "bold", color = "grey15"))

ggsave(file.path(OUTPUT_DIR, paste0("Figure3_QuartileCont_", year_range, ".png")),
       fig3_combined, width = 26, height = 9.5, dpi = 300, bg = "white")

cat("✓ Figure 3_QuartileCont (combined one-row figure) saved\n")
cat("  Panel A: faceted boxplots by group; Panel B: stacked Q1 contribution by year\n")
cat("  Large labels/axes sized for publication\n\n")

# ==================== EXPORT ANALYSIS 1 DATA ====================
cat("--- Exporting Analysis 1 data to CSV ---\n")

# Export cluster average production by quartile (averaged across years)
write_csv(cluster_avg_production, 
          file.path(CSV_DIR, "analysis1_cluster_avg_production_by_quartile.csv"))

# Export yearly cluster proportions (for verification that they sum to 1)
write_csv(cluster_yearly_production %>%
            select(year, cluster, quartile_clean, cluster_proportion) %>%
            arrange(year, quartile_clean, cluster),
          file.path(CSV_DIR, "analysis1_cluster_proportions_by_year.csv"))

# Export boxplot data (raw data with all years)
write_csv(quartile_boxplot_data %>% 
            select(year, cluster, quartile_clean, within_quartile_pct, cluster_proportion),
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

# Export cluster median values used for spatial mapping
write_csv(cluster_median_values,
          file.path(CSV_DIR, "analysis1_cluster_median_values_for_maps.csv"))

cat("=== ANALYSIS 1 COMPLETE ===\n")
cat("Created: Figure 3a (faceted by group, free y) + Figure 3b (stacked by year)\n")
cat("Exported: 5 CSV files\n")
cat("FIGURE DESIGN:\n")
cat("  - Figure 3a: facet per group (3 x 2), x = seasonal quartile, free y-axes, group colors\n")
cat("  - Figure 3b: stacked bar per year of each group's share of the annual run\n\n")


################################################################################
# ANALYSIS 2: CUMULATIVE DISTRIBUTION OVER TIME (WITH CPUE CURVES)
################################################################################
# MODIFICATION: Adds black curves showing remaining CPUE proportion
# Creates 6-panel figure (one per year) with both cluster cumulative distributions
# and CPUE curves overlaid
################################################################################

cat("=== STARTING ANALYSIS 2: CUMULATIVE DISTRIBUTION (CORRECTED) ===\n")

# Define date range: June 1st (DOY 152) to July 23rd (DOY 205)
MIN_DOY <- 149  # June 1st
MAX_DOY <- 205  # July 23rd

# Load cumulative distribution data
cumulative_data <- read_csv(CUMULATIVE_DATA_PATH) %>%
  filter(mgmt_river != "Johnson", mgmt_river != "Lower Kusko") %>%
  left_join(cluster_assignments, by = "mgmt_river") %>%
  filter(!is.na(cluster)) %>%
  filter(doy >= MIN_DOY & doy <= MAX_DOY)  # Trim to June 1st - July 23rd

cat("✓ Loaded cumulative distribution data\n")
cat("  Data trimmed to DOY", MIN_DOY, "-", MAX_DOY, "(June 1st - July 23rd)\n")
cat("  Years in cumulative data:", paste(sort(unique(cumulative_data$year)), collapse = ", "), "\n")

# Load CPUE data and filter to only years present in cumulative data
CPUE_DATA_PATH <- DAILY_CPUE_CSV             # written by R/04_cpue_remaining.R

cpue_data_full <- read_csv(CPUE_DATA_PATH)
cat("  Years in CPUE data (before filtering):", paste(sort(unique(cpue_data_full$year)), collapse = ", "), "\n")

# Get years available in cumulative data
available_years <- unique(cumulative_data$year)

cpue_data <- cpue_data_full %>%
  filter(year %in% available_years) %>%  # Only keep years that match cluster data
  filter(doy >= MIN_DOY & doy <= MAX_DOY) %>%  # Trim to June 1st - July 23rd
  mutate(remaining_proportion_pct = remaining_proportion * 100)

cat("  Years in CPUE data (after filtering):", paste(sort(unique(cpue_data$year)), collapse = ", "), "\n")
cat("  DOY range after trimming:", min(cumulative_data$doy), "to", max(cumulative_data$doy), "\n\n")

if (length(setdiff(unique(cpue_data_full$year), available_years)) > 0) {
  excluded_years <- setdiff(unique(cpue_data_full$year), available_years)
  cat("  NOTE: Excluding CPUE data for years:", paste(excluded_years, collapse = ", "), "\n")
  cat("        (these years not present in cluster cumulative data)\n\n")
}

################################################################################
# SECTION 2.1: DATA AGGREGATION (CORRECTED)
################################################################################

cat("--- Aggregating management units to cluster level (CORRECTED METHOD) ---\n")

cluster_annual_totals <- cumulative_data %>%
  group_by(year, cluster, mgmt_river) %>%
  summarise(final_production = max(cumulative_production, na.rm = TRUE), .groups = "drop") %>%
  group_by(year, cluster) %>%
  summarise(cluster_annual_total = sum(final_production, na.rm = TRUE), .groups = "drop")

# CORRECTED AGGREGATION: Calculate interval production first, then aggregate
cluster_aggregated <- cumulative_data %>%
  # First, calculate production at each interval for each unit
  group_by(year, mgmt_river, cluster) %>%
  arrange(doy) %>%
  mutate(interval_production = cumulative_production - lag(cumulative_production, default = 0)) %>%
  ungroup() %>%
  
  # Then sum interval production across units within each cluster
  group_by(year, doy, cluster) %>%
  summarise(cluster_interval_production = sum(interval_production, na.rm = TRUE), .groups = "drop") %>%
  
  # Now calculate cumulative production for the cluster
  group_by(year, cluster) %>%
  arrange(doy) %>%
  mutate(cluster_cumulative_production = cumsum(cluster_interval_production)) %>%
  ungroup() %>%
  
  # Join with annual totals and calculate percentages
  left_join(cluster_annual_totals, by = c("year", "cluster")) %>%
  mutate(
    cluster_cumulative_percent = (cluster_cumulative_production / cluster_annual_total) * 100,
    cluster = factor(cluster, levels = cluster_order)
  )

# Define color scheme (should already be defined earlier, but include here for completeness)
if (!exists("cluster_colors")) {
  cluster_colors <- colorRampPalette(c("#8B0000", "#FF0000", "#FF8C00", "#87CEEB", "#1E90FF", "#000080"))(6)
  names(cluster_colors) <- cluster_order
}

years <- sort(unique(cluster_aggregated$year))
cat("Years to plot:", paste(years, collapse = ", "), "\n")
cat("Aggregation method: Interval-based (CORRECTED)\n\n")


################################################################################
# SECTION 2.2: INDIVIDUAL YEAR PLOTS (DISABLED - NOT EXPORTED)
################################################################################

cat("--- Skipping individual year plots (only creating combined figure) ---\n\n")


################################################################################
# SECTION 2.3: AVERAGE PLOT WITH STANDARD ERROR BANDS (DISABLED - NOT EXPORTED)
################################################################################

cat("--- Skipping average plot (only creating combined figure) ---\n\n")

# Calculate average data (still computed for CSV export)
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


################################################################################
# SECTION 2.4: MULTI-PANEL PLOT WITH ALL YEARS (MODIFIED WITH CPUE)
################################################################################

cat("--- Creating multi-panel plot with all years and CPUE curves ---\n")

# Calculate plot height based on number of years
n_years <- length(years)
plot_height <- max(15, n_years * 2.5)

p_multi <- ggplot() +
  # Horizontal reference lines
  geom_hline(yintercept = c(25, 50, 75), color = "gray70", linetype = "dashed", alpha = 0.5) +
  
  # Vertical line at DOY 162 (June 11)
  geom_vline(xintercept = 162, color = "grey50", linetype = "dashed", linewidth = 1.2, alpha = 0.8) +
  
  # Cluster cumulative distribution lines (colored)
  geom_line(data = cluster_aggregated,
            aes(x = doy, y = cluster_cumulative_percent, color = cluster),
            linewidth = 1.2, alpha = 0.9) +
  
  geom_point(data = cluster_aggregated,
             aes(x = doy, y = cluster_cumulative_percent, color = cluster),
             size = 0.8, alpha = 0.7) +
  
  # CPUE remaining proportion curve (black)
  geom_line(data = cpue_data,
            aes(x = doy, y = remaining_proportion_pct),
            color = "grey40", linewidth = 1.2, alpha = 0.6) +
  
  # Scales
  scale_color_manual(
    values = cluster_colors, 
    breaks = cluster_order,
    labels = function(x) gsub("Cluster", "Group", x),
    guide = guide_legend(override.aes = list(linewidth = 2))
  ) +
  
  scale_y_continuous(
    labels = function(x) paste0(round(x, 0), "%"),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  
  scale_x_continuous(
    breaks = seq(152, 205, by = 10),  # June 1st (152) to July 23rd (205)
    labels = function(x) {
      date_str <- format(as.Date(x - 1, origin = "2020-01-01"), "%b %d")
      date_str
    }
  ) +
  
  facet_wrap(~ year, ncol = 1, scales = "fixed", strip.position = "top") +
  
  labs(
    title = "Cumulative distribution vs. Time",
    x = "Date",
    y = "Cumulative Percent of Total Relative Abundance by Group",
    color = "Group"
  ) +
  
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
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
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14, color = "grey20", face = "bold"),
    axis.text.y = element_text(size = 14, color = "grey20", face = "bold"),
    axis.title.x = element_text(size = 16, color = "grey20", face = "bold"),
    axis.title.y = element_text(size = 16, color = "grey20", face = "bold")
  )

ggsave(
  file.path(OUTPUT_DIR, "cluster_cumulative_distribution_MULTIPANEL.png"),
  p_multi, 
  width = 7, 
  height = plot_height, 
  dpi = 300, 
  bg = "white",
  limitsize = FALSE
)

cat("Created multi-panel plot with all years and CPUE curves\n")
cat("  Plot dimensions: 10 x", plot_height, "inches\n")
cat("  Number of panels:", n_years, "\n")
cat("  X-axis: Fixed scale across all panels\n")
cat("  Axis labels: Larger and bolder\n")
cat("  Black curves: Remaining CPUE proportion\n\n")


# ==================== EXPORT ANALYSIS 2 DATA ====================
cat("--- Exporting Analysis 2 data to CSV ---\n")

# Export cluster-aggregated cumulative data (all years)
write_csv(cluster_aggregated %>%
            select(year, doy, cluster, cluster_interval_production, 
                   cluster_cumulative_production, cluster_annual_total, 
                   cluster_cumulative_percent),
          file.path(CSV_DIR, "analysis2_cluster_cumulative_by_year.csv"))

# Export average cumulative data with SE
write_csv(avg_cumulative,
          file.path(CSV_DIR, "analysis2_average_cumulative_with_se.csv"))

# Export CPUE data (for reference)
write_csv(cpue_data,
          file.path(CSV_DIR, "analysis2_cpue_remaining_proportion.csv"))

cat("=== ANALYSIS 2 COMPLETE ===\n")
cat("Created: 1 multi-panel plot with CPUE curves\n")
cat("Exported: 3 CSV files\n")
cat("MODIFICATIONS:\n")
cat("  - Added black CPUE remaining proportion curves\n")
cat("  - Maintains original multi-panel (1 column) format\n")
cat("  - All other styling unchanged\n\n")


################################################################################
# ANALYSIS 3: DAYS TO 50% RETURN BY CLUSTER (WITH CIRCLE MARKER)
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

# --- Create horizontal range plot with circle marker (MODIFIED) ---

p_days_to_50 <- ggplot(days_to_50_summary, aes(y = cluster)) +
  # Range lines (min to max)
  geom_segment(aes(x = min_days, xend = max_days, color = cluster), 
               linewidth = 6, alpha = 0.7) +
  
  # MODIFIED: Mean markers (white-filled circles instead of black lines)
  geom_point(aes(x = mean_days, fill = cluster), 
             color = "white", size = 4, shape = 21, stroke = 2) +
  
  # Colors
  scale_color_manual(values = cluster_colors, breaks = cluster_order, guide = "none") +
  scale_fill_manual(values = cluster_colors, breaks = cluster_order, guide = "none") +
  
  scale_x_continuous(
    name = "Days to 50% return",
    breaks = function(x) pretty(x, n = 6),
    expand = expansion(mult = c(0.02, 0.02)),
    labels = function(x) paste0(x, " days")
  ) +
  
  scale_y_discrete(
    name = NULL,
    limits = rev(cluster_order),
    labels = function(x) gsub("Cluster", "Group", x)
  ) +
  
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

# [trimmed: not a paper figure] days-to-50% plot not exported (CSVs below kept).
# ggsave(file.path(OUTPUT_DIR, "cluster_days_to_50_percent.png"),
#        p_days_to_50, width = 10, height = 6, dpi = 300, bg = "white")

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
# ANALYSIS 4: FRONT-END CLOSURE PROTECTION (MODIFIED - POINT PLOT)
################################################################################
# Uses ACTUAL CLOSURE WINDOW DATES (not Q1) to match management unit analysis
# Default: DOY 153-162 (approximately June 1-11)
# MODIFICATION: Shows individual points for each year instead of violin plot
################################################################################

# Define cluster colors (if not already defined earlier in script)
if (!exists("cluster_colors")) {
  cluster_colors <- colorRampPalette(c("#8B0000", "#FF0000", "#FF8C00", "#87CEEB", "#1E90FF", "#000080"))(6)
  names(cluster_colors) <- cluster_order
}

cat("=== STARTING ANALYSIS 4: FRONT-END CLOSURE PROTECTION (MODIFIED) ===\n")
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
cat("\n--- Closure Protection by Group ---\n")
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
# STEP 4.3: CREATE VISUALIZATION (MODIFIED - POINT PLOT WITH RANGE BARS)
################################################################################

cat("\n--- Creating visualization (point plot with range bars and median) ---\n")

# Calculate summary statistics for range bars and median
closure_summary_plot <- cluster_closure_protection %>%
  group_by(cluster) %>%
  summarise(
    min_pct = min(closure_protection_pct, na.rm = TRUE),
    max_pct = max(closure_protection_pct, na.rm = TRUE),
    median_pct = median(closure_protection_pct, na.rm = TRUE),
    .groups = "drop"
  )

p_closure <- ggplot(cluster_closure_protection, 
                    aes(y = cluster, x = closure_protection_pct)) +
  # Add semi-transparent range bars underneath (min to max)
  geom_segment(data = closure_summary_plot,
               aes(x = min_pct, xend = max_pct, y = cluster, yend = cluster,
                   color = cluster),
               linewidth = 8, alpha = 0.3) +
  
  # Show individual points for each year (colored by cluster)
  geom_point(aes(color = cluster, fill = cluster), 
             size = 4, alpha = 0.8) +
  
  # Add black vertical line for median
  geom_segment(data = closure_summary_plot,
               aes(x = median_pct, xend = median_pct, 
                   y = as.numeric(cluster) - 0.3, yend = as.numeric(cluster) + 0.3),
               color = "black", linewidth = 1.2) +
  
  scale_color_manual(values = rev(cluster_colors), guide = "none") +
  scale_fill_manual(values = rev(cluster_colors), guide = "none") +
  scale_x_continuous(
    name = "Protection Effectiveness (%)",
    labels = function(x) paste0(round(x, 1), "%"),
    expand = expansion(mult = c(0.02, 0.05)),
    limits = c(0, NA)
  ) +
  scale_y_discrete(
    name = NULL,
    labels = function(x) gsub("Cluster", "Group", x)
  ) +
  labs(
    title = "Front-End Closure Protection by Group") +
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

# [trimmed: superseded by Fig 5 offset comparison] single-window closure plot
# not exported. The cluster_closure_protection computation above is still used
# by Analysis 5 / 6, so nothing downstream changes.
# ggsave(file.path(OUTPUT_DIR, "cluster_front_end_closure_protection.png"),
#        p_closure, width = 12, height = 10, dpi = 300, bg = "white")

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
cat("Created: 1 point plot showing front-end closure protection\n")
cat("Exported: 3 CSV files\n\n")

################################################################################
# ANALYSIS 5: FRONT-END CLOSURE TIME OFFSET ANALYSIS
################################################################################
# PURPOSE: Analyze protection percentages at different time offsets from the
#          current front end closure window (+3, +6, -3 days)
# USES: Same calculation methodology as Analysis 4 (front end closure protection)
# OUTPUT: Horizontal comparison plot showing protection % by cluster for all time offsets
################################################################################

# Load required library for pivot_wider
if (!require("tidyr", quietly = TRUE)) {
  library(tidyr)
}

cat("=== STARTING ANALYSIS 5: CLOSURE TIME OFFSET ANALYSIS ===\n")
cat("Base closure window: DOY", CLOSURE_START_DOY, "to", CLOSURE_END_DOY, "\n")
cat("Analyzing offsets: +3 days, +6 days, -3 days\n\n")

################################################################################
# DEFINE TIME OFFSETS
################################################################################

# Define the time offsets to analyze
# Note: 0 days corresponds to June 11 (DOY 162, the end of the current closure)
time_offsets <- list(
  "minus_3" = list(offset = -3, label = "-3 days (June 8)", doy_label = paste0("DOY ", CLOSURE_END_DOY - 3)),
  "plus_3" = list(offset = 3, label = "+3 days (June 14)", doy_label = paste0("DOY ", CLOSURE_END_DOY + 3)),
  "plus_6" = list(offset = 6, label = "+6 days (June 17)", doy_label = paste0("DOY ", CLOSURE_END_DOY + 6))
)

################################################################################
# FUNCTION: CALCULATE CLOSURE PROTECTION FOR A GIVEN TIME OFFSET
################################################################################

calculate_offset_protection <- function(cumulative_data, offset_days) {
  # Calculate adjusted closure window
  adj_start_doy <- CLOSURE_START_DOY + offset_days
  adj_end_doy <- CLOSURE_END_DOY + offset_days
  
  cat(paste0("  Calculating for offset window: DOY ", adj_start_doy, " to ", adj_end_doy, "\n"))
  
  # Step 1: Calculate closure window production for each unit
  offset_protection_by_unit <- cumulative_data %>%
    group_by(year, mgmt_river, cluster) %>%
    arrange(doy) %>%
    summarise(
      # Total annual production (same for all offsets)
      total_annual_production = max(cumulative_production, na.rm = TRUE),
      
      # Production during OFFSET closure window
      closure_window_production = {
        # Find production at the last day of offset closure window (or closest day)
        closure_days <- doy[doy >= adj_start_doy & doy <= adj_end_doy]
        if (length(closure_days) > 0) {
          max_closure_doy <- max(closure_days)
          cumulative_production[doy == max_closure_doy][1]
        } else {
          0  # No fish during closure window
        }
      },
      
      # Check if we have data for closure window
      has_closure_data = any(doy >= adj_start_doy & doy <= adj_end_doy),
      
      .groups = "drop"
    ) %>%
    mutate(
      # Handle cases where closure window production is NA
      closure_window_production = ifelse(is.na(closure_window_production), 0, closure_window_production),
      
      # Calculate protection percentage for this unit
      unit_protection_pct = (closure_window_production / total_annual_production) * 100
    )
  
  # Step 2: Aggregate to cluster level
  cluster_offset_protection <- offset_protection_by_unit %>%
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
  
  return(list(
    by_unit = offset_protection_by_unit,
    by_cluster = cluster_offset_protection
  ))
}

################################################################################
# CALCULATE PROTECTION FOR EACH TIME OFFSET
################################################################################

cat("\n--- Calculating protection for each time offset ---\n")

# Store results for each offset
offset_results <- list()

for (offset_name in names(time_offsets)) {
  offset_info <- time_offsets[[offset_name]]
  cat(paste0("\nProcessing ", offset_info$label, " offset...\n"))
  
  # Calculate protection for this offset
  results <- calculate_offset_protection(cumulative_data, offset_info$offset)
  
  # Add offset information to the results
  results$by_cluster <- results$by_cluster %>%
    mutate(offset_label = offset_info$label)
  
  results$by_unit <- results$by_unit %>%
    mutate(offset_label = offset_info$label)
  
  # Store results
  offset_results[[offset_name]] <- results
  
  # Print summary
  summary_stats <- results$by_cluster %>%
    group_by(cluster) %>%
    summarise(
      mean_pct = round(mean(closure_protection_pct), 1),
      min_pct = round(min(closure_protection_pct), 1),
      max_pct = round(max(closure_protection_pct), 1),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_pct))
  
  cat(paste0("\n  Summary for ", offset_info$label, " offset:\n"))
  print(summary_stats)
}

################################################################################
# CREATE COMBINED DATA FOR PLOTTING
################################################################################

cat("\n--- Combining data for plotting ---\n")

# Combine all offset data plus the original closure data
combined_data <- bind_rows(
  # Original closure data (offset = 0)
  cluster_closure_protection %>%
    mutate(offset_label = "Current (June 11)",
           offset_order = 2),
  
  # -3 days offset
  offset_results$minus_3$by_cluster %>%
    mutate(offset_order = 1),
  
  # +3 days offset
  offset_results$plus_3$by_cluster %>%
    mutate(offset_order = 3),
  
  # +6 days offset
  offset_results$plus_6$by_cluster %>%
    mutate(offset_order = 4)
) %>%
  mutate(
    offset_label = factor(offset_label, 
                          levels = c("-3 days (June 8)", "Current (June 11)", 
                                     "+3 days (June 14)", "+6 days (June 17)"))
  )

################################################################################
# CREATE HORIZONTAL COMPARISON PLOT
################################################################################

cat("\n--- Creating horizontal comparison plot ---\n")

# Use cluster colors (already defined earlier in script)
# cluster_colors should already be defined with 6 colors for Cluster 1-6

# Calculate summary statistics for plotting (min, max, median)
horizontal_plot_summary <- combined_data %>%
  group_by(cluster, offset_label) %>%
  summarise(
    min_pct = min(closure_protection_pct, na.rm = TRUE),
    max_pct = max(closure_protection_pct, na.rm = TRUE),
    median_pct = median(closure_protection_pct, na.rm = TRUE),
    .groups = "drop"
  )

# Add position offsets for horizontal separation of points
horizontal_plot_data <- combined_data %>%
  mutate(
    cluster_num = as.numeric(factor(cluster, levels = cluster_order)),
    offset_position = case_when(
      offset_label == "-3 days (June 8)" ~ cluster_num - 0.2,
      offset_label == "Current (June 11)" ~ cluster_num - 0.07,
      offset_label == "+3 days (June 14)" ~ cluster_num + 0.07,
      offset_label == "+6 days (June 17)" ~ cluster_num + 0.2
    )
  )

horizontal_summary_positioned <- horizontal_plot_summary %>%
  mutate(
    cluster_num = as.numeric(factor(cluster, levels = cluster_order)),
    offset_position = case_when(
      offset_label == "-3 days (June 8)" ~ cluster_num - 0.2,
      offset_label == "Current (June 11)" ~ cluster_num - 0.07,
      offset_label == "+3 days (June 14)" ~ cluster_num + 0.07,
      offset_label == "+6 days (June 17)" ~ cluster_num + 0.2
    )
  )

# Create a mapping from cluster to color
cluster_color_mapping <- setNames(cluster_colors, cluster_order)

# Create horizontal comparison plot
p_comparison <- ggplot() +
  # Add transparent range bars (min to max) colored by cluster
  geom_segment(data = horizontal_summary_positioned,
               aes(x = offset_position, xend = offset_position,
                   y = min_pct, yend = max_pct,
                   color = cluster),
               linewidth = 1.5, alpha = 0.25) +
  
  # Add individual year points (colored by cluster)
  geom_point(data = horizontal_plot_data,
             aes(x = offset_position, y = closure_protection_pct,
                 color = cluster, fill = cluster),
             size = 3.5, alpha = 0.5, shape = 21, stroke = 0.8) +
  
  # Add black median markers (horizontal lines)
  geom_segment(data = horizontal_summary_positioned,
               aes(x = offset_position - 0.05, xend = offset_position + 0.05,
                   y = median_pct, yend = median_pct),
               color = "black", linewidth = 2.5, alpha = 1) +
  
  # Add scenario labels below each set of points
  geom_text(data = data.frame(
    cluster_num = rep(1:6, each = 4),
    offset_position = rep(c(0.8, 0.93, 1.07, 1.2), 6) + rep(0:5, each = 4),
    label = rep(c("-3", "0", "+3", "+6"), 6)
  ),
  aes(x = offset_position, y = -2, label = label),
  size = 3.5, color = "gray30", fontface = "bold") +
  
  scale_color_manual(values = cluster_color_mapping, name = "Group",
                     labels = function(x) gsub("Cluster", "Group", x)) +
  scale_fill_manual(values = cluster_color_mapping, name = "Group",
                    labels = function(x) gsub("Cluster", "Group", x)) +
  scale_x_continuous(
    name = "Group",
    breaks = 1:6,
    labels = paste0("Group ", 1:6),
    limits = c(0.5, 6.5),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    name = "Percent of Relative Abundance Protected",
    labels = function(x) paste0(round(x, 0), "%"),
    expand = expansion(mult = c(0.1, 0.05)),
    limits = c(-3, NA)
  ) +
  labs(
    title = "Front-End Closure Protection by Group and Timing Offset",
    subtitle = paste0("Comparison across -3, 0, +3, and +6 day offsets (", year_range, ")")
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "gray90", linewidth = 0.4),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    axis.text = element_text(size = 11, color = "gray20"),
    axis.title = element_text(face = "bold", size = 12, color = "gray20"),
    axis.line.x = element_line(color = "gray40", linewidth = 0.5),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 30, 10)
  )

ggsave(file.path(OUTPUT_DIR, "cluster_front_end_closure_protection_comparison.png"),
       p_comparison, width = 11, height = 7, dpi = 300, bg = "white")

cat("✓ Horizontal comparison plot saved\n")

################################################################################
# EXPORT ANALYSIS 5 DATA TO CSV
################################################################################

cat("\n--- Exporting Analysis 5 data to CSV ---\n")

# Export cluster-level data for each offset
for (offset_name in names(time_offsets)) {
  offset_info <- time_offsets[[offset_name]]
  cluster_data <- offset_results[[offset_name]]$by_cluster
  
  # Export cluster-level year-by-year data
  write_csv(cluster_data %>% 
              select(year, cluster, closure_protection_pct, n_units,
                     cluster_total_production, cluster_closure_production, offset_label),
            file.path(CSV_DIR, paste0("analysis5_cluster_closure_protection_by_year_",
                                      gsub(" ", "_", offset_info$label), ".csv")))
  
  # Export unit-level data
  unit_data <- offset_results[[offset_name]]$by_unit
  write_csv(unit_data %>% 
              select(year, mgmt_river, cluster, total_annual_production, 
                     closure_window_production, unit_protection_pct, offset_label),
            file.path(CSV_DIR, paste0("analysis5_unit_closure_protection_by_year_",
                                      gsub(" ", "_", offset_info$label), ".csv")))
}

# Export combined summary statistics across all offsets
combined_summary <- combined_data %>%
  group_by(cluster, offset_label, offset_order) %>%
  summarise(
    mean_pct = mean(closure_protection_pct, na.rm = TRUE),
    median_pct = median(closure_protection_pct, na.rm = TRUE),
    min_pct = min(closure_protection_pct, na.rm = TRUE),
    max_pct = max(closure_protection_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(offset_order, desc(mean_pct)) %>%
  select(cluster, offset_label, mean_pct, median_pct, min_pct, max_pct)

write_csv(combined_summary,
          file.path(CSV_DIR, "analysis5_closure_protection_summary_all_offsets.csv"))

# Export combined data for all offsets (long format)
write_csv(combined_data %>%
            select(year, cluster, closure_protection_pct, offset_label, n_units,
                   cluster_total_production, cluster_closure_production),
          file.path(CSV_DIR, "analysis5_cluster_closure_protection_all_offsets.csv"))

cat("=== ANALYSIS 5 COMPLETE ===\n")
cat("Created: 1 horizontal comparison plot\n")
cat("Exported: 9 CSV files\n\n")

################################################################################
# SUMMARY TABLE
################################################################################

cat("\n=== SUMMARY: PROTECTION PERCENTAGES BY OFFSET ===\n\n")

# Create a summary table showing mean protection by cluster for each offset
summary_table <- combined_summary %>%
  select(cluster, offset_label, mean_pct) %>%
  pivot_wider(names_from = offset_label, values_from = mean_pct) %>%
  arrange(cluster)

print(summary_table, n = Inf)

cat("\n✓ ANALYSIS 5 COMPLETE ✓\n\n")


################################################################################
# ANALYSIS 6: CPUE FOREGONE VS. PROTECTION EFFECTIVENESS
################################################################################
# PURPOSE: Plot the relationship between CPUE foregone and front-end closure
#          protection effectiveness for each cluster, year, and timing scenario
# CREATES: Combined scatter plot showing mean CPUE foregone % vs. mean protection %
#          with error bars showing ±1 standard deviation across years
################################################################################

cat("=== STARTING ANALYSIS 6: CPUE vs. PROTECTION ANALYSIS ===\n")

################################################################################
# STEP 6.1: GET CPUE FOREGONE AT CLOSURE END FOR EACH SCENARIO
################################################################################

cat("\n--- Calculating CPUE foregone at closure end for each scenario ---\n")

# Function to get CPUE foregone at a specific DOY
get_cpue_at_doy <- function(cpue_data, target_doy) {
  cpue_data %>%
    filter(doy == target_doy) %>%
    select(year, remaining_proportion_pct) %>%
    mutate(
      target_doy = target_doy,
      foregone_proportion_pct = 100 - remaining_proportion_pct  # Calculate foregone
    )
}

# Calculate CPUE foregone for each scenario
cpue_scenarios <- bind_rows(
  # -3 days (June 8): DOY 159 (162 - 3)
  get_cpue_at_doy(cpue_data, CLOSURE_END_DOY - 3) %>%
    mutate(offset_label = "-3 days (June 8)", offset_order = 1, offset_short = "-3 days"),
  
  # Current (June 11): DOY 162
  get_cpue_at_doy(cpue_data, CLOSURE_END_DOY) %>%
    mutate(offset_label = "Current (June 11)", offset_order = 2, offset_short = "0 days"),
  
  # +3 days (June 14): DOY 165 (162 + 3)
  get_cpue_at_doy(cpue_data, CLOSURE_END_DOY + 3) %>%
    mutate(offset_label = "+3 days (June 14)", offset_order = 3, offset_short = "+3 days"),
  
  # +6 days (June 17): DOY 168 (162 + 6)
  get_cpue_at_doy(cpue_data, CLOSURE_END_DOY + 6) %>%
    mutate(offset_label = "+6 days (June 17)", offset_order = 4, offset_short = "+6 days")
) %>%
  mutate(
    offset_label = factor(offset_label, 
                          levels = c("-3 days (June 8)", "Current (June 11)", 
                                     "+3 days (June 14)", "+6 days (June 17)"))
  )

cat("CPUE scenarios data:\n")
print(head(cpue_scenarios))
cat("Years in CPUE data:", paste(sort(unique(cpue_scenarios$year)), collapse = ", "), "\n")

################################################################################
# STEP 6.2: COMBINE CPUE WITH PROTECTION DATA
################################################################################

cat("\n--- Combining CPUE foregone with protection effectiveness ---\n")

# Use the combined_data from Analysis 5 which has all scenarios
protection_cpue_combined <- combined_data %>%
  left_join(cpue_scenarios, by = c("year", "offset_label")) %>%
  filter(!is.na(foregone_proportion_pct)) %>%  # Remove any years without CPUE data
  mutate(
    cluster = factor(cluster, levels = cluster_order)
  )

cat("Combined data rows:", nrow(protection_cpue_combined), "\n")
cat("Years with complete data:", paste(sort(unique(protection_cpue_combined$year)), collapse = ", "), "\n")
cat("Clusters:", paste(levels(protection_cpue_combined$cluster), collapse = ", "), "\n")

# Verify we have data for all scenarios
scenario_check <- protection_cpue_combined %>%
  group_by(offset_label) %>%
  summarise(
    n_observations = n(),
    n_years = n_distinct(year),
    n_clusters = n_distinct(cluster)
  )

cat("\nData availability by scenario:\n")
print(scenario_check)

################################################################################
# STEP 6.3: CALCULATE MEAN VALUES BY CLUSTER AND SCENARIO
################################################################################

cat("\n--- Calculating mean values across years ---\n")

# Calculate mean protection and CPUE for each cluster-scenario combination
mean_data <- protection_cpue_combined %>%
  group_by(cluster, offset_label) %>%
  summarise(
    n_years = n(),
    mean_cpue_foregone = mean(foregone_proportion_pct, na.rm = TRUE),
    sd_cpue_foregone = sd(foregone_proportion_pct, na.rm = TRUE),
    min_cpue_foregone = min(foregone_proportion_pct, na.rm = TRUE),
    max_cpue_foregone = max(foregone_proportion_pct, na.rm = TRUE),
    mean_protection = mean(closure_protection_pct, na.rm = TRUE),
    sd_protection = sd(closure_protection_pct, na.rm = TRUE),
    min_protection = min(closure_protection_pct, na.rm = TRUE),
    max_protection = max(closure_protection_pct, na.rm = TRUE),
    offset_short = first(offset_short),
    .groups = "drop"
  ) %>%
  mutate(
    cluster = factor(cluster, levels = cluster_order),
    # Add offset_order back based on offset_label
    offset_order = case_when(
      offset_label == "-3 days (June 8)" ~ 1,
      offset_label == "Current (June 11)" ~ 2,
      offset_label == "+3 days (June 14)" ~ 3,
      offset_label == "+6 days (June 17)" ~ 4,
      TRUE ~ NA_real_
    )
  ) %>%
  arrange(cluster, offset_order)

cat("Mean data summary:\n")
print(mean_data)

################################################################################
# STEP 6.4: CREATE COMBINED SCATTER PLOTS
################################################################################

cat("\n--- Creating combined scatter plots ---\n")

# Define scenario colors for consistency
scenario_colors <- c(
  "-3 days (June 8)" = "#E74C3C",    # Red
  "Current (June 11)" = "#3498DB",    # Blue
  "+3 days (June 14)" = "#2ECC71",    # Green
  "+6 days (June 17)" = "#F39C12"     # Orange
)

# Version 1: All clusters combined in one panel
cat("Creating combined plot (all groups in one panel)...\n")

p_combined <- ggplot(mean_data,
                     aes(x = mean_cpue_foregone, 
                         y = mean_protection,
                         color = cluster,
                         group = cluster)) +
  geom_line(linewidth = 1.5, alpha = 0.6) +
  geom_point(size = 6, alpha = 0.6,
             stroke = 2, shape = 21, aes(fill = cluster), color = "black") +
  geom_text(data = mean_data %>% 
              group_by(offset_label) %>% 
              filter(mean_protection == max(mean_protection)),  # Label topmost point in each scenario
            aes(label = offset_short), 
            vjust = -1.2, hjust = 0.5, size = 6, 
            show.legend = FALSE, fontface = "bold",
            color = "black") +
  scale_color_manual(values = cluster_colors, name = "Group",
                     labels = paste0("Group ", 1:6)) +
  scale_fill_manual(values = cluster_colors, name = "Group",
                    labels = paste0("Group ", 1:6)) +
  guides(color = guide_legend(title = "Group"),
         fill = guide_legend(title = "Group")) +
  scale_x_continuous(
    name = "Mean foregone harvest opportunity (%)",
    labels = function(x) paste0(round(x, 0), "%"),
    breaks = seq(0, 30, by = 5),
    limits = c(0, 30),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Mean Protection Effectiveness (%)",
    labels = function(x) paste0(round(x, 0), "%"),
    breaks = seq(0, 55, by = 10),
    limits = c(0, 55),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(face = "bold", size = 18),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 18),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(OUTPUT_DIR, "cpue_vs_protection_combined.png"),
       p_combined, width = 10, height = 10, dpi = 300, bg = "white")

cat("✓ Combined plot saved\n")

# Version 1-ALT: Identical to the combined plot above (same layout, legend
# position, and 10x10 dimensions) but with a grey dashed 1:1 reference line added.
cat("Creating combined ALT plot (same as original, with 1:1 line)...\n")

p_combined_alt <- ggplot(mean_data,
                         aes(x = mean_cpue_foregone,
                             y = mean_protection,
                             color = cluster,
                             group = cluster)) +
  geom_abline(slope = 1, intercept = 0,
              color = "grey60", linetype = "dashed", linewidth = 1) +
  geom_line(linewidth = 1.5, alpha = 0.6) +
  geom_point(size = 6, alpha = 0.6,
             stroke = 2, shape = 21, aes(fill = cluster), color = "black") +
  geom_text(data = mean_data %>%
              group_by(offset_label) %>%
              filter(mean_protection == max(mean_protection)),  # Label topmost-protection point in each scenario
            aes(label = offset_short),
            vjust = -1.2, hjust = 0.5, size = 6,
            show.legend = FALSE, fontface = "bold",
            color = "black") +
  scale_color_manual(values = cluster_colors, name = "Group",
                     labels = paste0("Group ", 1:6)) +
  scale_fill_manual(values = cluster_colors, name = "Group",
                    labels = paste0("Group ", 1:6)) +
  guides(color = guide_legend(title = "Group"),
         fill = guide_legend(title = "Group")) +
  scale_x_continuous(
    name = "Mean foregone harvest opportunity (%)",
    labels = function(x) paste0(round(x, 0), "%"),
    breaks = seq(0, 30, by = 5),
    limits = c(0, 30),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "Mean Protection Effectiveness (%)",
    labels = function(x) paste0(round(x, 0), "%"),
    breaks = seq(0, 55, by = 10),
    limits = c(0, 55),
    expand = c(0, 0)
  ) +
  theme_minimal(base_size = 18) +
  theme(
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(face = "bold", size = 18),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
    axis.title = element_text(face = "bold", size = 20),
    axis.text = element_text(face = "bold", size = 18),
    axis.line = element_line(color = "black", linewidth = 1),
    panel.grid.major = element_line(color = "grey80", linewidth = 0.5),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(file.path(OUTPUT_DIR, "cpue_vs_protection_combined_ALT.png"),
       p_combined_alt, width = 10, height = 10, dpi = 300, bg = "white")

cat("✓ Combined ALT plot saved\n")

# Version 2: Faceted by cluster
cat("Creating faceted plot by cluster...\n")

# Create a labeled version of the data for faceting
mean_data_labeled <- mean_data %>%
  mutate(cluster_label = factor(
    paste0("Group ", gsub("Cluster ", "", cluster)),
    levels = paste0("Group ", 1:6)
  ))

p_faceted <- ggplot(mean_data_labeled,
                    aes(x = mean_cpue_foregone, 
                        y = mean_protection,
                        color = cluster,
                        group = cluster)) +
  geom_line(linewidth = 1, alpha = 0.6) +
  geom_point(size = 4, alpha = 0.6, position = position_jitter(width = 0.3, height = 0, seed = 42),
             stroke = 2, shape = 21, aes(fill = cluster), color = "black") +
  geom_text(aes(label = offset_short), 
            vjust = -1.2, hjust = 0.5, size = 3, 
            show.legend = FALSE, fontface = "bold",
            color = "black") +
  scale_color_manual(values = cluster_colors, name = "Group") +
  scale_fill_manual(values = cluster_colors, name = "Group") +
  scale_x_continuous(
    name = "Mean CPUE Foregone (%)",
    labels = function(x) paste0(round(x, 0), "%"),
    breaks = scales::pretty_breaks(n = 8)
  ) +
  scale_y_continuous(
    name = "Mean Protection Effectiveness (%)",
    labels = function(x) paste0(round(x, 0), "%"),
    breaks = scales::pretty_breaks(n = 8)
  ) +
  facet_wrap(~ cluster_label, ncol = 2, scales = "free") +
  labs(
    title = "Mean CPUE Foregone vs. Front-End Closure Protection",
    subtitle = paste0("By Stock Group Across Timing Scenarios (", year_range, ")")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray95", color = NA),
    legend.position = "none",  # Remove legend since each panel shows one cluster
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray80"),
    plot.background = element_rect(fill = "white", color = NA)
  )

# [trimmed: variant of Fig 6] faceted CPUE-vs-protection plot not exported.
# ggsave(file.path(OUTPUT_DIR, "cpue_vs_protection_by_cluster.png"),
#        p_faceted, width = 12, height = 10, dpi = 300, bg = "white")

cat("✓ Faceted plot saved\n")

################################################################################
# STEP 6.5: CALCULATE CORRELATIONS
################################################################################

cat("\n--- Calculating correlations using mean data ---\n")

# Overall correlation (all mean data points)
overall_cor <- cor.test(mean_data$mean_cpue_foregone,
                        mean_data$mean_protection)

cat("\nOverall correlation (mean data):\n")
cat("  Pearson's r =", round(overall_cor$estimate, 3), "\n")
cat("  p-value =", format.pval(overall_cor$p.value, digits = 3), "\n")

# Correlations by cluster (using mean data)
cluster_correlations <- mean_data %>%
  group_by(cluster) %>%
  summarise(
    n_scenarios = n(),
    correlation = cor(mean_cpue_foregone, mean_protection, 
                      use = "complete.obs"),
    .groups = "drop"
  ) %>%
  arrange(desc(abs(correlation)))

cat("\nCorrelations by cluster (mean data):\n")
print(cluster_correlations, n = Inf)

################################################################################
# STEP 6.6: EXPORT ANALYSIS 6 DATA
################################################################################

cat("\n--- Exporting Analysis 6 data to CSV ---\n")

# Export mean data
write_csv(mean_data %>%
            select(cluster, offset_label, offset_short, n_years,
                   mean_cpue_foregone, sd_cpue_foregone, 
                   min_cpue_foregone, max_cpue_foregone,
                   mean_protection, sd_protection,
                   min_protection, max_protection),
          file.path(CSV_DIR, "analysis6_mean_cpue_vs_protection.csv"))

# Export correlation statistics
write_csv(cluster_correlations,
          file.path(CSV_DIR, "analysis6_correlations_by_cluster.csv"))

# Also export the full dataset for reference
write_csv(protection_cpue_combined %>%
            select(year, cluster, offset_label, 
                   foregone_proportion_pct, remaining_proportion_pct, 
                   closure_protection_pct,
                   cluster_total_production, cluster_closure_production),
          file.path(CSV_DIR, "analysis6_full_cpue_vs_protection_data.csv"))

cat("=== ANALYSIS 6 COMPLETE ===\n")
cat("Created: 2 scatter plots\n")
cat("  1. Combined plot - all groups in one panel\n")
cat("  2. Faceted plot - one panel per group\n")
cat("Exported: 3 CSV files\n")
cat("  - Mean CPUE foregone vs. protection data (with SD and range)\n")
cat("  - Correlations by cluster\n")
cat("  - Full dataset for reference\n\n")

cat("\n================================\n")
cat("ALL ANALYSES COMPLETE!\n")
cat("================================\n\n")
cat("PUBLICATION-READY MODIFICATIONS:\n")
cat("  ✓ Analysis 1: 4-panel boxplot figure (maps removed)\n")
cat("    - Figure size: 28 x 10 inches (wider and taller)\n")
cat("    - Boxplot outline thickness: 0.8 (slightly bolder for definition)\n")
cat("    - Single centered X-axis label: \"Production (%)\" at bottom\n")
cat("    - Boxplot appearance: Delicate boxes with slightly bolder outlines\n")
cat("    - Y-axis tick labels (10%, 30%, 50%): 24pt\n")
cat("    - Group labels (left side): 26pt\n")
cat("  ✓ Analyses 2-6: Unchanged from original\n")