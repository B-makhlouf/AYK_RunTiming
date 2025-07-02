################################################################################
#                    AVERAGE MANAGEMENT UNIT PRODUCTION MAPS                  #
################################################################################
# PURPOSE: Create maps and boxplots with standardized watershed ordering
# NOTE: Now uses watershed ordering (upstream to downstream) for all plots
################################################################################

################################################################################
#                              SETUP & LIBRARIES                              #
################################################################################

# Load required libraries
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
  library(gridExtra)
  library(readr)
  library(grid)
  library(png)
})

# Source spatial utils for watershed ordering functions
source("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/code/utils/spatial_utils.R")

# Set file paths
DATA_PATH <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv"
SPATIAL_PATH <- "/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC_joined.shp"
BASIN_PATH <- "/Users/benjaminmakhlouf/Desktop/Research/isoscapes_new/Kusko/Kusko_basin.shp"
OUTPUT_DIR <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Figures/Average_Management_Production_2017_2022_WATERSHED_ORDERED"

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== AVERAGE MANAGEMENT UNIT PRODUCTION ANALYSIS (WATERSHED ORDERED) ===\n")
cat("Creating maps and boxplots with standardized watershed ordering\n")
cat("Output directory:", OUTPUT_DIR, "\n")

################################################################################
#                              DATA PREPARATION                               #
################################################################################

cat("\n--- LOADING AND PREPARING DATA ---\n")

# Load management unit production data
mgmt_data <- read_csv(DATA_PATH) %>%
  filter(mgmt_river != "Johnson") %>%  # Remove Johnson as in other analyses
  mutate(
    quartile_clean = case_when(
      quartile == "Q1" ~ "Q1",
      quartile == "Q2" ~ "Q2", 
      quartile == "Q3" ~ "Q3",
      quartile == "Q4" ~ "Q4",
      TRUE ~ quartile
    )
  )

# Check what years are available
available_years <- sort(unique(mgmt_data$year))
cat("Years available in data:", paste(available_years, collapse = ", "), "\n")

# Check if 2022 data is present
if (2022 %in% available_years) {
  cat("✓ 2022 data found in dataset\n")
  year_range <- "2017-2022"
  years_included <- available_years
} else {
  cat("⚠ 2022 data not found - using available years only\n")
  year_range <- paste0(min(available_years), "-", max(available_years))
  years_included <- available_years
}

cat("Loaded data for", length(unique(mgmt_data$mgmt_river)), "management units and", 
    length(years_included), "years\n")

# Load spatial data
edges <- st_read(SPATIAL_PATH, quiet = TRUE)
basin <- st_read(BASIN_PATH, quiet = TRUE)

cat("Spatial data loaded successfully\n")

################################################################################
#                     CALCULATE AVERAGE PRODUCTION BY QUARTILE                #
################################################################################

cat("\n--- CALCULATING AVERAGES BY QUARTILE ---\n")

# Calculate average within-quartile production for each management unit and quartile
avg_production_by_quartile <- mgmt_data %>%
  group_by(mgmt_river, quartile_clean) %>%
  summarise(
    avg_within_quartile_prop = mean(within_quartile_prop, na.rm = TRUE),
    n_years = n(),
    .groups = "drop"
  ) %>%
  # Convert to the same scale as original maps (0-1 proportion)
  mutate(production_proportion = avg_within_quartile_prop)

# Apply watershed ordering
avg_production_by_quartile <- apply_watershed_order(avg_production_by_quartile, "mgmt_river", reverse_for_plots = FALSE)

cat("Calculated averages for", nrow(avg_production_by_quartile), "mgmt_unit × quartile combinations\n")

################################################################################
#              CREATE MAPS IDENTICAL TO EXISTING MANAGEMENT MAPS              #
################################################################################

cat("\n--- CREATING QUARTILE MAPS WITH WATERSHED ORDERING ---\n")

#' Create management map identical to original DOY_Quartile/Management style
create_average_mgmt_map <- function(quartile_data, quartile_label, output_filepath, year_range) {
  
  # Create the PNG with dimensions for single map (no bar plot)
  png(file = output_filepath, width = 10, height = 8, units = "in", res = 300, bg = "white")
  
  # Join with spatial data and set up linewidths like DFA maps
  edges_with_data <- edges %>%
    left_join(quartile_data, by = "mgmt_river") %>%
    filter(!is.na(mgmt_river) & mgmt_river != "") %>%
    mutate(
      # Set production proportion (handle NA values)
      production_proportion = ifelse(is.na(production_proportion), 0, production_proportion),
      # Set stream order and linewidth exactly like DFA maps
      stream_order = ifelse(is.na(Str_Order), 3, Str_Order),
      line_width = pmax(0.3, pmin(3.0, 0.3 + (stream_order - min(stream_order, na.rm = TRUE)) * 
                                    (3.0 - 0.3) / (max(stream_order, na.rm = TRUE) - min(stream_order, na.rm = TRUE))))
    )
  
  # Ensure consistent CRS
  if (st_crs(basin) != st_crs(edges_with_data)) {
    basin <- st_transform(basin, st_crs(edges_with_data))
  }
  
  # Main map plot - single map without bar chart
  main_plot <- ggplot() +
    geom_sf(data = basin, fill = "gray95", color = "gray70", 
            linewidth = 0.5, alpha = 0.3) +
    geom_sf(data = edges_with_data, 
            aes(color = production_proportion, linewidth = stream_order), 
            alpha = 0.8) +
    scale_color_gradientn(
      colors = brewer.pal(9, "YlOrRd"),
      name = "Production\nProportion",
      na.value = "grey60",
      labels = scales::percent_format(accuracy = 1),
      guide = guide_colorbar(
        barwidth = 1, barheight = 15,
        frame.colour = "grey40", ticks.colour = "grey40",
        show.limits = TRUE
      )
    ) +
    scale_linewidth_continuous(
      range = c(0.3, 3.0), 
      name = "Stream\nOrder"
    ) +
    coord_sf(datum = NA) +
    labs(
      title = paste0("Average ", quartile_label, ": Management Rivers - Kusko Watershed"),
      subtitle = paste("Average across all years (", year_range, ")", sep = "")
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5, color = "grey30"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey50"),
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold", color = "grey30"),
      legend.text = element_text(color = "grey30"),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(10, 10, 10, 10, "mm")
    )
  
  print(main_plot)
  
  dev.off()
  
  cat("Created average map:", basename(output_filepath), "\n")
}

# Create maps for each quartile
quartile_colors <- c("Q1" = "#E31A1C", "Q2" = "#1F78B4", 
                     "Q3" = "#33A02C", "Q4" = "#FF7F00")

for (q in c("Q1", "Q2", "Q3", "Q4")) {
  cat("Creating average map for", q, "\n")
  
  # Get average data for this quartile
  quartile_data <- avg_production_by_quartile %>%
    filter(quartile_clean == q) %>%
    filter(!is.na(production_proportion))
  
  # Create map with exact original styling
  map_filename <- paste0("Average_", q, "_Management_", year_range, "_WATERSHED_ORDERED.png")
  map_filepath <- file.path(OUTPUT_DIR, map_filename)
  
  create_average_mgmt_map(quartile_data, q, map_filepath, year_range)
}

cat("✓ All average quartile maps created with watershed ordering\n")

################################################################################
#                       CREATE CLEAN QUARTILE BOXPLOTS                        #
################################################################################

cat("\n--- CREATING CLEAN QUARTILE BOXPLOTS WITH WATERSHED ORDERING ---\n")

# Prepare data for boxplots - show variability in EACH QUARTILE across years
quartile_boxplot_data <- mgmt_data %>%
  select(year, mgmt_river, quartile_clean, within_quartile_prop) %>%
  # Convert to percentage for easier interpretation
  mutate(within_quartile_pct = within_quartile_prop * 100)

# Apply watershed ordering (reverse for coord_flip so upstream appears at top)
quartile_boxplot_data <- apply_watershed_order(quartile_boxplot_data, "mgmt_river", reverse_for_plots = TRUE)

cat("\nWatershed ordering applied (top to bottom in plots):\n")
final_levels <- levels(quartile_boxplot_data$mgmt_river)
for (i in 1:length(final_levels)) {
  cat(paste(i, ".", final_levels[i], "\n"))
}

# Create individual CLEAN boxplot for each quartile
for (q in c("Q1", "Q2", "Q3", "Q4")) {
  cat("Creating clean boxplot for", q, "\n")
  
  # Filter data for this quartile
  q_data <- quartile_boxplot_data %>%
    filter(quartile_clean == q)
  
  # Create CLEAN boxplot for this quartile - NO jittered points
  q_boxplot <- ggplot(q_data, aes(x = mgmt_river, y = within_quartile_pct)) +
    geom_boxplot(
      fill = "grey90", 
      alpha = 0.7, 
      outlier.alpha = 0.8,          # Make outliers clearly visible
      outlier.size = 2.5,            # Larger outliers for visibility
      outlier.shape = 16,            # Solid circles for outliers
      linewidth = 0.6                # Slightly thicker boxplot lines
    ) +
    # Y-axis formatting
    scale_y_continuous(
      labels = function(x) paste0(round(x, 1), "%"),
      limits = c(0, max(quartile_boxplot_data$within_quartile_pct, na.rm = TRUE) * 1.05),
      expand = expansion(mult = c(0.02, 0.05))  # Small expansion for clean look
    ) +
    # Flip coordinates for better readability
    coord_flip() +
    # Labels and theme - UPDATE THE SUBTITLE TO INCLUDE ACTUAL YEARS
    labs(
      title = paste("Management Unit Production Variability:", q),
      subtitle = paste("Within-", q, " production for each management unit across years (", year_range, ") | Watershed order: upstream → downstream", sep = ""),
      x = "Management Unit (Watershed Position)",
      y = paste(q, "Production (%)"),
      caption = "Boxes show median, quartiles, and whiskers; dots show outliers | Ordered by watershed position (upstream at top)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "grey50"),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
      plot.caption = element_text(size = 10, hjust = 0.5, face = "italic", color = "gray50"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),      # Remove horizontal grid lines
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),  # Keep faint vertical lines
      axis.title = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 11),
      axis.text.x = element_text(size = 10),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(15, 15, 15, 15),
      panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)  # Add subtle border
    )
  
  # Save individual quartile boxplot
  ggsave(file.path(OUTPUT_DIR, paste0("management_unit_", q, "_boxplot_WATERSHED_ORDERED.png")), 
         q_boxplot, width = 12, height = 10, dpi = 300, bg = "white")
}

cat("✓ Clean quartile boxplots created with watershed ordering\n")

################################################################################
#              CREATE COMBINED FACETED BOXPLOT (ALL QUARTILES)                #
################################################################################

cat("\n--- CREATING COMBINED FACETED BOXPLOT WITH WATERSHED ORDERING ---\n")

# Create combined faceted plot showing all quartiles
combined_boxplot <- ggplot(quartile_boxplot_data, aes(x = mgmt_river, y = within_quartile_pct)) +
  geom_boxplot(
    fill = "grey90", 
    alpha = 0.7, 
    outlier.alpha = 0.8,
    outlier.size = 2,
    outlier.shape = 16,
    linewidth = 0.5
  ) +
  facet_wrap(~quartile_clean, scales = "free_y", ncol = 2) +
  scale_y_continuous(
    labels = function(x) paste0(round(x, 1), "%"),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  coord_flip() +
  labs(
    title = "Management Unit Production Variability Across All Quartiles",
    subtitle = paste("Within-quartile production for each management unit (", year_range, ") | Watershed order: upstream → downstream", sep = ""),
    x = "Management Unit (Watershed Position)",
    y = "Production (%)",
    caption = "Boxes show median, quartiles, and whiskers; dots show outliers | Ordered by watershed position (upstream at top)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "grey30"),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic", color = "gray50"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    strip.text = element_text(face = "bold", size = 12),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(15, 15, 15, 15)
  )

ggsave(file.path(OUTPUT_DIR, "management_unit_all_quartiles_boxplot_WATERSHED_ORDERED.png"), 
       combined_boxplot, width = 14, height = 12, dpi = 300, bg = "white")

################################################################################
#                          SAVE SUMMARY DATA                                  #
################################################################################

cat("\n--- SAVING SUMMARY DATA WITH WATERSHED ORDERING ---\n")

# Save average production data with updated filename
write_csv(avg_production_by_quartile, 
          file.path(OUTPUT_DIR, paste0("average_production_by_quartile_", year_range, "_WATERSHED_ORDERED.csv")))

# Calculate summary statistics for quartiles with watershed ordering
quartile_summary_stats <- quartile_boxplot_data %>%
  group_by(mgmt_river, quartile_clean) %>%
  summarise(
    mean_quartile_production = mean(within_quartile_pct, na.rm = TRUE),
    median_quartile_production = median(within_quartile_pct, na.rm = TRUE),
    sd_quartile_production = sd(within_quartile_pct, na.rm = TRUE),
    min_quartile_production = min(within_quartile_pct, na.rm = TRUE),
    max_quartile_production = max(within_quartile_pct, na.rm = TRUE),
    cv_quartile_production = sd_quartile_production / mean_quartile_production,
    n_years = n(),
    .groups = "drop"
  ) %>%
  arrange(mgmt_river, quartile_clean)  # This maintains watershed ordering

write_csv(quartile_summary_stats, 
          file.path(OUTPUT_DIR, "management_unit_quartile_summary_statistics_WATERSHED_ORDERED.csv"))

# Calculate detailed statistics with proper scaling and watershed ordering
detailed_stats <- quartile_summary_stats %>%
  rename(
    mean_pct = mean_quartile_production,
    median_pct = median_quartile_production,
    sd_pct = sd_quartile_production,
    min_pct = min_quartile_production,
    max_pct = max_quartile_production,
    cv_pct = cv_quartile_production
  ) %>%
  arrange(quartile_clean, mgmt_river)  # Maintains watershed ordering within quartiles

# Save detailed statistics
write_csv(detailed_stats, 
          file.path(OUTPUT_DIR, "management_unit_detailed_statistics_WATERSHED_ORDERED.csv"))

cat("✓ Summary data saved with watershed ordering maintained\n")

################################################################################
#                               FINAL SUMMARY                                 #
################################################################################

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("✓ Created 4 average quartile maps with watershed ordering\n")
cat("✓ Added stream order-based line widths (DFA style)\n")
cat("✓ Created 4 individual quartile boxplot figures with watershed ordering\n")
cat("✓ Created 1 combined faceted boxplot with watershed ordering\n")
cat("✓ Generated summary statistics maintaining watershed order\n")
cat("✓ Year range:", year_range, "\n")
cat("✓ Years included:", paste(years_included, collapse = ", "), "\n")
cat("✓ WATERSHED ORDERING: Upstream (top/left) → Downstream (bottom/right)\n")
cat("\nAll outputs saved to:", OUTPUT_DIR, "\n")

cat("\nFiles created:\n")
files_created <- c(
  paste0("Average_Q1_Management_", year_range, "_WATERSHED_ORDERED.png"),
  paste0("Average_Q2_Management_", year_range, "_WATERSHED_ORDERED.png"), 
  paste0("Average_Q3_Management_", year_range, "_WATERSHED_ORDERED.png"),
  paste0("Average_Q4_Management_", year_range, "_WATERSHED_ORDERED.png"),
  "management_unit_Q1_boxplot_WATERSHED_ORDERED.png",
  "management_unit_Q2_boxplot_WATERSHED_ORDERED.png",
  "management_unit_Q3_boxplot_WATERSHED_ORDERED.png", 
  "management_unit_Q4_boxplot_WATERSHED_ORDERED.png",
  "management_unit_all_quartiles_boxplot_WATERSHED_ORDERED.png",
  paste0("average_production_by_quartile_", year_range, "_WATERSHED_ORDERED.csv"),
  "management_unit_quartile_summary_statistics_WATERSHED_ORDERED.csv",
  "management_unit_detailed_statistics_WATERSHED_ORDERED.csv"
)

for (file in files_created) {
  cat("  -", file, "\n")
}

cat("\n=== WATERSHED ORDERING APPLIED TO ALL VISUALIZATIONS ===\n")
cat("✓ Maps: Spatial representation (not affected by ordering)\n")
cat("✓ Boxplots: Upstream management units at TOP, downstream at BOTTOM\n")
cat("✓ Bar plots: Would follow same upstream → downstream pattern\n")
cat("✓ CSV files: Maintain watershed ordering for consistency\n")
cat("✓ All analysis maintains biological/spatial relevance\n")