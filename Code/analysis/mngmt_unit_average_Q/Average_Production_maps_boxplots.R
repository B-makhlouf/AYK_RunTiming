################################################################################
#                    AVERAGE MANAGEMENT UNIT PRODUCTION MAPS                  #
################################################################################
# PURPOSE: Create maps that look EXACTLY like the existing DOY_Quartile/Management
#          maps, but showing average production across all years (2017-2021)
#
# OUTPUTS: 
#   - 4 maps (Q1-Q4) with identical styling to existing management maps
#   - 1 boxplot showing production variability across years
#   - Summary statistics CSV
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

# Set file paths
DATA_PATH <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv"
SPATIAL_PATH <- "/Users/benjaminmakhlouf/Spatial Data/KuskoUSGS_HUC_joined.shp"
BASIN_PATH <- "/Users/benjaminmakhlouf/Desktop/Research/isoscapes_new/Kusko/Kusko_basin.shp"
OUTPUT_DIR <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Figures/Average_Management_Production"

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("=== AVERAGE MANAGEMENT UNIT PRODUCTION ANALYSIS ===\n")
cat("Creating maps identical to DOY_Quartile/Management style\n")
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

cat("Loaded data for", length(unique(mgmt_data$year)), "years and", 
    length(unique(mgmt_data$mgmt_river)), "management units\n")

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

cat("Calculated averages for", nrow(avg_production_by_quartile), "mgmt_unit × quartile combinations\n")

################################################################################
#              CREATE MAPS IDENTICAL TO EXISTING MANAGEMENT MAPS              #
################################################################################

cat("\n--- CREATING QUARTILE MAPS WITH ORIGINAL STYLING ---\n")

#' Create management map identical to original DOY_Quartile/Management style
create_average_mgmt_map <- function(quartile_data, quartile_label, output_filepath) {
  
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
      subtitle = "Average across all years (2017-2021)"
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
  map_filename <- paste0("Average_", q, "_Management_2017_2021.png")
  map_filepath <- file.path(OUTPUT_DIR, map_filename)
  
  create_average_mgmt_map(quartile_data, q, map_filepath)
}

cat("✓ All average quartile maps created with original styling\n")

################################################################################
#                       CREATE CLEAN QUARTILE BOXPLOTS                        #
################################################################################

cat("\n--- CREATING CLEAN QUARTILE BOXPLOTS ---\n")

# Prepare data for boxplots - show variability in EACH QUARTILE across years
quartile_boxplot_data <- mgmt_data %>%
  select(year, mgmt_river, quartile_clean, within_quartile_prop) %>%
  # Convert to percentage for easier interpretation
  mutate(within_quartile_pct = within_quartile_prop * 100)

# First, let's see what management units we actually have in the data
cat("Management units found in data:\n")
unique_mgmt_units <- sort(unique(quartile_boxplot_data$mgmt_river))
for (i in seq_along(unique_mgmt_units)) {
  cat(paste(i, ".", unique_mgmt_units[i], "\n"))
}

# Define the EXACT order from your image (REVERSED for coord_flip to show N. Fork Kusko at top)
desired_order <- c(
  "N. Fork Kusko",
  "E. Fork Kuskokwim", 
  "S. Fork Kusko", 
  "Takotna and Nixon Fork",
  "Big River",
  "Upper Kusko Main",
  "Tatlawiksuk",
  "Kwethluk",
  "Stony",
  "Swift",
  "Holitna and Hoholitna", 
  "George",
  "Oskakawlik",
  "Middle Kusko Main",
  "Holokuk",
  "Aniak",
  "Tuluksak",
  "Kisaralik",
  "Hoholitna",
  "Johnson",
  "Lower Kusko"  # Changed from "Lower Kusko Main"
)

# Reverse the order for coord_flip
mgmt_order <- rev(desired_order)

# Only include units that actually exist in the data, in the specified order
existing_units <- intersect(mgmt_order, unique_mgmt_units)

# Add any units from data that weren't in our specified list
additional_units <- setdiff(unique_mgmt_units, desired_order)
final_order <- c(existing_units, additional_units)

quartile_boxplot_data <- quartile_boxplot_data %>%
  mutate(mgmt_river = factor(mgmt_river, levels = final_order))

cat("\nFinal ordering will be (top to bottom after coord_flip):\n")
final_levels <- levels(quartile_boxplot_data$mgmt_river)
for (i in length(final_levels):1) {  # Print in reverse to show actual display order
  cat(paste(length(final_levels) - i + 1, ".", final_levels[i], "\n"))
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
    # Labels and theme
    labs(
      title = paste("Management Unit Production Variability:", q),
      subtitle = paste("Within-", q, " production for each management unit across years (2017-2021)", sep = ""),
      x = "Management Unit",
      y = paste(q, "Production (%)"),
      caption = "Boxes show median, quartiles, and whiskers; dots show outliers"
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
  ggsave(file.path(OUTPUT_DIR, paste0("management_unit_", q, "_boxplot.png")), 
         q_boxplot, width = 12, height = 10, dpi = 300, bg = "white")
}

cat("✓ Clean quartile boxplots created and saved\n")

################################################################################
#                          SAVE SUMMARY DATA                                  #
################################################################################

cat("\n--- SAVING SUMMARY DATA ---\n")

# Save average production data
write_csv(avg_production_by_quartile, 
          file.path(OUTPUT_DIR, "average_production_by_quartile_2017_2021.csv"))

# Calculate summary statistics for quartiles
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
  arrange(mgmt_river, quartile_clean)

write_csv(quartile_summary_stats, 
          file.path(OUTPUT_DIR, "management_unit_quartile_summary_statistics.csv"))

cat("✓ Summary data saved\n")

################################################################################
#                               FINAL SUMMARY                                 #
################################################################################

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("✓ Created 4 average quartile maps (no bar plots)\n")
cat("✓ Added stream order-based line widths (DFA style)\n")
cat("✓ Created 4 individual quartile boxplot figures\n")
cat("✓ Generated summary statistics\n")
cat("\nAll outputs saved to:", OUTPUT_DIR, "\n")

cat("\nFiles created:\n")
files_created <- c(
  "Average_Q1_Management_2017_2021.png",
  "Average_Q2_Management_2017_2021.png", 
  "Average_Q3_Management_2017_2021.png",
  "Average_Q4_Management_2017_2021.png",
  "management_unit_Q1_boxplot.png",
  "management_unit_Q2_boxplot.png",
  "management_unit_Q3_boxplot.png", 
  "management_unit_Q4_boxplot.png",
  "average_production_by_quartile_2017_2021.csv",
  "management_unit_quartile_summary_statistics.csv"
)

for (file in files_created) {
  cat("  -", file, "\n")
}

cat("\n=== Maps: No bar plots, DFA-style line widths ===\n")
cat("=== Boxplots: Clean style with outliers only ===\n")

# Save detailed statistics
write_csv(detailed_stats, 
          file.path(OUTPUT_DIR, "management_unit_detailed_statistics_clean.csv"))

cat("✓ Detailed summary statistics saved\n")

# Print top units by quartile
cat("\nTop 3 management units by mean production for each quartile:\n")
for (q in c("Q1", "Q2", "Q3", "Q4")) {
  cat(paste("\n", q, ":\n"))
  top_units <- detailed_stats %>%
    filter(quartile_clean == q) %>%
    slice_head(n = 3) %>%
    select(mgmt_river, mean_pct, median_pct, cv_pct)
  print(top_units)
}

################################################################################
#                               FINAL SUMMARY                                 #
################################################################################

cat("\n=== CLEAN BOXPLOT ANALYSIS COMPLETE ===\n")
cat("✓ Created 4 individual clean quartile boxplots\n")
cat("✓ Created 1 combined quartile boxplot (faceted)\n")
cat("✓ Removed jittered points for cleaner visualization\n")
cat("✓ Enhanced outlier visibility\n")
cat("✓ Generated detailed summary statistics\n")
cat("\nAll outputs saved to:", OUTPUT_DIR, "\n")

cat("\nClean boxplot files created:\n")
files_created <- c(
  "management_unit_Q1_boxplot_clean.png",
  "management_unit_Q2_boxplot_clean.png",
  "management_unit_Q3_boxplot_clean.png", 
  "management_unit_Q4_boxplot_clean.png",
  "management_unit_all_quartiles_boxplot_clean.png",
  "management_unit_detailed_statistics_clean.csv"
)

for (file in files_created) {
  cat("  -", file, "\n")
}

cat("\n=== Clean Boxplots: Just boxes, whiskers, and outliers ===\n")
cat("=== No jittered points for cleaner, more professional look ===\n")