################################################################################
# FIXED FRONT-END CLOSURE BOXPLOT WITH CORRECT WATERSHED ORDERING
################################################################################
# PURPOSE: Shows what % of EACH management unit's total annual production 
#          falls within the Q1 closure window with CORRECT watershed ordering
################################################################################

# Load libraries
library(dplyr)
library(ggplot2)
library(scales)

################################################################################
# CORRECTED WATERSHED ORDERING FUNCTIONS
################################################################################

#' Get standardized watershed ordering for all plots and analyses
#' 
#' @return Character vector of management unit names ordered from upstream to downstream
get_watershed_order <- function() {
  # Exact order from your document - upstream to downstream
  watershed_order <- c(
    "N. Fork Kusko",
    "E. Fork Kuskokwim River", 
    "S. Fork Kusko",
    "Upper Kusko Main",
    "Big River",
    "Takotna and Nixon Fork",
    "Tatlawiksuk",
    "Swift",
    "Stony", 
    "Holitna",
    "Hoholitna",
    "Middle Kusko Main",
    "George",
    "Oskakawlik", 
    "Holokuk",
    "Aniak",
    "Tuluksak",
    "Kisaralik",
    "Kwethluk",
    "Johnson",
    "Lower Kusko"
  )
  
  return(watershed_order)
}

#' Apply watershed ordering to a data frame
#' 
#' @param data Data frame containing management unit column
#' @param mgmt_col Name of the management unit column (default: "mgmt_river")
#' @param reverse_for_plots Logical, whether to reverse order for coord_flip() plots (default: FALSE)
#' @return Data frame with management unit column converted to ordered factor
apply_watershed_order <- function(data, mgmt_col = "mgmt_river", reverse_for_plots = FALSE) {
  
  # Get the standard ordering
  standard_order <- get_watershed_order()
  
  # Reverse if needed for coord_flip plots (so upstream appears at top)
  if (reverse_for_plots) {
    standard_order <- rev(standard_order)
  }
  
  # Filter to only units present in the data
  units_in_data <- unique(data[[mgmt_col]])
  final_order <- standard_order[standard_order %in% units_in_data]
  
  # Handle name variations/mismatches
  # Check for common variations and map them
  name_mapping <- c(
    "E. Fork Kuskokwim" = "E. Fork Kuskokwim River",
    "E. Fork Kuskokwim River" = "E. Fork Kuskokwim River"
  )
  
  # Apply any name mappings if needed
  for (old_name in names(name_mapping)) {
    if (old_name %in% units_in_data && !(name_mapping[old_name] %in% units_in_data)) {
      # Replace the old name with the new name in the data
      data[[mgmt_col]][data[[mgmt_col]] == old_name] <- name_mapping[old_name]
      units_in_data <- unique(data[[mgmt_col]])
      final_order <- standard_order[standard_order %in% units_in_data]
    }
  }
  
  # Add any units from data that aren't in our standard list
  missing_units <- setdiff(units_in_data, standard_order)
  if (length(missing_units) > 0) {
    warning("Found management units not in standard order: ", paste(missing_units, collapse = ", "))
    cat("Missing units that will be added at the end:\n")
    for(unit in missing_units) cat("  -", unit, "\n")
    final_order <- c(final_order, missing_units)
  }
  
  # Apply factor ordering
  data[[mgmt_col]] <- factor(data[[mgmt_col]], levels = final_order)
  
  return(data)
}

################################################################################
# MAIN ANALYSIS WITH CORRECTED ORDERING
################################################################################

# Output directory
OUTPUT_DIR <- "/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Figures/Front_End_Closure_Boxplots"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load ALL quartile data (Q1, Q2, Q3, Q4) to calculate totals
mgmt_data <- read.csv("/Users/benjaminmakhlouf/Research_repos/AYK_RunTiming/Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv")

# Remove Johnson river and clean quartile names
mgmt_data_clean <- mgmt_data %>%
  filter(mgmt_river != "Johnson") %>%
  mutate(quartile_clean = case_when(
    quartile == "Q1" ~ "Q1",
    quartile == "Q2" ~ "Q2", 
    quartile == "Q3" ~ "Q3",
    quartile == "Q4" ~ "Q4",
    TRUE ~ quartile
  ))

cat("Units found in data:\n")
unique_units <- sort(unique(mgmt_data_clean$mgmt_river))
for(i in 1:length(unique_units)) {
  cat(paste(i, ".", unique_units[i], "\n"))
}

# Calculate total annual production for each management unit in each year
annual_totals <- mgmt_data_clean %>%
  group_by(year, mgmt_river) %>%
  summarise(
    total_annual_production = sum(total_run_prop, na.rm = TRUE),
    .groups = "drop"
  )

# Get Q1 production values
q1_production <- mgmt_data_clean %>%
  filter(quartile_clean == "Q1") %>%
  select(year, mgmt_river, q1_production = total_run_prop)

# Join Q1 with annual totals and calculate closure percentage
boxplot_data <- q1_production %>%
  left_join(annual_totals, by = c("year", "mgmt_river")) %>%
  mutate(
    closure_percentage = (q1_production / total_annual_production) * 100
  ) %>%
  select(year, mgmt_river, q1_production, total_annual_production, closure_percentage)

################################################################################
# APPLY CORRECTED WATERSHED ORDERING
################################################################################

# Apply watershed ordering (reverse for coord_flip so upstream appears at top)
boxplot_data <- apply_watershed_order(boxplot_data, "mgmt_river", reverse_for_plots = TRUE)

# Print the final ordering to verify
cat("\nCorrected watershed order applied (will appear top to bottom in plot):\n")
final_levels <- levels(boxplot_data$mgmt_river)
for (i in 1:length(final_levels)) {
  cat(paste(i, ".", final_levels[i], "\n"))
}

################################################################################
# CREATE CORRECTED BOXPLOT
################################################################################

year_range <- paste0(min(boxplot_data$year), "-", max(boxplot_data$year))

boxplot <- ggplot(boxplot_data, aes(x = mgmt_river, y = closure_percentage)) +
  geom_boxplot(
    fill = "lightblue", 
    alpha = 0.7, 
    outlier.size = 2.5,
    outlier.shape = 16,
    linewidth = 0.6,
    color = "darkblue"
  ) +
  coord_flip() +
  scale_y_continuous(
    labels = function(x) paste0(round(x, 1), "%"),
    limits = c(0, max(boxplot_data$closure_percentage, na.rm = TRUE) * 1.05),
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  labs(
    title = "Front-End Closure Protection by Management Unit",
    subtitle = paste("% of EACH UNIT'S total annual production within Q1 closure window | Watershed order: upstream → downstream | Years:", year_range),
    x = "Management Unit (Watershed Position: Upstream → Downstream)",
    y = "% of Unit's Total Annual Production in Q1 Closure Window",
    caption = "Shows what % of each unit's total annual production (Q1+Q2+Q3+Q4) occurs within Q1 closure period\nOrdered by position in watershed from headwaters (top) to mouth (bottom)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, color = "grey20"),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    plot.caption = element_text(size = 10, hjust = 0.5, face = "italic", color = "gray50"),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
    axis.title = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 11),
    axis.text.x = element_text(size = 10),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(15, 15, 15, 15),
    panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.5)
  )

################################################################################
# SAVE AND SUMMARIZE
################################################################################

ggsave(file.path(OUTPUT_DIR, "front_end_closure_protection_boxplot_CORRECTED_ORDER.png"), 
       boxplot, width = 12, height = 10, dpi = 300, bg = "white")

print(boxplot)

# Summary statistics with corrected watershed ordering maintained
summary_stats <- boxplot_data %>%
  group_by(mgmt_river) %>%
  summarise(
    mean_closure_pct = round(mean(closure_percentage, na.rm = TRUE), 1),
    median_closure_pct = round(median(closure_percentage, na.rm = TRUE), 1),
    min_closure_pct = round(min(closure_percentage, na.rm = TRUE), 1),
    max_closure_pct = round(max(closure_percentage, na.rm = TRUE), 1),
    n_years = n(),
    .groups = "drop"
  ) %>%
  # Maintain the watershed ordering in the summary
  arrange(mgmt_river)

write.csv(summary_stats, file.path(OUTPUT_DIR, "closure_protection_summary_CORRECTED_ORDER.csv"), row.names = FALSE)

cat("\n=== SUMMARY STATISTICS (CORRECTED WATERSHED ORDERING) ===\n")
cat("Management units in correct upstream → downstream order:\n")
for(i in 1:nrow(summary_stats)) {
  unit_name <- summary_stats$mgmt_river[i]
  mean_pct <- summary_stats$mean_closure_pct[i]
  cat(paste(i, ".", unit_name, "- Mean:", mean_pct, "%\n"))
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("✓ Applied CORRECTED standardized watershed ordering (upstream → downstream)\n")
cat("✓ N. Fork Kusko (most upstream) appears at TOP of plot\n")
cat("✓ Lower Kusko (most downstream) appears at BOTTOM of plot\n")
cat(paste("✓ Plot saved to:", file.path(OUTPUT_DIR, "front_end_closure_protection_boxplot_CORRECTED_ORDER.png")))

################################################################################
# DIAGNOSTIC: CHECK FOR NAME MISMATCHES
################################################################################

cat("\n=== DIAGNOSTIC INFORMATION ===\n")
cat("Standard watershed order:\n")
standard_order <- get_watershed_order()
for(i in 1:length(standard_order)) {
  cat(paste(i, ".", standard_order[i], "\n"))
}

cat("\nUnits actually found in data:\n")
data_units <- sort(unique(mgmt_data_clean$mgmt_river))
for(i in 1:length(data_units)) {
  in_standard <- data_units[i] %in% standard_order
  cat(paste(i, ".", data_units[i], ifelse(in_standard, "(✓ in standard)", "(⚠ NOT in standard)"), "\n"))
}

# Check for potential name mismatches
cat("\nPotential name mismatches to investigate:\n")
mismatched <- setdiff(data_units, standard_order)
if(length(mismatched) > 0) {
  for(unit in mismatched) {
    cat(paste("- '", unit, "' not found in standard order\n", sep=""))
  }
} else {
  cat("All units match the standard order perfectly!\n")
}