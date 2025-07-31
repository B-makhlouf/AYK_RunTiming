################################################################################
# FRONT-END CLOSURE BOXPLOT ANALYSIS - STANDALONE SCRIPT
################################################################################
# PURPOSE: Creates contemporary boxplot showing what % of each management unit's 
#          production falls within Q1 closure window with matching colors and ordering
################################################################################

# Load required libraries
library(sf)
library(dplyr)
library(here)
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(scales)
library(readr)

################################################################################
# FRONT-END CLOSURE BOXPLOT ANALYSIS
################################################################################

#' Run Front-end Closure Boxplot Analysis with CONTEMPORARY design and matching colors
run_closure_boxplot_analysis <- function(years, watershed = "Kusko") {
  
  message("=== Creating Contemporary Closure Boxplot with Matching Colors ===")
  
  # Output directory
  OUTPUT_DIR <- here("Figures/Front_End_Closure_Boxplots")
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Load management river analysis data
  DATA_PATH <- here("Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv")
  
  if (!file.exists(DATA_PATH)) {
    warning("Management river analysis data not found. Run DOY analysis first.")
    warning("Expected file: ", DATA_PATH)
    return(NULL)
  }
  
  mgmt_data <- read.csv(DATA_PATH) %>%
    filter(year %in% years) %>%  # Keep ALL management units including Johnson
    mutate(quartile_clean = case_when(
      quartile == "Q1" ~ "Q1",
      quartile == "Q2" ~ "Q2", 
      quartile == "Q3" ~ "Q3",
      quartile == "Q4" ~ "Q4",
      TRUE ~ quartile
    ))
  
  # Check what management units we have (should include Johnson)
  all_mgmt_units <- sort(unique(mgmt_data$mgmt_river))
  message(paste("ALL management units in data (including Johnson):", length(all_mgmt_units)))
  for (i in 1:length(all_mgmt_units)) {
    message(paste("  ", i, ". '", all_mgmt_units[i], "'", sep = ""))
  }
  
  # Check Johnson's data specifically
  johnson_data <- mgmt_data %>% filter(mgmt_river == "Johnson")
  if (nrow(johnson_data) > 0) {
    message(paste("\nJohnson data summary:"))
    message(paste("  Rows:", nrow(johnson_data)))
    message(paste("  Total_run_prop values:", paste(unique(johnson_data$total_run_prop), collapse = ", ")))
    message(paste("  Any NAs:", any(is.na(johnson_data$total_run_prop))))
    message(paste("  Any zeros:", any(johnson_data$total_run_prop == 0, na.rm = TRUE)))
  } else {
    message("\nWARNING: No Johnson data found!")
  }
  
  # Calculate total annual production (including Johnson with proper NA/0 handling)
  annual_totals <- mgmt_data %>%
    group_by(year, mgmt_river) %>%
    summarise(
      total_annual_production = sum(total_run_prop, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Keep Johnson even if total is 0 or NA
    mutate(total_annual_production = ifelse(is.na(total_annual_production), 0, total_annual_production))
  
  # Get Q1 production values (including Johnson)
  q1_production <- mgmt_data %>%
    filter(quartile_clean == "Q1") %>%
    select(year, mgmt_river, q1_production = total_run_prop) %>%
    # Handle Johnson's NA/0 values
    mutate(q1_production = ifelse(is.na(q1_production), 0, q1_production))
  
  # Create complete list of all management units and years to ensure Johnson is included
  all_combinations <- expand.grid(
    year = unique(mgmt_data$year),
    mgmt_river = unique(mgmt_data$mgmt_river),
    stringsAsFactors = FALSE
  )
  
  # Join with proper handling of missing values
  boxplot_data <- all_combinations %>%
    left_join(q1_production, by = c("year", "mgmt_river")) %>%
    left_join(annual_totals, by = c("year", "mgmt_river")) %>%
    mutate(
      # Handle cases where Johnson or other units have 0 or NA values
      q1_production = ifelse(is.na(q1_production), 0, q1_production),
      total_annual_production = ifelse(is.na(total_annual_production), 0, total_annual_production),
      # Calculate closure percentage, handling division by zero
      closure_percentage = case_when(
        total_annual_production == 0 ~ 0,  # If no annual production, closure % is 0
        TRUE ~ (q1_production / total_annual_production) * 100
      )
    ) %>%
    select(year, mgmt_river, q1_production, total_annual_production, closure_percentage)
  
  # Check what we have after the join operations (should be 21 including Johnson)
  final_mgmt_units <- sort(unique(boxplot_data$mgmt_river))
  message(paste("\nFinal management units in boxplot data:", length(final_mgmt_units)))
  for (i in 1:length(final_mgmt_units)) {
    unit <- final_mgmt_units[i]
    unit_data <- boxplot_data %>% filter(mgmt_river == unit)
    mean_closure <- round(mean(unit_data$closure_percentage, na.rm = TRUE), 1)
    message(paste("  ", i, ". '", unit, "' (mean closure: ", mean_closure, "%)", sep = ""))
  }
  
  # Special check for Johnson
  johnson_boxplot_data <- boxplot_data %>% filter(mgmt_river == "Johnson")
  if (nrow(johnson_boxplot_data) > 0) {
    message(paste("\nJohnson in final boxplot data:"))
    message(paste("  Rows:", nrow(johnson_boxplot_data)))
    message(paste("  Closure percentages:", paste(round(johnson_boxplot_data$closure_percentage, 1), collapse = ", ")))
    message(paste("  Mean closure %:", round(mean(johnson_boxplot_data$closure_percentage, na.rm = TRUE), 1)))
  } else {
    message("\nERROR: Johnson missing from final boxplot data!")
  }
  
  # GET EXACT SAME ORDERING AS CUMULATIVE PLOTS
  # Use the EXACT same approach as your cumulative distribution script
  
  # First, let's see what management units are actually in the boxplot data
  actual_mgmt_rivers <- sort(unique(boxplot_data$mgmt_river))
  message(paste("\nFINAL management units for plotting:", length(actual_mgmt_rivers)))
  for (i in 1:length(actual_mgmt_rivers)) {
    message(paste("  ", i, ". '", actual_mgmt_rivers[i], "'", sep = ""))
  }
  
  # Your cumulative script should have 20 units (21 minus Johnson)
  # Let's make sure we match exactly what your cumulative script produces
  expected_count <- 20  # 21 total minus Johnson
  if (length(actual_mgmt_rivers) != expected_count) {
    message(paste("\nWARNING: Expected", expected_count, "units but found", length(actual_mgmt_rivers)))
  }
  
  # Check for name variations (especially East Fork)
  east_fork_variations <- c(
    "E. Fork Kuskokwim River",
    "E. Fork Kuskokwim", 
    "East Fork Kuskokwim",
    "East Fork Kuskokwim River"
  )
  
  found_east_fork <- intersect(east_fork_variations, actual_mgmt_rivers)
  if (length(found_east_fork) > 0) {
    message(paste("Found East Fork as: '", found_east_fork[1], "'", sep = ""))
    
    # If the actual name differs from desired, update desired_order
    if (found_east_fork[1] != "E. Fork Kuskokwim") {
      desired_order[desired_order == "E. Fork Kuskokwim"] <- found_east_fork[1]
      message(paste("Updated desired order to use actual name: '", found_east_fork[1], "'", sep = ""))
    }
  } else {
    message("WARNING: East Fork not found with expected names!")
    east_candidates <- actual_mgmt_rivers[grepl("East|E\\.|Fork|Kusko", actual_mgmt_rivers, ignore.case = TRUE)]
    for (candidate in east_candidates) {
      message(paste("  Candidate: '", candidate, "'", sep = ""))
    }
  }
  
  # EXACT DESIRED ORDER - matching your specification and actual data names
  desired_order <- c(
    "N. Fork Kusko",
    "E. Fork Kuskokwim River",  # Updated to match actual data name
    "S. Fork Kusko", 
    "Upper Kusko Main",
    "Big River",
    "Takotna and Nixon Fork",
    "Tatlawiksuk",
    "Swift",
    "Stony", 
    "Holitna and Hoholitna",    # This appears to be the actual combined name
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
  
  # Try to match management units exactly with desired order
  matched_units <- intersect(desired_order, actual_mgmt_rivers)
  unmatched_units <- setdiff(actual_mgmt_rivers, desired_order)
  
  message(paste("\nMATCHED units from desired order:", length(matched_units)))
  for (i in 1:length(matched_units)) {
    message(paste("  ", i, ". MATCHED: '", matched_units[i], "'", sep = ""))
  }
  
  message(paste("\nUNMATCHED units (will be added at end):", length(unmatched_units)))
  for (unit in unmatched_units) {
    message(paste("  UNMATCHED: '", unit, "'", sep = ""))
  }
  
  # Create final ordering: matched units in desired order + unmatched units at end
  # Use the desired_order sequence for matched units to preserve the exact order
  final_order <- c()
  for (unit in desired_order) {
    if (unit %in% actual_mgmt_rivers) {
      final_order <- c(final_order, unit)
    }
  }
  # Add any unmatched units at the end
  final_order <- c(final_order, unmatched_units)
  
  message("\nFINAL ORDER FOR BOXPLOT (upstream → downstream):")
  for (i in 1:length(final_order)) {
    message(paste("  ", i, ". '", final_order[i], "'", sep = ""))
  }
  
  # Double-check we have all the original units
  if (length(final_order) != length(actual_mgmt_rivers)) {
    message("\nWARNING: Mismatch in number of management units!")
    message(paste("Original:", length(actual_mgmt_rivers), "Final:", length(final_order)))
    missing_units <- setdiff(actual_mgmt_rivers, final_order)
    if (length(missing_units) > 0) {
      message("Missing units:")
      for (unit in missing_units) {
        message(paste("  MISSING: '", unit, "'", sep = ""))
      }
      # Add missing units to final order
      final_order <- c(final_order, missing_units)
    }
  }
  
  # Apply factor ordering (reverse for coord_flip so upstream appears at top)
  boxplot_data$mgmt_river <- factor(boxplot_data$mgmt_river, levels = rev(final_order))
  
  # CREATE EXACT SAME COLORS AS CUMULATIVE PLOTS
  # Red → orange → blue gradient function (matching your cumulative plots)
  create_watershed_colors <- function(mgmt_rivers_ordered) {
    n_watersheds <- length(mgmt_rivers_ordered)
    
    # Create gradient from red → orange → blue
    color_gradient <- colorRampPalette(c(
      "#8B0000", "#CC0000", "#FF0000", "#FF4500", "#FF8C00", "#FFA500", "#FFB347",
      "#87CEEB", "#4682B4", "#1E90FF", "#0000FF", "#000080"
    ))(n_watersheds)
    
    # Create named vector
    names(color_gradient) <- mgmt_rivers_ordered
    return(color_gradient)
  }
  
  watershed_colors <- create_watershed_colors(final_order)
  
  # Create contemporary boxplot
  year_range <- paste0(min(boxplot_data$year), "-", max(boxplot_data$year))
  
  boxplot <- ggplot(boxplot_data, aes(x = mgmt_river, y = closure_percentage, fill = mgmt_river)) +
    # Modern boxplot styling
    geom_boxplot(
      alpha = 0.8,
      outlier.size = 3,
      outlier.shape = 21,
      outlier.fill = "white",
      outlier.stroke = 1.5,
      linewidth = 0.6,
      color = "gray40",
      show.legend = FALSE  # Remove legend since colors encode position
    ) +
    
    # Use exact same colors as cumulative plots
    scale_fill_manual(values = watershed_colors) +
    
    coord_flip() +
    
    # Clean, modern axis formatting
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      limits = c(0, max(boxplot_data$closure_percentage, na.rm = TRUE) * 1.05),
      expand = expansion(mult = c(0.02, 0.05)),
      breaks = function(x) pretty(x, n = 6)
    ) +
    
    # Contemporary labels
    labs(
      title = "Front-End Closure Protection by Management Unit",
      subtitle = paste("Percentage of each unit's annual production within Q1 closure window •", year_range),
      x = NULL,  # Remove x-axis title for cleaner look
      y = "% of Unit's Annual Production in Q1 Closure",
      caption = "Ordered by watershed position: upstream (top) → downstream (bottom) • Colors match cumulative distribution analysis"
    ) +
    
    # Modern, clean theme
    theme_minimal(base_size = 12) +
    theme(
      # Modern typography
      plot.title = element_text(
        face = "bold", 
        size = 18, 
        hjust = 0.5, 
        color = "gray15",
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(
        size = 13, 
        hjust = 0.5, 
        color = "gray50",
        margin = margin(b = 20)
      ),
      plot.caption = element_text(
        size = 10, 
        hjust = 0.5, 
        color = "gray60",
        margin = margin(t = 15)
      ),
      
      # Clean axis styling
      axis.text.y = element_text(
        size = 11, 
        color = "gray30",
        hjust = 1
      ),
      axis.text.x = element_text(
        size = 11, 
        color = "gray30"
      ),
      axis.title.x = element_text(
        size = 12, 
        color = "gray20",
        face = "bold",
        margin = margin(t = 15)
      ),
      
      # Minimal grid
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(
        color = "gray92", 
        linewidth = 0.4
      ),
      
      # Clean background
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      
      # Modern spacing
      plot.margin = margin(20, 25, 20, 20, "mm"),
      
      # Remove panel border for ultra-clean look
      panel.border = element_blank(),
      axis.line.x = element_line(color = "gray60", linewidth = 0.5)
    )
  
  # Save with high quality
  ggsave(
    file.path(OUTPUT_DIR, "front_end_closure_protection_CONTEMPORARY.png"), 
    boxplot, 
    width = 14, height = 10, 
    dpi = 300, 
    bg = "white"
  )
  
  # Also save as PDF for publications
  ggsave(
    file.path(OUTPUT_DIR, "front_end_closure_protection_CONTEMPORARY.pdf"), 
    boxplot, 
    width = 14, height = 10, 
    dpi = 300, 
    bg = "white"
  )
  
  # Create summary statistics
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
    arrange(mgmt_river)
  
  # Save summary statistics
  write.csv(summary_stats, file.path(OUTPUT_DIR, "closure_protection_summary_CONTEMPORARY.csv"), row.names = FALSE)
  
  # Print summary with color verification
  message("\nContemporary boxplot created with:")
  message(paste("  • Management units:", length(unique(boxplot_data$mgmt_river))))
  message(paste("  • Color consistency: MATCHED to cumulative distribution plots"))
  message(paste("  • Ordering: Upstream (top) → Downstream (bottom)"))
  message(paste("  • Years:", year_range))
  message(paste("  • Saved as: front_end_closure_protection_CONTEMPORARY.png/.pdf"))
  message(paste("  • Summary stats: closure_protection_summary_CONTEMPORARY.csv"))
  
  # Verify color mapping (first 5 units)
  message("\nColor verification (first 5 units):")
  for (i in 1:min(5, length(final_order))) {
    unit <- final_order[i]
    color <- watershed_colors[unit]
    message(paste("  ", i, ".", unit, "→", color))
  }
  
  # Print summary statistics
  message("\nSummary Statistics (mean % in Q1 closure):")
  summary_display <- summary_stats %>%
    select(mgmt_river, mean_closure_pct, n_years) %>%
    arrange(desc(mean_closure_pct))
  
  print(summary_display)
  
  return(boxplot_data)
}

################################################################################
# MAIN EXECUTION
################################################################################

# Set analysis parameters
years <- c(2017, 2018, 2019, 2020, 2021)
watershed <- "Kusko"

message("=== Running Front-End Closure Boxplot Analysis ===")
message(paste("Years:", paste(years, collapse = ", ")))
message(paste("Watershed:", watershed))

# Check if required data file exists
required_file <- here("Analysis_Results/Management_River_Analysis/management_river_analysis_tidy.csv")
if (!file.exists(required_file)) {
  message("ERROR: Required data file not found!")
  message("Expected: ", required_file)
  message("Please run the DOY quartile analysis first to generate this file.")
} else {
  message("Found required data file. Proceeding with analysis...")
  
  # Run the analysis
  results <- run_closure_boxplot_analysis(years, watershed)
  
  if (!is.null(results)) {
    message("\n=== Analysis Complete ===")
    message("✓ Contemporary boxplot created")
    message("✓ Colors match cumulative distribution plots")
    message("✓ Proper watershed ordering applied")
    message("✓ Summary statistics exported")
    message("\nOutput location: Figures/Front_End_Closure_Boxplots/")
  } else {
    message("❌ Analysis failed. Check error messages above.")
  }
}

