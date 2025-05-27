# Enhanced Management Unit Analysis
# Adds proportional contribution within each quartile (0-1 scaling)

# run_enhanced_management_analysis.R
# Complete script to run management unit analysis for Kuskokwim watershed
# Clear environment to avoid conflicts
rm(list = ls())

# Load required libraries
library(sf)
library(dplyr)
library(here)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(readr)

cat("=== Enhanced Kuskokwim Management Unit Analysis ===\n")
cat("Loading functions...\n")

# Source all required files in correct order
source(here("code/utils/spatial_utils.R"))       # Modified version with management units
source(here("code/assignment.R"))                # Existing assignment functions
source(here("code/export_management_analysis.R")) # New management export function

# Define analysis parameters
years_to_analyze <- c("2017","2018", "2019", "2020", "2021")
watershed <- "Kusko"

cat("Analysis parameters:\n")
cat("- Watershed:", watershed, "\n")
cat("- Years:", paste(years_to_analyze, collapse = ", "), "\n")
cat("- Output directory: Analysis_Results/Management_Units/\n\n")

# Run the complete management unit analysis
cat("Starting management unit analysis...\n")

results <- export_management_analysis_results(
  years = years_to_analyze,
  watersheds = watershed,
  output_dir = here("Analysis_Results/Management_Units")
)

# Create visualization if we have data
if (!is.null(results$mgmt_file) && file.exists(results$mgmt_file)) {
  cat("\nCreating visualization...\n")
  plot_path <- create_management_unit_plot(
    mgmt_data_file = results$mgmt_file,
    output_dir = here("Analysis_Results/Management_Units")
  )
  results$plot_file <- plot_path
}

# ENHANCED ANALYSIS: Add within-quartile proportional contributions
if (!is.null(results$mgmt_file) && file.exists(results$mgmt_file)) {
  cat("\n=== Calculating Within-Quartile Proportional Contributions ===\n")
  
  # Read the management unit data
  mgmt_data <- read.csv(results$mgmt_file)
  
  # Calculate proportional contributions within each year-quartile combination
  enhanced_mgmt_data <- mgmt_data %>%
    group_by(year, quartile) %>%
    mutate(
      # Total managed production within this year-quartile
      total_managed_in_quartile = sum(percent_of_total_run, na.rm = TRUE),
      
      # Proportional contribution within this quartile (0-1 scale)
      proportion_within_quartile = percent_of_total_run / total_managed_in_quartile,
      
      # Percentage contribution within this quartile (0-100 scale)
      percent_within_quartile = proportion_within_quartile * 100,
      
      # Number of management units active in this quartile
      mgmt_units_in_quartile = n()
    ) %>%
    ungroup() %>%
    # Handle any NaN values (if total_managed_in_quartile is 0)
    mutate(
      proportion_within_quartile = ifelse(is.nan(proportion_within_quartile), 0, proportion_within_quartile),
      percent_within_quartile = ifelse(is.nan(percent_within_quartile), 0, percent_within_quartile)
    )
  
  # Save the enhanced data
  enhanced_file_path <- file.path(here("Analysis_Results/Management_Units"), "enhanced_management_unit_analysis.csv")
  write.csv(enhanced_mgmt_data, enhanced_file_path, row.names = FALSE)
  cat("Enhanced analysis saved to:", enhanced_file_path, "\n")
  
  # Create summary showing both metrics
  cat("\n=== Summary with Both Metrics ===\n")
  summary_both_metrics <- enhanced_mgmt_data %>%
    group_by(mgmt_river) %>%
    summarise(
      # Original metric (% of total watershed run)
      total_contribution_to_watershed = sum(percent_of_total_run),
      avg_percent_of_watershed = mean(percent_of_total_run),
      
      # New metric (average proportional contribution within quartiles)
      avg_proportion_within_quartiles = mean(proportion_within_quartile),
      avg_percent_within_quartiles = mean(percent_within_quartile),
      max_proportion_within_quartile = max(proportion_within_quartile),
      
      # General stats
      years_present = n_distinct(year),
      quartiles_active = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(avg_proportion_within_quartiles))
  
  print(summary_both_metrics)
  
  # Show detailed comparison for a few examples
  cat("\n=== Detailed Example: Comparison of Metrics ===\n")
  example_comparison <- enhanced_mgmt_data %>%
    select(year, quartile, mgmt_river, 
           percent_of_total_run, percent_within_quartile, 
           total_managed_in_quartile, mgmt_units_in_quartile) %>%
    arrange(year, quartile, desc(percent_within_quartile)) %>%
    head(15)
  
  print(example_comparison)
  
  # Create visualizations comparing both metrics
  cat("\n=== Creating Enhanced Visualizations ===\n")
  
  # 1. Comparison plot: % of watershed vs % within quartile
  p1 <- ggplot(enhanced_mgmt_data, aes(x = percent_of_total_run, y = percent_within_quartile)) +
    geom_point(aes(color = mgmt_river), alpha = 0.7, size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    scale_x_continuous(labels = function(x) paste0(round(x, 1), "%")) +
    scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
    labs(
      title = "Management Unit Contributions: Two Different Metrics",
      subtitle = "Watershed % vs Within-Quartile % • Dashed line shows where metrics are equal",
      x = "Percent of Total Watershed Run",
      y = "Percent of Managed Rivers within Quartile",
      color = "Management River"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 2. Time series showing within-quartile proportions
  p2 <- enhanced_mgmt_data %>%
    mutate(
      quartile_decimal = case_when(
        quartile == "Q1" ~ 0.125,
        quartile == "Q2" ~ 0.375,
        quartile == "Q3" ~ 0.625,
        quartile == "Q4" ~ 0.875,
        TRUE ~ 0.5
      ),
      year_decimal = year + quartile_decimal
    ) %>%
    ggplot(aes(x = year_decimal, y = proportion_within_quartile, color = mgmt_river)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    scale_x_continuous(
      breaks = seq(min(enhanced_mgmt_data$year), max(enhanced_mgmt_data$year), by = 1),
      labels = scales::number_format(accuracy = 1)
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      title = "Management Unit Proportional Contributions Within Each Quartile",
      subtitle = "Shows relative dominance among managed rivers within each time period",
      x = "Year (with quarterly increments)",
      y = "Proportion of Managed Rivers (0-1 scale)",
      color = "Management River"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # 3. Faceted comparison by management river
  p3 <- enhanced_mgmt_data %>%
    select(year, quartile, mgmt_river, percent_of_total_run, percent_within_quartile) %>%
    pivot_longer(cols = c(percent_of_total_run, percent_within_quartile), 
                 names_to = "metric_type", values_to = "percentage") %>%
    mutate(
      metric_label = case_when(
        metric_type == "percent_of_total_run" ~ "% of Total Watershed",
        metric_type == "percent_within_quartile" ~ "% Within Quartile",
        TRUE ~ metric_type
      ),
      time_period = paste0(year, "-", quartile)
    ) %>%
    ggplot(aes(x = time_period, y = percentage, fill = metric_label)) +
    geom_col(position = "dodge", alpha = 0.8) +
    facet_wrap(~ mgmt_river, scales = "free_y", ncol = 3) +
    scale_fill_brewer(palette = "Set2", name = "Metric") +
    scale_y_continuous(labels = function(x) paste0(round(x, 1), "%")) +
    labs(
      title = "Comparison of Metrics by Management River",
      subtitle = "Blue: % of total watershed run • Orange: % of managed rivers within quartile",
      x = "Time Period",
      y = "Percentage"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  # Save plots
  plot_dir <- here("Analysis_Results/Management_Units")
  
  ggsave(file.path(plot_dir, "metric_comparison_scatter.png"), p1, width = 10, height = 8, dpi = 300)
  ggsave(file.path(plot_dir, "within_quartile_timeseries.png"), p2, width = 12, height = 8, dpi = 300)
  ggsave(file.path(plot_dir, "metrics_comparison_faceted.png"), p3, width = 14, height = 10, dpi = 300)
  
  # Show quartile-by-quartile breakdown
  cat("\n=== Within-Quartile Breakdown by Time Period ===\n")
  quartile_breakdown <- enhanced_mgmt_data %>%
    group_by(year, quartile) %>%
    arrange(year, quartile, desc(percent_within_quartile)) %>%
    summarise(
      total_managed_percent = first(total_managed_in_quartile),
      dominant_river = first(mgmt_river),
      dominant_river_within_quartile_percent = first(percent_within_quartile),
      num_active_rivers = first(mgmt_units_in_quartile),
      .groups = "drop"
    )
  
  print(quartile_breakdown)
  
  # Update results list
  results$enhanced_file <- enhanced_file_path
  results$enhanced_plots <- list(
    scatter = file.path(plot_dir, "metric_comparison_scatter.png"),
    timeseries = file.path(plot_dir, "within_quartile_timeseries.png"),
    faceted = file.path(plot_dir, "metrics_comparison_faceted.png")
  )
}

# Print final results summary
cat("\n=== Enhanced Analysis Complete! ===\n")
cat("Files created:\n")
if (!is.null(results$mgmt_file)) {
  cat("- Management unit data:", results$mgmt_file, "\n")
}
if (!is.null(results$enhanced_file)) {
  cat("- Enhanced management data:", results$enhanced_file, "\n")
}
if (!is.null(results$huc_file)) {
  cat("- HUC data:", results$huc_file, "\n")
}
if (!is.null(results$trib_file)) {
  cat("- Tributary data:", results$trib_file, "\n")
}
if (!is.null(results$summary_file)) {
  cat("- Summary data:", results$summary_file, "\n")
}
if (!is.null(results$plot_file)) {
  cat("- Original visualization:", results$plot_file, "\n")
}
if (!is.null(results$enhanced_plots)) {
  cat("- Enhanced visualizations:\n")
  for (plot_name in names(results$enhanced_plots)) {
    cat("  *", plot_name, ":", results$enhanced_plots[[plot_name]], "\n")
  }
}

cat("\n=== Key Metrics Explained ===\n")
cat("• 'percent_of_total_run': % of ENTIRE watershed run each management unit represents\n")
cat("• 'percent_within_quartile': % of MANAGED rivers each unit represents within that quartile\n")
cat("• 'proportion_within_quartile': Same as above but 0-1 scale instead of 0-100\n")
cat("\nThe within-quartile metrics show relative dominance among managed rivers in each time period.\n")