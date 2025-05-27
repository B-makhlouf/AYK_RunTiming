# run_management_analysis.R
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

cat("=== Kuskokwim Management Unit Analysis ===\n")
cat("Loading functions...\n")

# Source all required files in correct order
source(here("code/utils/spatial_utils.R"))       # Modified version with management units
source(here("code/assignment.R"))                # Existing assignment functions
source(here("code/export_management_analysis.R")) # New management export function

# Define analysis parameters
years_to_analyze <- c("2017", "2019", "2020", "2021")
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

# Print results summary
cat("\n=== Analysis Complete! ===\n")
cat("Files created:\n")
if (!is.null(results$mgmt_file)) {
  cat("- Management unit data:", results$mgmt_file, "\n")
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
  cat("- Visualization:", results$plot_file, "\n")
}

# Preview the main results
if (!is.null(results$mgmt_file) && file.exists(results$mgmt_file)) {
  cat("\n=== Preview of Management Unit Results ===\n")
  mgmt_data <- read.csv(results$mgmt_file)
  
  # Show summary by management unit
  cat("Summary by Management Unit (% of total watershed run):\n")
  summary_by_unit <- mgmt_data %>%
    group_by(mgmt_river) %>%
    summarise(
      total_contribution = sum(percent_of_total_run),
      avg_per_quartile = mean(percent_of_total_run),
      max_quartile = max(percent_of_total_run),
      years_present = n_distinct(year),
      .groups = "drop"
    ) %>%
    arrange(desc(total_contribution))
  
  print(summary_by_unit)
  
  # Show example of detailed data
  cat("\nExample detailed data (first 10 rows):\n")
  example_data <- mgmt_data %>%
    select(year, quartile, mgmt_river, percent_of_total_run, total_production, edge_count) %>%
    arrange(desc(percent_of_total_run)) %>%
    head(10)
  
  print(example_data)
  
  # Show overall managed vs unmanaged breakdown
  cat("\nManaged vs Unmanaged Summary by Year:\n")
  yearly_summary <- mgmt_data %>%
    group_by(year) %>%
    summarise(
      total_managed_percent = sum(percent_of_total_run),
      unmanaged_percent = 100 - total_managed_percent,
      .groups = "drop"
    )
  
  print(yearly_summary)
}

cat("\n=== Analysis completed successfully! ===\n")
cat("Key metric: 'percent_of_total_run' shows what % of the ENTIRE watershed run each management unit represents\n")

